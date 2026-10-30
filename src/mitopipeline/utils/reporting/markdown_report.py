"""
markdown_report.py

Renders one sample's real statistics, biological validation
verdict, and circular genome map as a per-sample
Markdown report, including a live methods paragraph
(methods.py) describing the actual tools/settings/rescue history used.
"""

# Imports
import os
import re
from datetime import datetime, timezone
from pathlib import Path
from mitopipeline.utils.archival.field_hints import RERUN_NOTE, field_hint
from mitopipeline.utils.reporting.compute_metrics import MachineInfo
from mitopipeline.utils.reporting.methods import render_methods_text
from mitopipeline.utils.reporting.models import SampleReportData
from mitopipeline.utils.reporting.typst_style import (
    STATUS_BADGE_STYLES,
    render_callout_box_close,
    render_callout_box_open,
    render_figure,
    render_title_block,
    render_typst_preamble,
)

# Every CommonMark/GFM ASCII punctuation character that can change
# meaning if left unescaped in a table cell or inline value -- a
# literal "*" in a real value (a note, an accession, a free-text
# field) is otherwise read as emphasis markup and silently dropped
# from the rendered output (confirmed directly: "*bold*" rendered as
# italic "bold", the asterisks themselves vanishing). Backslash must
# be escaped first, or escaping any of the others would double-escape
# through it.
_MARKDOWN_SPECIAL_CHARS = "\\`*_{}[]()#+!.$@"

def escape_markdown_text(value: object) -> str:
    """
    Escapes a real data value for safe literal use anywhere inside
    generated Markdown prose or a table cell -- every CommonMark
    special character is backslash-escaped so pandoc renders it as
    the literal input character, never as formatting.

    Args:
        value (object): The value to escape.

    Returns:
        str: The escaped text.
    """
    text = str(value)
    for char in _MARKDOWN_SPECIAL_CHARS:
        text = text.replace(char, "\\" + char)
    return text

def escape_markdown_cell(value: object) -> str:
    """Escapes a value for safe use inside a Markdown table cell."""
    text = str(value).replace("\n", " ").strip()
    text = escape_markdown_text(text)
    return text.replace("\\|", "|").replace("|", "\\|")

def wrap_cell_text(value: object, width: int = 22, extra_separators: tuple[str, ...] = ()) -> str:
    """
    Breaks a long cell value into multiple short lines so it wraps
    inside a table column instead of overflowing into the next one --
    Typst's own pandoc-generated tables size each column to its
    widest unbroken line with no automatic wrapping, confirmed
    directly (a 5-column table with "MANUAL_REVIEW_REQUIRED" cells
    overlapped the adjacent column's text). Breaks at the nearest
    space/comma/semicolon at or before width, joined by a real raw
    Typst line break (`<br>` does not survive pandoc's gfm reader,
    confirmed directly -- it is silently dropped, concatenating the
    surrounding words with no break and no space).

    Args:
        value (object): The value to wrap.
        width (int): The target maximum characters per line. Defaults to 22.
        extra_separators (tuple[str, ...]): Additional break points beyond
            space/comma/semicolon -- e.g. "_" for an unspaced identifier
            like a sample ID, confirmed directly to otherwise overflow
            into the next column in an already-tight multi-column table.

    Returns:
        str: The value, already Markdown-escaped, with inline raw Typst line breaks inserted as needed.
    """
    text = escape_markdown_cell(value)
    if len(text) <= width:
        return text

    separators = (" ", ", ", "; ", *extra_separators)
    lines: list[str] = []
    remaining = text
    while len(remaining) > width:
        break_at = -1
        for separator in separators:
            idx = remaining.rfind(separator, 0, width + 1)
            if idx > break_at:
                break_at = idx + len(separator.rstrip())
        if break_at <= 0:
            break_at = width
        lines.append(remaining[:break_at].rstrip(", ;"))
        remaining = remaining[break_at:].lstrip()
    lines.append(remaining)
    return "`#linebreak()`{=typst}".join(lines)

def humanize_key(key: str) -> str:
    """Converts a snake_case field/key name into a human-readable label (e.g. "expected_assembly_status" -> "Expected assembly status")."""
    words = str(key).replace("_", " ").replace("-", " ").split()
    if not words:
        return str(key)
    return " ".join([words[0].capitalize(), *words[1:]])

def humanize_enum(value: object) -> str:
    """Converts a raw enum-like value (SCREAMING_SNAKE_CASE or snake_case) into a human-readable label (e.g. "MANUAL_REVIEW_REQUIRED" -> "Manual review required")."""
    text = str(value)
    if not text:
        return text
    words = text.replace("_", " ").replace("-", " ").split()
    return " ".join([words[0].capitalize(), *[w.lower() for w in words[1:]]]) if words else text

def humanize_list(raw: str, separator: str = ";") -> str:
    """
    Converts a delimiter-joined list of raw field names into a
    de-duplicated, humanized, comma-separated string (e.g.
    "collection_date;geo_loc_name" -> "Collection date, Geo loc
    name"). Order-preserving de-duplication, since the same field
    name can legitimately appear from two independent sources (a
    missing file and a missing metadata field sharing one name).

    Args:
        raw (str): The raw delimiter-joined list.
        separator (str): The delimiter. Defaults to ";".

    Returns:
        str: A humanized, comma-separated, de-duplicated string, or "" if raw is empty.
    """
    if not raw:
        return ""
    seen: dict[str, None] = {}
    for item in raw.split(separator):
        item = item.strip()
        if item:
            seen.setdefault(item, None)
    return ", ".join(humanize_key(item) for item in seen)

def status_badge(value: object) -> str:
    """
    Renders a short status value (PASS/REVIEW/FAIL/READY/
    MANUAL_REVIEW_REQUIRED/etc.) as a small colored badge -- an
    inline raw Typst span, confirmed to work inside a GFM table cell
    via pandoc's raw_attribute extension. Colors are a fixed,
    colorblind-safe set (distinct hue AND lightness per tier, so
    grayscale printing still distinguishes them) and the label itself
    always carries the meaning, color is never the only signal.

    Args:
        value (object): The raw status value.

    Returns:
        str: An inline raw Typst badge, or the escaped humanized text if the value isn't a recognized status.
    """
    key = str(value).strip().upper()
    style = STATUS_BADGE_STYLES.get(key)
    label = humanize_enum(value)
    if style is None:
        return escape_markdown_cell(label)
    fill, text_color = style
    return (
        f'`#box(fill: rgb("{fill}"), inset: (x: 4pt, y: 2pt), radius: 2pt, '
        f'outset: (y: 1pt))[#text(fill: rgb("{text_color}"), size: 9.5pt)[{label}]]`{{=typst}}'
    )

def markdown_table(headers: list[str], rows: list[list[object]]) -> str:
    """Renders a GitHub-flavored Markdown table."""
    header = "| " + " | ".join(escape_markdown_cell(v) for v in headers) + " |"
    separator = "| " + " | ".join("---" for _ in headers) + " |"
    body = ["| " + " | ".join(str(v) for v in row) + " |" for row in rows]
    return "\n".join([header, separator, *body])

def format_integer(value: int | float) -> str:
    """Formats an integer-like value with thousands separators."""
    return f"{value:,.0f}"

def format_bp(value: int | float) -> str:
    """Formats a base-pair length with thousands separators and a space before the unit (e.g. "16,596 bp")."""
    return f"{format_integer(value)} bp"

def format_percent(value: float, fraction: bool = False) -> str:
    """Formats a percentage value; fraction=True scales a 0-1 ratio to a percentage first."""
    if fraction: value *= 100.0
    return f"{value:.2f}%"

def format_duration(seconds: float) -> str:
    """Formats a real runtime duration as "Hh Mm Ss", dropping any leading all-zero unit (e.g. "45.2s", "3m 2.1s", "1h 4m 0.0s")."""
    seconds = float(seconds)
    hours, remainder = divmod(seconds, 3600)
    minutes, remainder_seconds = divmod(remainder, 60)
    if hours >= 1:
        return f"{int(hours)}h {int(minutes)}m {remainder_seconds:.1f}s"
    if minutes >= 1:
        return f"{int(minutes)}m {remainder_seconds:.1f}s"
    return f"{remainder_seconds:.1f}s"

def format_memory_mb(megabytes: float) -> str:
    """Formats a real peak-memory value in MB, switching to GB above 1024 MB (e.g. "512.3 MB", "2.10 GB")."""
    megabytes = float(megabytes)
    if megabytes >= 1024:
        return f"{megabytes / 1024:.2f} GB"
    return f"{megabytes:.1f} MB"

def format_number_in_text(text: str) -> str:
    """
    Reformats bare integers immediately followed by "bp" inside a
    free-text string (a validation check's own "reason" sentence) to
    carry thousands separators and a space before the unit (e.g.
    "16596bp is within..." -> "16,596 bp is within...") -- the same
    convention every table value already uses, applied to prose that
    quotes a raw tool-computed number rather than a typed field.

    Args:
        text (str): The free-text string to reformat.

    Returns:
        str: The reformatted string.
    """
    def _replace(match: re.Match) -> str:
        return f"{int(match.group(1)):,} bp"
    text = re.sub(r"(\d+)bp\b", _replace, text)

    # The length-bounds gate's own "[min, max]bp" interval phrasing
    # (utils/validation/gates.py) puts the unit after the closing
    # bracket, not immediately after either number, so the pattern
    # above never reaches min_length/max_length on its own. This
    # always runs on already-escaped text (every caller applies
    # escape_markdown_cell() first), so the literal brackets arrive
    # backslash-escaped (`\[`/`\]`) -- confirmed directly via a real
    # compiled report still showing unformatted bounds despite this
    # function running, before the optional `\\?` was added here.
    def _replace_interval(match: re.Match) -> str:
        return f"{match.group(1)}{int(match.group(2)):,}, {int(match.group(3)):,}{match.group(4)} bp"
    return re.sub(r"(\\?\[)(\d+), (\d+)(\\?\])bp", _replace_interval, text)

def pluralize(count: int, singular: str, plural: str | None = None) -> str:
    """Formats a count with the correct singular/plural noun (e.g. pluralize(1, "sample") -> "1 sample")."""
    word = singular if count == 1 else (plural or f"{singular}s")
    return f"{format_integer(count)} {word}"

def relative_path(target_path: str | Path, report_path: str | Path) -> str:
    """Returns a Markdown-safe path to target_path, relative to report_path's own directory."""
    target_path, report_path = Path(target_path).resolve(), Path(report_path).resolve()
    return Path(os.path.relpath(target_path, start = report_path.parent)).as_posix()

def file_link(label: str, target_path: str | Path | None, report_path: str | Path) -> str | None:
    """Builds a Markdown link to target_path, or None when it doesn't exist."""
    if target_path is None or not Path(target_path).exists():
        return None
    return f"[{label}]({relative_path(target_path, report_path)})"

def basename_code(path: str) -> str:
    """Renders a file path as its basename, in monospace -- never the full path, which would otherwise expose the real execution host's absolute home-directory layout."""
    return f"`{os.path.basename(str(path))}`"

def format_scientific_name(metadata: dict[str, str]) -> str:
    """Formats genus/species/subspecies from manifest metadata into a scientific name."""
    genus, species = metadata.get("genus", "").strip(), metadata.get("species", "").strip()
    subspecies = metadata.get("subspecies", "").strip()
    parts = [genus, species]
    if subspecies.lower() not in {"", "na", "n/a", "none", "null"}:
        parts.append(subspecies)
    name = " ".join(p for p in parts if p)
    return name or metadata.get("scientific_name", "Not provided")

# Metadata fields whose real value is a filesystem path -- rendered as
# a monospace basename (never the full absolute path).
_PATH_METADATA_FIELDS = {"r1", "r2"}

def render_metadata_section(report_data: SampleReportData) -> list[str]:
    """Renders the sample's manifest metadata as a table."""
    rows = []
    for key, value in report_data.metadata.items():
        if key in _PATH_METADATA_FIELDS:
            rows.append([humanize_key(key), basename_code(value)])
        else:
            rows.append([humanize_key(key), wrap_cell_text(value, width = 48)])
    return ["## Sample Metadata", "", markdown_table(["Field", "Value"], rows), ""]

def render_aggregated_data_section(report_data: SampleReportData, report_path: Path) -> list[str]:
    """
    Links directly to this sample's own consolidated files in the
    run's aggregated_data/ tree (assembly FASTA, annotation GFF/map,
    phylogeny Newick/figure, real NCBI archival materials) -- quick
    access without hunting across a dozen per-stage directories for
    the files that actually matter, and a single place to review
    before the real manual NCBI submission.
    """
    lines = ["## Aggregated Data", ""]
    aggregated_dir = report_data.aggregated_data_dir
    if aggregated_dir is None:
        return [*lines, "_Aggregated data was not available for this run._", ""]

    sample_id = report_data.sample_id
    links = [
        link for link in [
            file_link("Assembly FASTA", aggregated_dir / "assemblies" / f"{sample_id}.assembly.fasta", report_path),
            file_link("Annotation GFF", aggregated_dir / "annotations" / "gff" / f"{sample_id}.annotation.gff", report_path),
            file_link("Annotated genome map", aggregated_dir / "annotations" / f"{sample_id}.annotated_map.png", report_path),
            file_link("Phylogenetic tree (Newick)", aggregated_dir / "phylogenies" / f"{sample_id}.treefile.nwk", report_path),
            file_link("Tree figure", aggregated_dir / "phylogenies" / f"{sample_id}.tree_figure.png", report_path),
        ]
        if link is not None
    ]
    if links:
        lines += ["Quick links to this sample's own consolidated output files:", "", "- " + "\n- ".join(links), ""]
    else:
        lines += ["_No aggregated files were available for this sample._", ""]

    archival_links = [
        link for link in [
            file_link("SRA metadata table", aggregated_dir / "archival" / "sra_metadata.tsv", report_path),
            file_link("BioSample metadata table", aggregated_dir / "archival" / "biosample_metadata.tsv", report_path),
            file_link("GenBank assembly FASTA", aggregated_dir / "archival" / "genbank" / sample_id / f"{sample_id}.fasta", report_path),
            file_link("GenBank annotation GFF", aggregated_dir / "archival" / "genbank" / sample_id / f"{sample_id}.gff", report_path),
            file_link("Submission readiness report", aggregated_dir / "archival" / "archival_validation.tsv", report_path),
            file_link("Submission checklist", aggregated_dir / "archival" / "submission_checklist.md", report_path),
            file_link("MD5 checksums", aggregated_dir / "archival" / "md5_checksums.tsv", report_path),
        ]
        if link is not None
    ]
    if archival_links:
        lines += ["NCBI archival materials, staged for manual review before submission:", "", "- " + "\n- ".join(archival_links), ""]

    return lines

def render_fastp_section(report_data: SampleReportData) -> list[str]:
    """Renders read trimming/filtering statistics."""
    stats = report_data.fastp_stats
    lines = ["## Read Trimming and Filtering", ""]
    if stats is None:
        return [*lines, "_Read trimming statistics were not available._", ""]
    rows = [
        ["Reads before filtering", format_integer(stats["reads_in"])],
        ["Reads after filtering", format_integer(stats["reads_out"])],
        ["Read retention rate", format_percent(stats["read_retention_rate"], fraction = True)],
        ["Q20 rate (after)", format_percent(stats["q20_rate_out"], fraction = True)],
        ["Q30 rate (after)", format_percent(stats["q30_rate_out"], fraction = True)],
        ["GC content (after)", format_percent(stats["gc_content_out"], fraction = True)],
        ["Low-quality reads removed", format_integer(stats["low_quality_reads"])],
    ]
    return [*lines, markdown_table(["Metric", "Value"], rows), ""]

def render_read_qc_section(report_data: SampleReportData) -> list[str]:
    """
    Renders FastQC's own real quality control for every read-set this
    sample's assembly actually depended on: the raw sequencing reads,
    the trimmed reads, and -- the reads GetOrganelle itself recruited
    as mitochondrial-like and actually handed to SPAdes for assembly,
    not just the raw input the rest of the pipeline started from.
    """
    lines = ["## Read Quality Control", ""]
    readsets = [
        ("Raw Reads", report_data.raw_read_qc),
        ("Trimmed Reads", report_data.trimmed_read_qc),
        ("Recruited Reads (used for assembly)", report_data.recruited_read_qc),
    ]
    if all(qc is None for _, qc in readsets):
        return [*lines, "_FastQC quality control was not available._", ""]

    for label, qc in readsets:
        lines += [f"### {label}", ""]
        if qc is None:
            lines += ["_Not available._", ""]
            continue

        rows = []
        if qc.total_sequences is not None:
            rows.append(["Total sequences", format_integer(qc.total_sequences)])
        if qc.sequence_length:
            rows.append(["Sequence length", f"{qc.sequence_length} bp"])
        if qc.percent_gc is not None:
            rows.append(["GC content", format_percent(qc.percent_gc)])
        rows.append(["FastQC verdict", status_badge(qc.overall_verdict)])
        lines += [markdown_table(["Metric", "Value"], rows), ""]
        if qc.failed_modules:
            lines += [f"Failed: {escape_markdown_cell(', '.join(qc.failed_modules))}", ""]
        if qc.warned_modules:
            lines += [f"Warned: {escape_markdown_cell(', '.join(qc.warned_modules))}", ""]

        if qc.representative_image_png is not None and Path(qc.representative_image_png).exists():
            image_path = str(Path(qc.representative_image_png).resolve())
            caption = f"Per-base sequence quality for {report_data.sample_id}'s {label.lower()} (R1)."
            lines += [render_figure(image_path, caption), ""]

    return lines

def render_assembly_section(report_data: SampleReportData) -> list[str]:
    """Renders assembly statistics."""
    stats = report_data.assembly
    lines = ["## Assembly", ""]
    if stats is None:
        return [
            *lines,
            render_callout(
                "No assembly: contig selection did not find a usable circular or "
                "fragmented sequence across the primary attempt and any rescue "
                "attempts for this sample."
            ),
            "",
        ]
    rows = [
        ["Status", humanize_enum(stats.status)],
        ["Source attempt", "Primary" if stats.source_attempt == 0 else f"Rescue attempt {stats.source_attempt}"],
        ["Contig count", format_integer(stats.contig_count)],
        ["Total length", format_bp(stats.length_bp)],
        ["GC content", format_percent(stats.gc_content_percent)],
    ]
    if stats.recruited_read_count is not None:
        rows.append(["Recruited reads", format_integer(stats.recruited_read_count)])
    if stats.runtime_seconds is not None:
        rows.append(["Assembly runtime", f"{stats.runtime_seconds:.2f} seconds"])
    return [*lines, markdown_table(["Metric", "Value"], rows), ""]

def render_annotation_section(report_data: SampleReportData) -> list[str]:
    """Renders MITOS2 annotation statistics."""
    stats = report_data.annotation
    lines = ["## Annotation", ""]
    if stats is None:
        if report_data.assembly is None:
            return [*lines, "_Not run: no assembly._", ""]
        return [*lines, "_Annotation statistics were not available._", ""]
    rows = [
        ["Protein-coding genes", f"{stats.core_gene_count} / 13"],
    ]
    if stats.trna_count is not None:
        rows.append(["tRNAs", f"{stats.trna_count} / 22"])
    if stats.rrna_count is not None:
        rows.append(["rRNAs", f"{stats.rrna_count} / 2"])
    section = [*lines, markdown_table(["Metric", "Value"], rows), ""]
    if stats.missing_core_genes:
        section += [f"Missing core genes: {humanize_list(';'.join(stats.missing_core_genes))}", ""]
    return section

def render_blast_section(report_data: SampleReportData) -> list[str]:
    """Renders the top BLAST identification match."""
    match = report_data.blast
    lines = ["## Closest BLAST Match", ""]
    if match is None:
        if report_data.assembly is None:
            return [*lines, "_Not run: no assembly._", ""]
        return [*lines, "_BLAST identification results were not available._", ""]
    scientific_name = match.top_scientific_name
    name_display = scientific_name if scientific_name and scientific_name.upper() != "N/A" else "Taxon unresolved"
    rows = [
        ["Accession", basename_code(match.top_accession.strip("|").removeprefix("ref|").removeprefix("ref"))
         if match.top_accession.startswith("ref|") else escape_markdown_cell(match.top_accession)],
        ["Scientific name", escape_markdown_cell(name_display)],
        ["Identity", format_percent(match.top_identity_percent)],
        ["Query coverage", format_percent(match.top_query_coverage_percent)],
        ["E-value", f"{match.top_evalue:.3g}"],
    ]
    if match.decision_status is not None:
        rows.append(["Decision status", status_badge(match.decision_status)])
    return [*lines, markdown_table(["Metric", "Value"], rows), ""]

def render_callout(text: str) -> str:
    """A single-line gray callout note, for a short standalone message (e.g. "no assembly was produced")."""
    return f'`#block(fill: rgb("#f2f2f2"), inset: 8pt, radius: 4pt, width: 100%)[{escape_markdown_text(text)}]`{{=typst}}'

def render_validation_section(report_data: SampleReportData) -> list[str]:
    """Renders the biological validation gate's own verdict."""
    lines = ["## Biological Validation", ""]
    if report_data.overall_verdict is None:
        if report_data.assembly is None:
            return [*lines, "_Not run: no assembly._", ""]
        return [*lines, "_Biological validation was not available._", ""]
    lines.append(f"**Plausibility for Complete Genome:** {status_badge(report_data.overall_verdict)}")
    lines.append("")
    rows = [
        [humanize_key(c.get("name", "")), status_badge(c.get("verdict", "")), format_number_in_text(escape_markdown_cell(c.get("reason", "")))]
        for c in report_data.validation_checks
    ]
    if rows:
        lines += [markdown_table(["Check", "Verdict", "Reason"], rows), ""]
    return lines

def render_machine_info_rows(info: MachineInfo) -> list[list[str]]:
    """
    Builds the real machine-spec table rows shared by both the
    per-sample and run-level reports -- deliberately only
    non-identifying specs (no hostname/username/IP). A GPU row always
    appears when one was detected, carrying a fixed disclaimer that
    this pipeline doesn't use it (recorded for reproducibility and any
    future GPU-accelerated tooling, not because anything here runs on
    it today).

    Args:
        info (MachineInfo): The real collected machine specs.

    Returns:
        list[list[str]]: One [field, value] row per spec.
    """
    cpu_cores_value = "Not available"
    if info.cpu_physical_cores is not None or info.cpu_logical_threads is not None:
        physical = format_integer(info.cpu_physical_cores) if info.cpu_physical_cores is not None else "?"
        logical = format_integer(info.cpu_logical_threads) if info.cpu_logical_threads is not None else "?"
        cpu_cores_value = f"{physical} physical / {logical} logical (threads)"

    memory_spec_value = "Not available"
    if info.memory_type or info.memory_speed_mhz is not None:
        memory_spec_value = " ".join(filter(None, [
            info.memory_type, f"{info.memory_speed_mhz} MHz" if info.memory_speed_mhz is not None else None,
        ]))

    rows = [
        ["Operating system", f"{info.os_name} {info.os_release}"],
        ["Architecture", info.architecture],
        ["CPU model", info.cpu_model or "Not available"],
        ["CPU cores (physical / logical)", cpu_cores_value],
        ["CPU clock speed", f"{info.cpu_clock_mhz:.0f} MHz" if info.cpu_clock_mhz is not None else "Not available"],
        ["Total RAM", f"{info.total_memory_gb} GiB" if info.total_memory_gb is not None else "Not available"],
        ["RAM speed / type", memory_spec_value],
        [
            "GPU",
            f"{info.gpu_model} (not utilized by any part of this pipeline; reported for reproducibility and any future GPU-accelerated tooling)"
            if info.gpu_model else "No GPU detected",
        ],
        ["Python version", info.python_version],
    ]
    return rows

def render_compute_resources_section(report_data: SampleReportData) -> list[str]:
    """Renders this sample's own real compute environment and per-stage compute cost, plus exactly where to find the underlying raw metrics/log files."""
    lines = ["## Compute Resources", ""]
    if report_data.machine_info is not None:
        lines += ["### Compute Environment", "", markdown_table(["Field", "Value"], render_machine_info_rows(report_data.machine_info)), ""]
    if not report_data.compute_summary:
        return [*lines, "_Compute metrics were not available for this sample._", ""]

    rows = [
        [
            basename_code(s.stage), s.tools or "—", format_integer(s.invocation_count),
            format_duration(s.total_runtime_seconds), format_duration(s.total_cpu_seconds),
            format_memory_mb(s.peak_memory_mb),
        ]
        for s in report_data.compute_summary
    ]
    lines += [markdown_table(["Stage", "Tool(s)", "Invocations", "Runtime", "CPU time", "Peak memory"], rows), ""]

    access_lines = []
    if report_data.log_dir is not None:
        access_lines.append(f"- **Log files:** `{report_data.log_dir}/<stage>/{report_data.sample_id}.log` (one per real stage this sample went through)")
    if report_data.metrics_dir is not None:
        access_lines.append(f"- **Raw compute metrics:** `{report_data.metrics_dir}/<stage>/{report_data.sample_id}/` (one real JSON record per tool invocation -- runtime, CPU, peak memory)")
    if access_lines:
        lines += ["**Accessing the underlying data:**", "", *access_lines, ""]
    return lines

def render_archival_section(report_data: SampleReportData) -> list[str]:
    """Renders this sample's real NCBI archival preparation status."""
    lines = ["## NCBI Archival Readiness", ""]
    archival = report_data.archival
    if archival is None:
        return [*lines, "_Archival preparation was not available._", ""]

    rows = [
        ["SRA", status_badge(archival.sra_review_status) if archival.sra_review_status else "Not available",
         wrap_cell_text(humanize_list(archival.sra_missing_fields) or "—", width = 24)],
        ["BioSample", status_badge(archival.biosample_review_status) if archival.biosample_review_status else "Not available",
         wrap_cell_text(humanize_list(archival.biosample_missing_fields) or "—", width = 24)],
        ["GenBank", status_badge(archival.genbank_review_status) if archival.genbank_review_status else "Not available",
         wrap_cell_text(humanize_list(archival.genbank_missing) or "—", width = 24)],
    ]
    lines += [markdown_table(["Record", "Review status", "Needs review"], rows), ""]
    ready_badge = status_badge("READY") if archival.ready_for_review else status_badge("MANUAL_REVIEW_REQUIRED")
    lines += [f"**Ready for submission review:** {ready_badge}", ""]
    lines += [
        "This pipeline prepares submission-ready files only; it never submits to NCBI directly. "
        "A human reviewer should resolve any fields listed above before uploading through the "
        "NCBI Submission Portal.", "",
    ]

    missing_fields = dict.fromkeys(f for raw in (
        archival.sra_missing_fields, archival.biosample_missing_fields, archival.genbank_missing,
    ) for f in raw.split(";") if f)
    if missing_fields:
        lines += [f"**How to complete:** {RERUN_NOTE}", ""]
        lines += [f"- `{field}`: {field_hint(field)}" for field in missing_fields]
        lines += [""]

    return lines

def render_visualization_section(report_data: SampleReportData, report_path: Path) -> list[str]:
    """Renders the circular mitochondrial genome map as a numbered Typst figure, kept with this heading."""
    lines = ["## Circular Mitochondrial Genome Map", ""]
    if report_data.visualization_png is None or not Path(report_data.visualization_png).exists():
        return [*lines, "_Circular genome visualization was not available._", ""]

    # An absolute path, not relative_path()'s report-relative form: pandoc's
    # typst PDF engine compiles through an intermediate file in its own
    # temp directory, so any path relative to the report's own location
    # resolves against the wrong directory there (confirmed via a real
    # pandoc-to-pdf compile). pdf_export.py pairs this with --root=/ so
    # Typst's sandbox permits resolving the real absolute location.
    image_path = str(Path(report_data.visualization_png).resolve())
    caption = f"Annotated mitochondrial genome map for {report_data.sample_id}"
    if report_data.assembly is not None:
        caption += f" ({format_bp(report_data.assembly.length_bp)}, GC {report_data.assembly.gc_content_percent:.1f}%)."
    else:
        caption += "."
    if report_data.assembly is not None and report_data.assembly.status != "circular":
        caption += (
            " The assembly did not close into a circle; the map is drawn circularly by the "
            "plotting tool and does not imply a closed genome."
        )
    lines += [render_figure(image_path, caption), ""]
    return lines

def render_sample_report(report_data: SampleReportData, report_path: str | Path) -> str:
    """
    Renders one sample's complete Markdown report.

    Args:
        report_data (SampleReportData): The assembled report data.
        report_path (str | Path): Where this report will be written -- used to compute relative links.

    Returns:
        str: The complete Markdown report text.
    """
    report_path = Path(report_path)
    species = format_scientific_name(report_data.metadata)
    generated_at = datetime.now(timezone.utc).strftime("%d %b %Y, %H:%M UTC")

    lines = [
        render_typst_preamble(sample_id = report_data.sample_id),
        render_title_block(
            title = "Mitochondrial Genome Report", species = species,
            sample_id = report_data.sample_id, generated_at = generated_at,
        ),
        *render_metadata_section(report_data),
        *render_aggregated_data_section(report_data, report_path),
        *render_fastp_section(report_data),
        *render_read_qc_section(report_data),
        *render_assembly_section(report_data),
        *render_annotation_section(report_data),
        *render_blast_section(report_data),
        *render_validation_section(report_data),
        *render_visualization_section(report_data, report_path),
        *render_phylogeny_section(report_data, report_path),
        *render_archival_section(report_data),
        *render_compute_resources_section(report_data),
        *render_methods_section(report_data),
    ]
    return "\n".join(lines)

def render_phylogeny_section(report_data: SampleReportData, report_path: Path) -> list[str]:
    """
    Renders this sample's own real phylogenetic tree -- its
    assembled genome placed among its own top-ranked BLAST reference
    matches (restoring v0.1.0's per-sample design: a sample's own
    tree never depends on how many other samples are in the same
    run).
    """
    lines = ["## Phylogenetic Context", ""]
    phylogeny = report_data.phylogeny
    if phylogeny is None:
        return [*lines, "_Phylogenetic placement was not available._", ""]
    if not phylogeny.success:
        return [*lines, render_callout(phylogeny.reason or "Phylogenetic placement could not be built for this sample."), ""]

    rows = []
    if phylogeny.model:
        rows.append(["Model", f"`{phylogeny.model}`"])
    if phylogeny.tip_count is not None:
        rows.append(["Tree tips", format_integer(phylogeny.tip_count)])
    if rows:
        lines += [markdown_table(["Metric", "Value"], rows), ""]

    if phylogeny.tree_figure_png is not None and Path(phylogeny.tree_figure_png).exists():
        image_path = str(Path(phylogeny.tree_figure_png).resolve())
        caption = (
            f"Maximum-likelihood mitochondrial phylogeny for {report_data.sample_id}, "
            "placed among its top BLAST reference matches."
        )
        lines += [render_figure(image_path, caption), ""]

    newick_link = file_link("Newick tree", phylogeny.newick_path, report_path)
    if newick_link is not None:
        lines += [f"**Download:** {newick_link}", ""]

    return lines

def render_methods_section(report_data: SampleReportData) -> list[str]:
    """
    Renders the "Methods" section as Markdown lines, with the
    paragraph itself set inside a gray callout box matching the
    project's own real document styling.

    Args:
        report_data (SampleReportData): The assembled report data.

    Returns:
        list[str]: The Markdown lines for this section.
    """
    return [
        "## Methods", "",
        render_callout_box_open(),
        render_methods_text(report_data),
        "",
        render_callout_box_close(),
        "",
    ]

def write_sample_report(report_data: SampleReportData, output_path: str | Path) -> Path:
    """
    Writes one sample's Markdown report atomically.

    Args:
        report_data (SampleReportData): The assembled report data.
        output_path (str | Path): Where to write the report.

    Returns:
        Path: The written report's path.
    """
    output_path = Path(output_path).resolve()
    output_path.parent.mkdir(parents = True, exist_ok = True)
    text = render_sample_report(report_data, output_path)
    temporary_path = output_path.with_suffix(output_path.suffix + ".tmp")
    temporary_path.write_text(text, encoding = "utf-8")
    temporary_path.replace(output_path)
    return output_path
