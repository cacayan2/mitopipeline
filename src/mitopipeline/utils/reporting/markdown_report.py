"""
markdown_report.py

Renders one sample's real statistics (#356), biological validation
verdict (#414), and circular genome map (#349-352) as a per-sample
Markdown report (#358/#359), including a live methods paragraph
(methods.py) describing the actual tools/settings/rescue history used.
"""

# Imports
import os
import re
from datetime import datetime, timezone
from pathlib import Path
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
    """Renders the biological validation gate's own verdict (#414)."""
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

def render_archival_section(report_data: SampleReportData) -> list[str]:
    """Renders this sample's real NCBI archival preparation status (#367-376)."""
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
        *render_fastp_section(report_data),
        *render_assembly_section(report_data),
        *render_annotation_section(report_data),
        *render_blast_section(report_data),
        *render_validation_section(report_data),
        *render_visualization_section(report_data, report_path),
        *render_archival_section(report_data),
        *render_methods_section(report_data),
    ]
    return "\n".join(lines)

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
