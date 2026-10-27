"""
markdown_report.py

Renders one sample's real statistics (#356), biological validation
verdict (#414), and circular genome map (#349-352) as a per-sample
Markdown report (#358/#359), including a live methods paragraph
(methods.py) describing the actual tools/settings/rescue history used.
"""

# Imports
import os
from datetime import datetime, timezone
from pathlib import Path
from mitopipeline.utils.reporting.methods import render_methods_section
from mitopipeline.utils.reporting.models import SampleReportData
from mitopipeline.utils.reporting.typst_style import render_title_block, render_typst_preamble

def escape_markdown_cell(value: object) -> str:
    """Escapes a value for safe use inside a Markdown table cell."""
    return str(value).replace("|", "\\|").replace("\n", " ").strip()

def markdown_table(headers: list[str], rows: list[list[object]]) -> str:
    """Renders a GitHub-flavored Markdown table."""
    header = "| " + " | ".join(escape_markdown_cell(v) for v in headers) + " |"
    separator = "| " + " | ".join("---" for _ in headers) + " |"
    body = ["| " + " | ".join(escape_markdown_cell(v) for v in row) + " |" for row in rows]
    return "\n".join([header, separator, *body])

def format_integer(value: int | float) -> str:
    """Formats an integer-like value with thousands separators."""
    return f"{value:,.0f}"

def format_percent(value: float, fraction: bool = False) -> str:
    """Formats a percentage value; fraction=True scales a 0-1 ratio to a percentage first."""
    if fraction: value *= 100.0
    return f"{value:.2f}%"

def relative_path(target_path: str | Path, report_path: str | Path) -> str:
    """Returns a Markdown-safe path to target_path, relative to report_path's own directory."""
    target_path, report_path = Path(target_path).resolve(), Path(report_path).resolve()
    return Path(os.path.relpath(target_path, start = report_path.parent)).as_posix()

def file_link(label: str, target_path: str | Path | None, report_path: str | Path) -> str | None:
    """Builds a Markdown link to target_path, or None when it doesn't exist."""
    if target_path is None or not Path(target_path).exists():
        return None
    return f"[{label}]({relative_path(target_path, report_path)})"

def format_scientific_name(metadata: dict[str, str]) -> str:
    """Formats genus/species/subspecies from manifest metadata into a scientific name."""
    genus, species = metadata.get("genus", "").strip(), metadata.get("species", "").strip()
    subspecies = metadata.get("subspecies", "").strip()
    parts = [genus, species]
    if subspecies.lower() not in {"", "na", "n/a", "none", "null"}:
        parts.append(subspecies)
    name = " ".join(p for p in parts if p)
    return name or metadata.get("scientific_name", "Not provided")

def render_metadata_section(report_data: SampleReportData) -> list[str]:
    """Renders the sample's manifest metadata as a table."""
    rows = [[k, v] for k, v in report_data.metadata.items()]
    return ["## Sample metadata", "", markdown_table(["Field", "Value"], rows), ""]

def render_fastp_section(report_data: SampleReportData) -> list[str]:
    """Renders read trimming/filtering statistics."""
    stats = report_data.fastp_stats
    lines = ["## Read trimming and filtering", ""]
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
        return [*lines, "_Assembly statistics were not available._", ""]
    rows = [
        ["Status", stats.status],
        ["Source attempt", "Primary" if stats.source_attempt == 0 else f"Rescue attempt {stats.source_attempt}"],
        ["Contig count", format_integer(stats.contig_count)],
        ["Total length", f"{format_integer(stats.length_bp)} bp"],
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
        return [*lines, "_Annotation statistics were not available._", ""]
    rows = [
        ["Total gene features", format_integer(stats.gene_count)],
        ["Core protein-coding genes found", f"{stats.core_gene_count} / 13"],
    ]
    section = [*lines, markdown_table(["Metric", "Value"], rows), ""]
    if stats.missing_core_genes:
        section += ["Missing core genes: " + ", ".join(f"`{g}`" for g in stats.missing_core_genes), ""]
    return section

def render_blast_section(report_data: SampleReportData) -> list[str]:
    """Renders the top BLAST identification match."""
    match = report_data.blast
    lines = ["## Closest BLAST match", ""]
    if match is None:
        return [*lines, "_BLAST identification results were not available._", ""]
    rows = [
        ["Accession", match.top_accession],
        ["Scientific name", match.top_scientific_name or "Not reported"],
        ["Identity", f"{match.top_identity_percent:.2f}%"],
        ["Query coverage", f"{match.top_query_coverage_percent:.2f}%"],
        ["E-value", f"{match.top_evalue:.3g}"],
    ]
    if match.decision_status is not None:
        rows.append(["Decision status", match.decision_status])
    return [*lines, markdown_table(["Metric", "Value"], rows), ""]

def render_validation_section(report_data: SampleReportData) -> list[str]:
    """Renders the biological validation gate's own verdict (#414)."""
    lines = ["## Biological validation", ""]
    if report_data.overall_verdict is None:
        return [*lines, "_Biological validation was not available._", ""]
    lines.append(f"**Overall verdict: {report_data.overall_verdict}**")
    lines.append("")
    rows = [[c.get("name"), c.get("verdict"), c.get("reason")] for c in report_data.validation_checks]
    if rows:
        lines += [markdown_table(["Check", "Verdict", "Reason"], rows), ""]
    return lines

def render_visualization_section(report_data: SampleReportData, report_path: Path) -> list[str]:
    """Renders the circular mitochondrial genome map."""
    lines = ["## Circular mitochondrial genome map", ""]
    link = file_link("PNG circular map", report_data.visualization_png, report_path)
    if link is None:
        return [*lines, "_Circular genome visualization was not available._", ""]
    lines += [f"- {link}", ""]
    image_path = relative_path(report_data.visualization_png, report_path)
    lines += ["![Circular mitochondrial genome map]" + f"({image_path})", ""]
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
    generated_at = datetime.now(timezone.utc).isoformat(timespec = "seconds")

    lines = [
        render_typst_preamble(),
        render_title_block(
            title = f"Mitochondrial Genome Report: {report_data.sample_id}",
            subtitle = species, prepared_for = f"Sample ID: {report_data.sample_id}",
        ),
        f"**Sample ID:** `{report_data.sample_id}`  ",
        f"**Species:** {species}  ",
        f"**Generated:** {generated_at}  ",
        "",
        *render_metadata_section(report_data),
        *render_fastp_section(report_data),
        *render_assembly_section(report_data),
        *render_annotation_section(report_data),
        *render_blast_section(report_data),
        *render_validation_section(report_data),
        *render_visualization_section(report_data, report_path),
        *render_methods_section(report_data),
    ]
    return "\n".join(lines)

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
