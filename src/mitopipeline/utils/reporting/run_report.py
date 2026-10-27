"""
run_report.py

Renders the run-level Markdown report (#360/#441): aggregates every
sample's own biological validation verdict (PASS/REVIEW/FAIL, #414's
own typed outcome -- not a generic tool-exit-code "success/failure"
binary), plus the batch-level phylogeny outcome (#344-348's own
job-level design).
"""

# Imports
from pathlib import Path
from mitopipeline.utils.reporting.markdown_report import file_link, format_integer, markdown_table, relative_path
from mitopipeline.utils.reporting.models import RunReportData
from mitopipeline.utils.reporting.typst_style import render_title_block, render_typst_preamble

def render_run_overview_section(run_data: RunReportData) -> list[str]:
    """Renders run-level metadata and verdict totals."""
    counts = run_data.verdict_counts
    rows = [
        ["Job ID", run_data.job_id],
        ["Generated", run_data.generated_at],
        ["Total samples", format_integer(len(run_data.samples))],
        ["PASS", format_integer(counts.get("PASS", 0))],
        ["REVIEW", format_integer(counts.get("REVIEW", 0))],
        ["FAIL", format_integer(counts.get("FAIL", 0))],
        ["Not assessed", format_integer(counts.get("not assessed", 0))],
    ]
    return ["## Run overview", "", markdown_table(["Field", "Value"], rows), ""]

def render_sample_status_section(run_data: RunReportData, report_path: Path) -> list[str]:
    """Renders every sample's own verdict and a link to its report."""
    rows = []
    for sample in run_data.samples:
        report_link = file_link("sample report", report_path.parent / f"{sample.sample_id}.report.md", report_path)
        rows.append([
            sample.sample_id, sample.overall_verdict or "Not assessed",
            "Yes" if sample.rescued else "No", report_link or "Not available",
        ])
    return ["## Sample status", "", markdown_table(["Sample", "Verdict", "Rescued", "Report"], rows), ""]

def render_archival_summary_section(run_data: RunReportData) -> list[str]:
    """Renders the run-level NCBI archival preparation summary (#367-376)."""
    lines = ["## NCBI archival summary", ""]
    samples_with_archival = [s for s in run_data.samples if s.archival is not None]
    if not samples_with_archival:
        return [*lines, "_Archival preparation was not available for this run._", ""]

    rows = [
        ["Samples ready for submission review", f"{run_data.archival_ready_count} / {len(run_data.samples)}"],
    ]
    per_sample_rows = [
        [s.sample_id, s.archival.sra_review_status or "Not available",
         s.archival.biosample_review_status or "Not available", s.archival.genbank_review_status or "Not available",
         "Yes" if s.archival.ready_for_review else "No"]
        for s in samples_with_archival
    ]
    return [
        *lines, markdown_table(["Metric", "Value"], rows), "",
        "### Per-sample archival status", "",
        markdown_table(["Sample", "SRA", "BioSample", "GenBank", "Ready"], per_sample_rows), "",
        "This pipeline prepares submission-ready files only; it never submits to NCBI directly.", "",
    ]

def render_phylogeny_section(run_data: RunReportData) -> list[str]:
    """Renders the batch-level phylogeny outcome."""
    lines = ["## Phylogenetic validation", ""]
    phylogeny = run_data.phylogeny
    if phylogeny is None:
        return [*lines, "_Phylogenetic validation was not available for this run._", ""]
    if not phylogeny.success:
        return [*lines, f"_Phylogeny could not be built: {phylogeny.reason}_", ""]
    rows = [["Tree tips", format_integer(phylogeny.tip_count)]]
    if phylogeny.model is not None:
        rows.append(["Substitution model", phylogeny.model])
    return [*lines, markdown_table(["Metric", "Value"], rows), ""]

def render_run_report(run_data: RunReportData, report_path: str | Path) -> str:
    """
    Renders the complete run-level Markdown report.

    Args:
        run_data (RunReportData): The assembled run-level data.
        report_path (str | Path): Where this report will be written -- used to compute relative links.

    Returns:
        str: The complete Markdown report text.
    """
    report_path = Path(report_path)
    lines = [
        render_typst_preamble(),
        render_title_block(
            title = f"MitoPipeline Run Report: {run_data.job_id}",
            subtitle = f"{len(run_data.samples)} sample(s)", prepared_for = f"Job ID: {run_data.job_id}",
        ),
        *render_run_overview_section(run_data),
        *render_sample_status_section(run_data, report_path),
        *render_phylogeny_section(run_data),
        *render_archival_summary_section(run_data),
    ]
    return "\n".join(lines)

def write_run_report(run_data: RunReportData, output_path: str | Path) -> Path:
    """
    Writes the run-level Markdown report atomically.

    Args:
        run_data (RunReportData): The assembled run-level data.
        output_path (str | Path): Where to write the report.

    Returns:
        Path: The written report's path.
    """
    output_path = Path(output_path).resolve()
    output_path.parent.mkdir(parents = True, exist_ok = True)
    text = render_run_report(run_data, output_path)
    temporary_path = output_path.with_suffix(output_path.suffix + ".tmp")
    temporary_path.write_text(text, encoding = "utf-8")
    temporary_path.replace(output_path)
    return output_path
