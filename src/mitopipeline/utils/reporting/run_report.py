"""
run_report.py

Renders the run-level Markdown report (#360/#441): aggregates every
sample's own biological validation verdict (PASS/REVIEW/FAIL, #414's
own typed outcome -- not a generic tool-exit-code "success/failure"
binary). Each sample's own real phylogenetic tree (restoring v0.1.0's
per-sample design) lives in that sample's own report, not here.
"""

# Imports
from pathlib import Path
from mitopipeline.utils.reporting.markdown_report import (
    basename_code,
    file_link,
    format_bp,
    format_duration,
    format_integer,
    format_memory_mb,
    format_percent,
    humanize_enum,
    markdown_table,
    pluralize,
    render_machine_info_rows,
    status_badge,
    wrap_cell_text,
)
from mitopipeline.utils.reporting.models import RunReportData
from mitopipeline.utils.reporting.typst_style import render_run_title_block, render_typst_preamble

def render_run_overview_section(run_data: RunReportData) -> list[str]:
    """Renders run-level metadata and verdict totals."""
    counts = run_data.verdict_counts
    rows = [
        ["Job ID", f"`{run_data.job_id}`"],
        ["Total samples", pluralize(len(run_data.samples), "sample")],
        ["Pass", status_badge("PASS") + f" {format_integer(counts.get('PASS', 0))}"],
        ["Review", status_badge("REVIEW") + f" {format_integer(counts.get('REVIEW', 0))}"],
        ["Fail", status_badge("FAIL") + f" {format_integer(counts.get('FAIL', 0))}"],
        ["Not assessed", status_badge("NOT ASSESSED") + f" {format_integer(counts.get('not assessed', 0))}"],
    ]
    return ["## Run Overview", "", markdown_table(["Field", "Value"], rows), ""]

def render_job_metadata_section(run_data: RunReportData) -> list[str]:
    """Renders this job's own real lifecycle metadata, compute environment, batch-wide per-stage resource cost, and exactly where to find the underlying logs/metrics -- all of it already being recorded on every real run, surfaced here rather than left unread."""
    lines = ["## Job Execution Metadata", ""]

    meta = run_data.job_metadata
    if meta is not None:
        # This report is itself one of the last rules in the same
        # Snakemake run it's describing, so current_metadata.json's
        # own terminal "completed"/"failed" entry -- written only
        # after the whole run finishes -- is never visible yet; its
        # status structurally always still reads "created" here. Shown
        # honestly as "Running" with a real elapsed-so-far duration
        # (collect_job_metadata() already substitutes this for the
        # stale runtime_seconds=0.0 placeholder), not a misleading
        # "Created" badge implying nothing happened.
        still_running = meta.status == "created"
        status_display = "Running" if still_running else (meta.status.upper() if meta.status else None)
        run_length_label = "Elapsed so far (job was still running when this report was generated)" if still_running else "Run length"
        rows = [
            ["Status", status_badge(status_display) if status_display else "Not available"],
            ["Started", meta.start_time or "Not available"],
            [run_length_label, format_duration(meta.runtime_seconds) if meta.runtime_seconds is not None else "Not available"],
            ["Samples processed", format_integer(meta.sample_count) if meta.sample_count is not None else "Not available"],
            ["Total attempts (all time)", format_integer(meta.attempt_count)],
            ["Cumulative runtime (all attempts)", format_duration(meta.cumulative_runtime_seconds)],
            ["Incomplete attempts", format_integer(meta.incomplete_attempts)],
        ]
        lines += [markdown_table(["Field", "Value"], rows), ""]
    else:
        lines += ["_Job lifecycle metadata was not available for this run._", ""]

    info = run_data.machine_info
    if info is not None:
        lines += ["### Compute Environment", "", markdown_table(["Field", "Value"], render_machine_info_rows(info)), ""]

    if run_data.compute_summary:
        rows = [
            [
                basename_code(s.stage), s.tools or "—", format_integer(s.invocation_count),
                format_duration(s.total_runtime_seconds), format_duration(s.total_cpu_seconds),
                format_memory_mb(s.peak_memory_mb),
            ]
            for s in run_data.compute_summary
        ]
        lines += [
            "### Per-Stage Compute Summary (every sample)", "",
            markdown_table(["Stage", "Tool(s)", "Invocations", "Total Runtime", "Total CPU Time", "Peak Memory"], rows), "",
        ]

    access_lines = []
    if run_data.job_dir is not None:
        access_lines.append(f"- **Job output root:** `{run_data.job_dir}`")
    if run_data.global_log_file is not None:
        access_lines.append(f"- **Shared global log:** `{run_data.global_log_file}`")
    if run_data.job_dir is not None:
        access_lines.append(f"- **Per-sample/per-stage logs:** `{run_data.job_dir}/logs/<stage>/<sample>.log`")
        access_lines.append(f"- **Raw compute metrics:** `{run_data.job_dir}/metrics/<stage>/<sample>/` (one real JSON record per tool invocation)")
        access_lines.append(f"- **Job status history:** `{run_data.job_dir}/current_metadata.json`, `{run_data.job_dir}/historical_metadata.json`")
    if access_lines:
        lines += ["**Accessing the underlying data:**", "", *access_lines, ""]

    return lines

def render_sample_results_table_section(run_data: RunReportData) -> list[str]:
    """Renders a compact per-sample results table: assembly status, length, contig count, GC, core gene count, BLAST identity, verdict."""
    lines = ["## Per-Sample Results", ""]
    rows = []
    for sample in run_data.samples:
        assembly, blast, annotation = sample.assembly, sample.blast, sample.annotation
        rows.append([
            # Plain wrapped text, not basename_code()'s monospace
            # backticks: this table already packs 8 columns into the
            # page width, and a Markdown code span can't carry the
            # embedded raw-Typst linebreak wrap_cell_text() needs to
            # break a long, space-free sample ID -- confirmed directly
            # that a long ID otherwise overflows into the Status
            # column with no gap at all between the two.
            wrap_cell_text(sample.sample_id, width = 14, extra_separators = ("_",)),
            humanize_enum(assembly.status) if assembly else "—",
            format_bp(assembly.length_bp) if assembly else "—",
            format_integer(assembly.contig_count) if assembly else "—",
            format_percent(assembly.gc_content_percent) if assembly else "—",
            f"{annotation.core_gene_count} / 13" if annotation else "—",
            format_percent(blast.top_identity_percent) if blast else "—",
            status_badge(sample.overall_verdict) if sample.overall_verdict else status_badge("NOT ASSESSED"),
        ])
    headers = ["Sample", "Status", "Length", "Contigs", "GC", "Core genes", "BLAST identity", "Verdict"]
    return [*lines, markdown_table(headers, rows), ""]

def render_sample_status_section(run_data: RunReportData, report_path: Path) -> list[str]:
    """Renders every sample's own verdict and a link to its real PDF report."""
    rows = []
    for sample in run_data.samples:
        # The real per-sample report is nested one directory deeper
        # than this file (report_path.parent/{sample}/{sample}.report.pdf),
        # not a sibling of it -- confirmed directly against the real
        # output layout (a prior version of this link pointed at the
        # wrong, non-nested path and always showed "Not available"
        # even though the PDF genuinely existed).
        pdf_path = report_path.parent / sample.sample_id / f"{sample.sample_id}.report.pdf"
        report_link = file_link(f"{sample.sample_id}.report.pdf", pdf_path, report_path)
        rows.append([
            basename_code(sample.sample_id),
            status_badge(sample.overall_verdict) if sample.overall_verdict else status_badge("NOT ASSESSED"),
            "Yes" if sample.rescued else "No", report_link or "Not available",
        ])
    return ["## Sample Status", "", markdown_table(["Sample", "Verdict", "Rescued", "Report"], rows), ""]

def render_archival_summary_section(run_data: RunReportData) -> list[str]:
    """Renders the run-level NCBI archival preparation summary (#367-376)."""
    lines = ["## NCBI Archival Summary", ""]
    samples_with_archival = [s for s in run_data.samples if s.archival is not None]
    if not samples_with_archival:
        return [*lines, "_Archival preparation was not available for this run._", ""]

    rows = [
        ["Samples ready for submission review", f"{run_data.archival_ready_count} / {len(run_data.samples)}"],
    ]
    per_sample_rows = [
        [basename_code(s.sample_id),
         status_badge(s.archival.sra_review_status) if s.archival.sra_review_status else "Not available",
         status_badge(s.archival.biosample_review_status) if s.archival.biosample_review_status else "Not available",
         status_badge(s.archival.genbank_review_status) if s.archival.genbank_review_status else "Not available",
         status_badge("READY") if s.archival.ready_for_review else status_badge("MANUAL_REVIEW_REQUIRED")]
        for s in samples_with_archival
    ]
    return [
        *lines, markdown_table(["Metric", "Value"], rows), "",
        "### Per-Sample Archival Status", "",
        markdown_table(["Sample", "SRA", "BioSample", "GenBank", "Ready"], per_sample_rows), "",
        "This pipeline prepares submission-ready files only; it never submits to NCBI directly.", "",
    ]

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
        render_run_title_block(
            title = "MitoPipeline Run Report", job_id = run_data.job_id,
            sample_count_label = pluralize(len(run_data.samples), "sample"), generated_at = run_data.generated_at,
        ),
        *render_run_overview_section(run_data),
        *render_job_metadata_section(run_data),
        *render_sample_results_table_section(run_data),
        *render_sample_status_section(run_data, report_path),
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
