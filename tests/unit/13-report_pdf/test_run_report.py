"""
test_run_report.py

Unit tests for the run-level Markdown report (#360/#441): PASS/
REVIEW/FAIL totals and per-sample status, including mixed outcomes.
Each sample's own real phylogenetic tree (restoring v0.1.0's
per-sample design) is tested in test_markdown_report.py instead --
there is no batch-level phylogeny artifact anymore.
"""

# Imports
from pathlib import Path
from mitopipeline.utils.reporting.compute_metrics import MachineInfo, StageComputeSummary
from mitopipeline.utils.reporting.models import JobMetadataSummary, RunReportData, SampleReportData
from mitopipeline.utils.reporting.run_report import render_run_report

def _sample(sample_id: str, verdict: str | None) -> SampleReportData:
    return SampleReportData(
        sample_id = sample_id, metadata = {}, fastp_stats = None,
        assembly = None, blast = None, annotation = None, overall_verdict = verdict,
    )

def test_run_overview_counts_mixed_verdicts(tmp_path: Path) -> None:
    run_data = RunReportData(
        job_id = "job1", generated_at = "27 Oct 2026, 00:00 UTC",
        samples = (_sample("s1", "PASS"), _sample("s2", "REVIEW"), _sample("s3", "FAIL"), _sample("s4", None)),
    )
    text = render_run_report(run_data, tmp_path / "run.report.md")
    assert "job1" in text
    assert "4 samples" in text
    assert " 1" in text  # each tier's own count appears somewhere in the overview table

def test_sample_status_section_lists_every_sample(tmp_path: Path) -> None:
    run_data = RunReportData(job_id = "job1", generated_at = "now", samples = (_sample("s1", "PASS"),))
    text = render_run_report(run_data, tmp_path / "run.report.md")
    assert "s1" in text
    assert "Pass" in text

def test_sample_status_section_links_the_real_nested_pdf_path(tmp_path: Path) -> None:
    """Regression: the per-sample PDF lives at report_path.parent/{sample}/{sample}.report.pdf,
    not as a sibling of the run report -- a prior version of this
    link always showed "Not available" even when the PDF existed."""
    report_path = tmp_path / "run.report.md"
    pdf_path = tmp_path / "s1" / "s1.report.pdf"
    pdf_path.parent.mkdir(parents = True)
    pdf_path.write_bytes(b"%PDF-1.7")
    run_data = RunReportData(job_id = "job1", generated_at = "now", samples = (_sample("s1", "PASS"),))
    text = render_run_report(run_data, report_path)
    assert "Not available" not in text.split("## Sample Status")[1].split("## NCBI Archival")[0]
    assert "s1.report.pdf" in text

def test_run_report_includes_real_typst_preamble(tmp_path: Path) -> None:
    run_data = RunReportData(job_id = "job1", generated_at = "now", samples = ())
    text = render_run_report(run_data, tmp_path / "run.report.md")
    assert 'margin: (x: 0.5in, y: 0.5in)' in text

def test_run_report_includes_the_per_sample_results_table(tmp_path: Path) -> None:
    run_data = RunReportData(job_id = "job1", generated_at = "now", samples = (_sample("s1", "PASS"),))
    text = render_run_report(run_data, tmp_path / "run.report.md")
    assert "## Per-Sample Results" in text

def test_per_sample_results_table_wraps_a_long_sample_id(tmp_path: Path) -> None:
    """Regression: this table packs 8 columns into the page width, and
    Typst's own pandoc-generated auto columns don't wrap long,
    space-free content -- confirmed via a real compiled report where
    a 21-character sample ID ("gargammel_fragmented") overflowed
    directly into the adjacent Status column with no gap at all
    ("gargammel_fragmentedFragmented")."""
    run_data = RunReportData(job_id = "job1", generated_at = "now", samples = (_sample("gargammel_fragmented", "REVIEW"),))
    text = render_run_report(run_data, tmp_path / "run.report.md")
    assert "gargammel_fragmentedFragmented" not in text
    assert "`#linebreak()`{=typst}" in text

def test_job_metadata_section_shows_not_available_when_missing(tmp_path: Path) -> None:
    run_data = RunReportData(job_id = "job1", generated_at = "now", samples = ())
    text = render_run_report(run_data, tmp_path / "run.report.md")
    assert "## Job Execution Metadata" in text
    assert "_Job lifecycle metadata was not available for this run._" in text

def test_job_metadata_section_shows_real_status_and_runtime(tmp_path: Path) -> None:
    job_metadata = JobMetadataSummary(
        status = "completed", start_time = "2026-10-28T09:28:39.920326", runtime_seconds = 217.66,
        sample_count = 3, attempt_count = 2, cumulative_runtime_seconds = 435.3, incomplete_attempts = 1,
    )
    run_data = RunReportData(job_id = "job1", generated_at = "now", samples = (), job_metadata = job_metadata)
    text = render_run_report(run_data, tmp_path / "run.report.md")
    assert "3m 37.7s" in text
    assert "Completed" in text
    assert "Incomplete attempts" in text

def test_job_metadata_section_shows_running_not_created_for_a_job_still_in_progress(tmp_path: Path) -> None:
    """Regression: this report is itself one of the last rules in the
    same Snakemake run it describes, so current_metadata.json's own
    terminal status is structurally never visible yet -- it always
    reads "created" with the original runtime_seconds=0.0 placeholder
    at report-generation time. Must read as "Running" with a real
    elapsed duration, never a bare, confusing "Created" / "0.0s"."""
    job_metadata = JobMetadataSummary(status = "created", start_time = "2026-10-28T09:28:39", runtime_seconds = 142.5)
    run_data = RunReportData(job_id = "job1", generated_at = "now", samples = (), job_metadata = job_metadata)
    text = render_run_report(run_data, tmp_path / "run.report.md")
    assert "Running" in text
    assert "Created" not in text
    assert "Elapsed so far" in text
    assert "2m 22.5s" in text

def test_job_metadata_section_shows_machine_info() -> None:
    from mitopipeline.utils.reporting.run_report import render_job_metadata_section
    info = MachineInfo(
        os_name = "Linux", os_release = "6.18.54.1", architecture = "x86_64",
        cpu_model = "Intel(R) Xeon(R) Platinum 8259CL CPU @ 2.50GHz",
        cpu_physical_cores = 4, cpu_logical_threads = 8, cpu_clock_mhz = 2500.0,
        total_memory_gb = 16.0, memory_speed_mhz = 3200, memory_type = "DDR4",
        gpu_model = None, python_version = "3.12.3",
    )
    run_data = RunReportData(job_id = "job1", generated_at = "now", samples = (), machine_info = info)
    text = "\n".join(render_job_metadata_section(run_data))
    assert "### Compute Environment" in text
    assert "Linux 6.18.54.1" in text
    assert "x86_64" in text
    assert "Intel(R) Xeon(R) Platinum 8259CL CPU @ 2.50GHz" in text
    assert "4 physical / 8 logical (threads)" in text
    assert "2500 MHz" in text
    assert "16.0 GiB" in text
    assert "DDR4 3200 MHz" in text
    assert "No GPU detected" in text
    assert "3.12.3" in text

def test_job_metadata_section_shows_gpu_with_not_utilized_disclaimer() -> None:
    from mitopipeline.utils.reporting.run_report import render_job_metadata_section
    info = MachineInfo(
        os_name = "Linux", os_release = "6.18.54.1", architecture = "x86_64",
        cpu_model = None, cpu_physical_cores = None, cpu_logical_threads = None, cpu_clock_mhz = None,
        total_memory_gb = None, memory_speed_mhz = None, memory_type = None,
        gpu_model = "NVIDIA GeForce RTX 3080", python_version = "3.12.3",
    )
    run_data = RunReportData(job_id = "job1", generated_at = "now", samples = (), machine_info = info)
    text = "\n".join(render_job_metadata_section(run_data))
    assert "NVIDIA GeForce RTX 3080" in text
    assert "not utilized by any part of this pipeline" in text

def test_job_metadata_section_shows_batch_wide_compute_summary_and_access_paths(tmp_path: Path) -> None:
    from mitopipeline.utils.reporting.run_report import render_job_metadata_section
    compute = (StageComputeSummary(stage = "04d-assembly", tools = "getorganelle", invocation_count = 3,
                                   total_runtime_seconds = 90.0, total_cpu_seconds = 80.0, peak_memory_mb = 2048.0),)
    run_data = RunReportData(job_id = "job1", generated_at = "now", samples = (), compute_summary = compute,
                             job_dir = Path("/outputs/job1"), global_log_file = Path("/outputs/job1/logs/00-job/job1.global.log"))
    text = "\n".join(render_job_metadata_section(run_data))
    assert "### Per-Stage Compute Summary" in text
    assert "getorganelle" in text
    assert "2.00 GB" in text
    assert "**Job output root:** `/outputs/job1`" in text
    assert "**Shared global log:**" in text
    assert "**Raw compute metrics:**" in text
    assert "**Job status history:**" in text
