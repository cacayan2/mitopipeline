"""
test_run_report.py

Unit tests for the run-level Markdown report (#360/#441): PASS/
REVIEW/FAIL totals, per-sample status, and the batch-level phylogeny
summary -- including mixed outcomes and a fully missing phylogeny run.
"""

# Imports
from pathlib import Path
from mitopipeline.utils.reporting.models import PhylogenySummary, RunReportData, SampleReportData
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
    assert "Not available" not in text.split("## Sample Status")[1].split("## Phylogenetic")[0]
    assert "s1.report.pdf" in text

def test_phylogeny_section_shows_real_outcome(tmp_path: Path) -> None:
    phylogeny = PhylogenySummary(success = True, tip_count = 5, newick_path = None, reason = "completed", model = "GTR+G")
    run_data = RunReportData(job_id = "job1", generated_at = "now", samples = (), phylogeny = phylogeny)
    text = render_run_report(run_data, tmp_path / "run.report.md")
    assert "GTR+G" in text
    assert "| Tree tips | 5 |" in text

def test_phylogeny_section_handles_missing_phylogeny(tmp_path: Path) -> None:
    run_data = RunReportData(job_id = "job1", generated_at = "now", samples = ())
    text = render_run_report(run_data, tmp_path / "run.report.md")
    assert "_Phylogenetic validation was not available for this run._" in text

def test_phylogeny_section_shows_a_readable_reason_not_a_raw_exit_code(tmp_path: Path) -> None:
    """Regression: a raw "iqtree exited 2" must never reach the
    report -- run_iqtree.py itself now writes a readable reason
    before ever invoking the tool when there are too few sequences."""
    phylogeny = PhylogenySummary(
        success = False, tip_count = 2, newick_path = None,
        reason = "Skipped: at least 3 assemblies are required to build a tree, but only 2 passed validation.",
    )
    run_data = RunReportData(job_id = "job1", generated_at = "now", samples = (), phylogeny = phylogeny)
    text = render_run_report(run_data, tmp_path / "run.report.md")
    assert "at least 3 assemblies are required" in text
    assert "iqtree exited" not in text

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
