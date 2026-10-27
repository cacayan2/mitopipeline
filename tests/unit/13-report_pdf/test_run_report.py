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
        job_id = "job1", generated_at = "2026-10-27T00:00:00Z",
        samples = (_sample("s1", "PASS"), _sample("s2", "REVIEW"), _sample("s3", "FAIL"), _sample("s4", None)),
    )
    text = render_run_report(run_data, tmp_path / "run.report.md")
    assert "job1" in text
    assert "| PASS | 1 |" in text
    assert "| REVIEW | 1 |" in text
    assert "| FAIL | 1 |" in text
    assert "| Not assessed | 1 |" in text

def test_sample_status_section_lists_every_sample(tmp_path: Path) -> None:
    run_data = RunReportData(job_id = "job1", generated_at = "now", samples = (_sample("s1", "PASS"),))
    text = render_run_report(run_data, tmp_path / "run.report.md")
    assert "s1" in text
    assert "PASS" in text

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

def test_phylogeny_section_handles_failed_tree(tmp_path: Path) -> None:
    phylogeny = PhylogenySummary(success = False, tip_count = None, newick_path = None, reason = "not enough taxa")
    run_data = RunReportData(job_id = "job1", generated_at = "now", samples = (), phylogeny = phylogeny)
    text = render_run_report(run_data, tmp_path / "run.report.md")
    assert "not enough taxa" in text

def test_phylogeny_section_embeds_the_tree_figure_when_present(tmp_path: Path) -> None:
    figure_png = tmp_path / "tree_figure.png"
    figure_png.write_bytes(b"fake png bytes")
    phylogeny = PhylogenySummary(success = True, tip_count = 5, newick_path = None, reason = "completed", tree_figure_png = figure_png)
    run_data = RunReportData(job_id = "job1", generated_at = "now", samples = (), phylogeny = phylogeny)
    text = render_run_report(run_data, tmp_path / "run.report.md")
    assert "![Maximum-likelihood phylogeny]" in text

def test_phylogeny_section_omits_the_image_embed_without_a_tree_figure(tmp_path: Path) -> None:
    phylogeny = PhylogenySummary(success = True, tip_count = 5, newick_path = None, reason = "completed")
    run_data = RunReportData(job_id = "job1", generated_at = "now", samples = (), phylogeny = phylogeny)
    text = render_run_report(run_data, tmp_path / "run.report.md")
    assert "![Maximum-likelihood phylogeny]" not in text

def test_run_report_includes_real_typst_preamble(tmp_path: Path) -> None:
    run_data = RunReportData(job_id = "job1", generated_at = "now", samples = ())
    text = render_run_report(run_data, tmp_path / "run.report.md")
    assert 'margin: (x: 0.5in, y: 0.5in)' in text
