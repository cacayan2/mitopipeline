"""
test_archival_in_report.py

Unit tests for NCBI archival metrics appearing in the report (the
explicit "add associated metrics to the final report" request):
collecting a sample's real SRA/BioSample/GenBank status, and
rendering it in both the per-sample and run-level reports.
"""

# Imports
from pathlib import Path
from mitopipeline.utils.reporting.collect import collect_archival_summary
from mitopipeline.utils.reporting.markdown_report import render_archival_section, render_sample_report
from mitopipeline.utils.reporting.models import ArchivalSummary, RunReportData, SampleReportData
from mitopipeline.utils.reporting.run_report import render_archival_summary_section, render_run_report

def _write(path: Path, text: str) -> Path:
    path.parent.mkdir(parents = True, exist_ok = True)
    path.write_text(text)
    return path

def _empty_sample(sample_id: str = "sample1", archival: ArchivalSummary | None = None) -> SampleReportData:
    return SampleReportData(
        sample_id = sample_id, metadata = {}, fastp_stats = None,
        assembly = None, blast = None, annotation = None, overall_verdict = None, archival = archival,
    )

def test_collect_archival_summary_from_real_tables(tmp_path: Path) -> None:
    sra_tsv = _write(tmp_path / "sra.tsv", "sample_name\treview_status\tmissing_fields\nsample1\tREADY\t\n")
    biosample_tsv = _write(tmp_path / "biosample.tsv", "sample_name\treview_status\tmissing_fields\nsample1\tMANUAL_REVIEW_REQUIRED\tcollection_date\n")
    genbank_json = tmp_path / "sample1.genbank_record.json"
    genbank_json.write_text('{"review_status": "READY", "missing_files": "", "missing_metadata": ""}')

    summary = collect_archival_summary("sample1", sra_tsv, biosample_tsv, genbank_json)

    assert summary.sra_review_status == "READY"
    assert summary.biosample_review_status == "MANUAL_REVIEW_REQUIRED"
    assert summary.biosample_missing_fields == "collection_date"
    assert summary.genbank_review_status == "READY"
    assert summary.ready_for_review is False

def test_collect_archival_summary_none_when_nothing_supplied() -> None:
    assert collect_archival_summary("sample1") is None

def test_archival_summary_ready_when_all_sources_ready() -> None:
    summary = ArchivalSummary(sra_review_status = "READY", biosample_review_status = "READY", genbank_review_status = "READY")
    assert summary.ready_for_review is True

def test_render_archival_section_shows_not_available_when_missing() -> None:
    lines = render_archival_section(_empty_sample())
    assert "_Archival preparation was not available._" in lines

def test_render_archival_section_shows_real_status() -> None:
    archival = ArchivalSummary(sra_review_status = "READY", biosample_review_status = "MANUAL_REVIEW_REQUIRED",
                               biosample_missing_fields = "collection_date", genbank_review_status = "READY")
    lines = render_archival_section(_empty_sample(archival = archival))
    text = "\n".join(lines)
    assert "READY" in text and "MANUAL_REVIEW_REQUIRED" in text
    assert "collection_date" in text
    assert "Ready for submission review: No" in text

def test_sample_report_includes_archival_section(tmp_path: Path) -> None:
    text = render_sample_report(_empty_sample(), tmp_path / "sample1.report.md")
    assert "## NCBI archival readiness" in text

def test_run_report_archival_summary_shows_real_counts(tmp_path: Path) -> None:
    ready = ArchivalSummary(sra_review_status = "READY", biosample_review_status = "READY", genbank_review_status = "READY")
    not_ready = ArchivalSummary(sra_review_status = "MANUAL_REVIEW_REQUIRED")
    run_data = RunReportData(
        job_id = "job1", generated_at = "now",
        samples = (_empty_sample("s1", ready), _empty_sample("s2", not_ready)),
    )
    assert run_data.archival_ready_count == 1
    text = render_run_report(run_data, tmp_path / "run.report.md")
    assert "1 / 2" in text

def test_run_report_archival_summary_not_available_when_no_samples_have_it(tmp_path: Path) -> None:
    run_data = RunReportData(job_id = "job1", generated_at = "now", samples = (_empty_sample(),))
    lines = render_archival_summary_section(run_data)
    assert "_Archival preparation was not available for this run._" in lines
