"""
test_markdown_report.py

Unit tests for per-sample Markdown report rendering (#358/#359):
every real section renders from typed statistics, and missing data
in any section renders a clear "not available" note rather than
crashing or fabricating a value.
"""

# Imports
from pathlib import Path
import pytest
from mitopipeline.utils.reporting.markdown_report import (
    file_link,
    format_scientific_name,
    markdown_table,
    render_sample_report,
)
from mitopipeline.utils.reporting.models import SampleReportData
from mitopipeline.utils.statistics.models import AnnotationStatistics, AssemblyStatistics, BlastStatistics

def _empty_sample(sample_id: str = "sample1") -> SampleReportData:
    return SampleReportData(
        sample_id = sample_id, metadata = {"genus": "Danio", "species": "rerio"},
        fastp_stats = None, assembly = None, blast = None, annotation = None, overall_verdict = None,
    )

def test_markdown_table_escapes_pipes_and_newlines() -> None:
    table = markdown_table(["A"], [["x|y\nz"]])
    assert "x\\|y z" in table

def test_format_scientific_name_from_metadata() -> None:
    assert format_scientific_name({"genus": "Danio", "species": "rerio"}) == "Danio rerio"

def test_format_scientific_name_falls_back_when_missing() -> None:
    assert format_scientific_name({}) == "Not provided"

def test_file_link_none_when_target_missing(tmp_path: Path) -> None:
    assert file_link("label", tmp_path / "missing.png", tmp_path / "report.md") is None

def test_file_link_present_when_target_exists(tmp_path: Path) -> None:
    target = tmp_path / "map.png"
    target.write_bytes(b"x")
    link = file_link("map", target, tmp_path / "report.md")
    assert link == "[map](map.png)"

def test_render_sample_report_with_no_data_shows_not_available_everywhere(tmp_path: Path) -> None:
    text = render_sample_report(_empty_sample(), tmp_path / "sample1.report.md")
    assert "sample1" in text
    assert "_Assembly statistics were not available._" in text
    assert "_BLAST identification results were not available._" in text
    assert "_Annotation statistics were not available._" in text
    assert "_Biological validation was not available._" in text
    assert "_Circular genome visualization was not available._" in text
    assert "_Read trimming statistics were not available._" in text

def test_render_sample_report_includes_real_typst_preamble(tmp_path: Path) -> None:
    text = render_sample_report(_empty_sample(), tmp_path / "sample1.report.md")
    assert '#set page(' in text
    assert 'margin: (x: 0.5in, y: 0.5in)' in text
    assert 'New Computer Modern' in text
    assert '#set par(justify: true)' in text

def test_render_sample_report_renders_real_assembly_stats(tmp_path: Path) -> None:
    sample = SampleReportData(
        sample_id = "sample1", metadata = {}, fastp_stats = None,
        assembly = AssemblyStatistics(sample_id = "sample1", status = "circular", source_attempt = 0,
                                      length_bp = 16596, contig_count = 1, gc_content_percent = 45.5),
        blast = None, annotation = None, overall_verdict = None,
    )
    text = render_sample_report(sample, tmp_path / "sample1.report.md")
    assert "16,596 bp" in text
    assert "45.50%" in text
    assert "circular" in text

def test_render_sample_report_includes_methods_section(tmp_path: Path) -> None:
    text = render_sample_report(_empty_sample(), tmp_path / "sample1.report.md")
    assert "## Methods" in text
    assert "GetOrganelle" in text or "was not configured" in text
