"""
test_pdf_export.py

Tests for Markdown-to-PDF export (#361/#442). Exercises the real
installed pandoc/typst binaries directly -- a mocked subprocess call
could never catch a real Typst syntax error in the embedded preamble/
title-block/callout-box markup, which is exactly the thing most
likely to break silently.
"""

# Imports
from pathlib import Path
import pytest
from mitopipeline.utils.reporting.markdown_report import render_sample_report
from mitopipeline.utils.reporting.models import SampleReportData
from mitopipeline.utils.reporting.pdf_export import export_markdown_to_pdf
from mitopipeline.utils.statistics.models import AssemblyStatistics

def test_real_report_renders_to_a_nonempty_pdf(tmp_path: Path) -> None:
    sample = SampleReportData(
        sample_id = "SAMPLE_01", metadata = {"genus": "Danio", "species": "rerio"}, fastp_stats = None,
        assembly = AssemblyStatistics(sample_id = "SAMPLE_01", status = "circular", source_attempt = 0,
                                      length_bp = 16596, contig_count = 1, gc_content_percent = 45.0),
        blast = None, annotation = None, overall_verdict = "PASS",
        validation_checks = ({"name": "length_bounds", "verdict": "PASS", "reason": "within bounds"},),
    )
    markdown_path = tmp_path / "SAMPLE_01.report.md"
    markdown_path.write_text(render_sample_report(sample, markdown_path))

    result = export_markdown_to_pdf(markdown_path)

    assert result.success is True
    assert result.pdf_path.is_file()
    assert result.pdf_path.stat().st_size > 0
    assert result.pdf_path.read_bytes()[:4] == b"%PDF"

def test_missing_markdown_file_raises(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        export_markdown_to_pdf(tmp_path / "does_not_exist.md")

def test_missing_executable_raises(tmp_path: Path) -> None:
    markdown_path = tmp_path / "report.md"
    markdown_path.write_text("# Hello\n")
    with pytest.raises(FileNotFoundError):
        export_markdown_to_pdf(markdown_path, pandoc_bin = "definitely-not-a-real-binary")
