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
import pymupdf
import pytest
from mitopipeline.utils.reporting.markdown_report import render_sample_report
from mitopipeline.utils.reporting.models import SampleReportData
from mitopipeline.utils.reporting.pdf_export import export_markdown_to_pdf
from mitopipeline.utils.statistics.models import AssemblyStatistics

def _extract_text(pdf_path: Path) -> str:
    with pymupdf.open(pdf_path) as doc:
        return "\n".join(page.get_text() for page in doc)

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

def test_special_characters_in_a_data_value_survive_unmangled_in_the_real_pdf(tmp_path: Path) -> None:
    """Regression (Emil's explicit audit-escaping request): a data
    value carrying every Markdown-special character Typst could
    misinterpret (_, *, #, $, @) must appear in the final compiled
    PDF's own text layer exactly as supplied -- not doubled, not
    stripped, not reinterpreted as emphasis/heading/math/raw-typst
    markup. Reads the real PDF back with pymupdf (no poppler-utils on
    this host) rather than asserting against the Markdown source,
    since the whole point is verifying what Typst actually rendered."""
    raw_value = "sample_ID *bold* #tag $5 @mention"
    sample = SampleReportData(
        sample_id = "SAMPLE_02", metadata = {"notes": raw_value}, fastp_stats = None,
        assembly = None, blast = None, annotation = None, overall_verdict = None,
    )
    markdown_path = tmp_path / "SAMPLE_02.report.md"
    markdown_path.write_text(render_sample_report(sample, markdown_path))

    result = export_markdown_to_pdf(markdown_path)

    assert result.success is True
    text = _extract_text(result.pdf_path)
    assert raw_value in text
    assert "sample__ID" not in text
    assert "**bold**" not in text

def test_sample_with_a_visualization_image_still_exports_a_pdf(tmp_path: Path) -> None:
    """Regression: pandoc's typst PDF engine compiles through an
    intermediate file in its own temp directory (pandoc's
    "toPdfViaTempFile"), so the map figure's path -- previously
    relative to the report's own directory -- resolved against the
    wrong directory there and Typst rejected it ("would escape the
    project root"), confirmed via a real compile against the
    gargammel_circular/gargammel_fragmented fixture outputs: both
    samples with a real map image silently lost their PDF (markdown
    still rendered, pdf_path was None) while the map-less incomplete
    sample exported fine. markdown_report.py now passes the image's
    real absolute path and pdf_export.py widens Typst's sandbox root
    to "/" to permit it."""
    png_dir = tmp_path / "visualization"
    png_dir.mkdir()
    # A 1x1 PNG is enough -- this test is about path resolution, not
    # the image's own content, which the map-generating tool owns.
    visualization_png = png_dir / "SAMPLE_03.circular_map.png"
    visualization_png.write_bytes(bytes.fromhex(
        "89504e470d0a1a0a0000000d4948445200000001000000010802000000"
        "907753de0000000c49444154789c63f8cfc0000003010100c9fe92ef00"
        "00000049454e44ae426082"
    ))
    sample = SampleReportData(
        sample_id = "SAMPLE_03", metadata = {}, fastp_stats = None,
        assembly = AssemblyStatistics(sample_id = "SAMPLE_03", status = "circular", source_attempt = 0,
                                      length_bp = 16596, contig_count = 1, gc_content_percent = 45.0),
        blast = None, annotation = None, overall_verdict = "PASS",
        visualization_png = visualization_png,
    )
    report_dir = tmp_path / "13-report_pdf" / "SAMPLE_03"
    report_dir.mkdir(parents = True)
    markdown_path = report_dir / "SAMPLE_03.report.md"
    markdown_path.write_text(render_sample_report(sample, markdown_path))

    result = export_markdown_to_pdf(markdown_path)

    assert result.success is True, result.stderr
    assert result.pdf_path.read_bytes()[:4] == b"%PDF"

def test_missing_markdown_file_raises(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        export_markdown_to_pdf(tmp_path / "does_not_exist.md")

def test_missing_executable_raises(tmp_path: Path) -> None:
    markdown_path = tmp_path / "report.md"
    markdown_path.write_text("# Hello\n")
    with pytest.raises(FileNotFoundError):
        export_markdown_to_pdf(markdown_path, pandoc_bin = "definitely-not-a-real-binary")
