"""
test_markdown_report.py

Unit tests for per-sample Markdown report rendering (#358/#359):
every real section renders from typed statistics, and missing data
in any section renders a clear "not available"/"not run" note rather
than crashing or fabricating a value.
"""

# Imports
from pathlib import Path
import pytest
from mitopipeline.utils.reporting.markdown_report import (
    escape_markdown_cell,
    escape_markdown_text,
    file_link,
    format_number_in_text,
    format_scientific_name,
    humanize_enum,
    humanize_key,
    humanize_list,
    markdown_table,
    pluralize,
    render_aggregated_data_section,
    render_phylogeny_section,
    render_read_qc_section,
    render_sample_report,
    status_badge,
    wrap_cell_text,
)
from mitopipeline.utils.reporting.models import PhylogenySummary, ReadQCSummary, SampleReportData
from mitopipeline.utils.statistics.models import AnnotationStatistics, AssemblyStatistics, BlastStatistics

def _empty_sample(sample_id: str = "sample1") -> SampleReportData:
    return SampleReportData(
        sample_id = sample_id, metadata = {"genus": "Danio", "species": "rerio"},
        fastp_stats = None, assembly = None, blast = None, annotation = None, overall_verdict = None,
    )

def test_markdown_table_escapes_header_but_passes_body_through() -> None:
    """markdown_table() escapes headers itself, but body cells are the
    caller's responsibility to pre-format/pre-escape (so raw Typst
    badge/linebreak markup can pass through unescaped)."""
    table = markdown_table(["A|B"], [[escape_markdown_cell("x|y")]])
    assert "A\\|B" in table
    assert "x\\|y" in table

def test_escape_markdown_text_escapes_every_special_character() -> None:
    raw = "a_b*c#d$e@f"
    escaped = escape_markdown_text(raw)
    assert escaped == "a\\_b\\*c\\#d\\$e\\@f"

def test_escape_markdown_cell_still_escapes_pipes() -> None:
    assert escape_markdown_cell("a|b") == "a\\|b"

def test_format_number_in_text_adds_thousands_separator_before_bp() -> None:
    assert format_number_in_text("16596bp is within range") == "16,596 bp is within range"

def test_format_number_in_text_formats_interval_bounds() -> None:
    """Regression: the length-bounds validation gate's own real reason
    string puts "bp" after the closing bracket, not immediately after
    either bound (e.g. "16596bp is within the configured interval
    [14000, 22000]bp") -- confirmed via a real compiled report that
    only the first number was getting thousands separators."""
    text = format_number_in_text("16596bp is within the configured interval [14000, 22000]bp")
    assert text == "16,596 bp is within the configured interval [14,000, 22,000] bp"

def test_format_number_in_text_formats_interval_bounds_after_escaping() -> None:
    """Every real caller runs escape_markdown_cell() before this
    function, which backslash-escapes `[`/`]` -- confirmed via a real
    compiled report that the interval was still unformatted until the
    regex above also tolerated the escaped brackets."""
    text = format_number_in_text(escape_markdown_cell("16596bp is within the configured interval [14000, 22000]bp"))
    assert text == "16,596 bp is within the configured interval \\[14,000, 22,000\\] bp"

def test_humanize_key_converts_snake_case() -> None:
    assert humanize_key("expected_assembly_status") == "Expected assembly status"
    assert humanize_key("sample_id") == "Sample id"

def test_humanize_enum_converts_screaming_snake_case() -> None:
    assert humanize_enum("MANUAL_REVIEW_REQUIRED") == "Manual review required"
    assert humanize_enum("READY") == "Ready"
    assert humanize_enum("PASS") == "Pass"

def test_humanize_list_deduplicates_and_humanizes() -> None:
    result = humanize_list("assembly_fasta;annotation_gff;assembly_fasta")
    assert result == "Assembly fasta, Annotation gff"

def test_humanize_list_empty_string() -> None:
    assert humanize_list("") == ""

def test_wrap_cell_text_leaves_short_values_alone() -> None:
    assert wrap_cell_text("short") == "short"

def test_wrap_cell_text_breaks_long_values_with_a_real_linebreak() -> None:
    wrapped = wrap_cell_text("Manual review required", width = 10)
    assert "`#linebreak()`{=typst}" in wrapped
    assert "Manual" in wrapped and "required" in wrapped

def test_status_badge_renders_a_colored_badge_for_known_statuses() -> None:
    badge = status_badge("PASS")
    assert "#box(" in badge
    assert "Pass" in badge

def test_status_badge_falls_back_to_escaped_text_for_unknown_values() -> None:
    badge = status_badge("something_else")
    assert "#box(" not in badge
    assert "Something else" in badge

def test_pluralize_singular_and_plural() -> None:
    assert pluralize(1, "sample") == "1 sample"
    assert pluralize(3, "sample") == "3 samples"
    assert pluralize(0, "sample") == "0 samples"

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
    assert "No assembly: contig selection did not find a usable" in text
    # With no assembly, BLAST/annotation/validation are all "not run",
    # not their own individual "not available" messages.
    assert "_Not run: no assembly._" in text
    assert "_Circular genome visualization was not available._" in text
    assert "_Read trimming statistics were not available._" in text

def test_render_sample_report_includes_real_typst_preamble(tmp_path: Path) -> None:
    text = render_sample_report(_empty_sample(), tmp_path / "sample1.report.md")
    assert '#set page(' in text
    assert 'margin: (x: 0.5in, y: 0.5in)' in text
    assert 'New Computer Modern' in text
    assert '#set par(justify: true' in text
    assert '#set table(' in text

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
    assert "Circular" in text

def test_render_sample_report_never_uses_the_word_candidate(tmp_path: Path) -> None:
    """Regression (Emil's explicit request, both for report text and
    report-rendering code): "candidate" must never appear in
    report-generated prose. Caught for real in the no-assembly
    callout ("...did not find a usable circular or fragmented
    candidate..."), which is report-generated text, not pass-through
    manifest data -- unlike a sample's own raw Notes field, which
    this check does not touch."""
    text = render_sample_report(_empty_sample(), tmp_path / "sample1.report.md")
    assert "candidate" not in text.lower()

def test_render_phylogeny_section_shows_not_available_when_missing(tmp_path: Path) -> None:
    lines = render_phylogeny_section(_empty_sample(), tmp_path / "sample1.report.md")
    assert "_Phylogenetic placement was not available._" in lines

def test_render_phylogeny_section_shows_a_readable_reason_not_a_raw_exit_code(tmp_path: Path) -> None:
    """Regression: a raw "iqtree exited 2" must never reach the
    report -- run_iqtree.py itself writes a readable reason before
    ever invoking the tool when there are too few sequences (now
    rare per-sample, since a sample's own dataset is its genome plus
    up to 6 live-fetched references, but still a real, possible
    outcome when a sample has very few BLAST hits)."""
    sample = SampleReportData(
        sample_id = "sample1", metadata = {}, fastp_stats = None,
        assembly = None, blast = None, annotation = None, overall_verdict = None,
        phylogeny = PhylogenySummary(
            success = False, tip_count = 2, newick_path = None,
            reason = "Skipped: at least 3 assemblies are required to build a tree, but only 2 passed validation.",
        ),
    )
    text = "\n".join(render_phylogeny_section(sample, tmp_path / "sample1.report.md"))
    assert "at least 3 assemblies are required" in text
    assert "iqtree exited" not in text

def test_render_phylogeny_section_shows_the_real_tip_count_and_model(tmp_path: Path) -> None:
    sample = SampleReportData(
        sample_id = "sample1", metadata = {}, fastp_stats = None,
        assembly = None, blast = None, annotation = None, overall_verdict = None,
        phylogeny = PhylogenySummary(success = True, tip_count = 7, newick_path = None, reason = "completed", model = "GTR+G"),
    )
    text = "\n".join(render_phylogeny_section(sample, tmp_path / "sample1.report.md"))
    assert "GTR+G" in text
    assert "| Tree tips | 7 |" in text

def test_render_sample_report_includes_a_phylogeny_section(tmp_path: Path) -> None:
    text = render_sample_report(_empty_sample(), tmp_path / "sample1.report.md")
    assert "## Phylogenetic Context" in text

def test_render_read_qc_section_shows_not_available_when_missing() -> None:
    lines = render_read_qc_section(_empty_sample())
    assert "_FastQC quality control was not available._" in lines

def test_render_read_qc_section_shows_real_stats_per_readset(tmp_path: Path) -> None:
    sample = SampleReportData(
        sample_id = "sample1", metadata = {}, fastp_stats = None,
        assembly = None, blast = None, annotation = None, overall_verdict = None,
        raw_read_qc = ReadQCSummary(total_sequences = 25014, percent_gc = 40.0, sequence_length = "100", overall_verdict = "PASS"),
        trimmed_read_qc = ReadQCSummary(total_sequences = 25014, percent_gc = 40.0, sequence_length = "100", overall_verdict = "WARN"),
        recruited_read_qc = None,
    )
    text = "\n".join(render_read_qc_section(sample))
    assert "### Raw Reads" in text
    assert "### Trimmed Reads" in text
    assert "### Recruited Reads (used for assembly)" in text
    assert "25,014" in text
    assert "100 bp" in text
    assert "40.00%" in text
    # One passes, one warns -- both badges render, never collapsed to a single value.
    assert "Not available" in text  # the missing recruited-reads row

def test_render_read_qc_section_embeds_the_representative_figure(tmp_path: Path) -> None:
    image_path = tmp_path / "raw_per_base_quality.png"
    image_path.write_bytes(b"fake png bytes")
    sample = SampleReportData(
        sample_id = "sample1", metadata = {}, fastp_stats = None,
        assembly = None, blast = None, annotation = None, overall_verdict = None,
        raw_read_qc = ReadQCSummary(
            total_sequences = 100, percent_gc = 40.0, sequence_length = "100",
            overall_verdict = "PASS", representative_image_png = image_path,
        ),
    )
    text = "\n".join(render_read_qc_section(sample))
    assert "#figure(" in text
    assert str(image_path) in text

def test_render_sample_report_includes_the_read_qc_section(tmp_path: Path) -> None:
    text = render_sample_report(_empty_sample(), tmp_path / "sample1.report.md")
    assert "## Read Quality Control" in text

def test_render_aggregated_data_section_shows_not_available_when_missing(tmp_path: Path) -> None:
    lines = render_aggregated_data_section(_empty_sample(), tmp_path / "sample1.report.md")
    assert "_Aggregated data was not available for this run._" in lines

def test_render_aggregated_data_section_links_real_existing_files(tmp_path: Path) -> None:
    aggregated_dir = tmp_path / "aggregated_data"
    assembly = aggregated_dir / "assemblies" / "sample1.assembly.fasta"
    assembly.parent.mkdir(parents = True)
    assembly.write_text(">sample1\nACGT\n")
    gff = aggregated_dir / "annotations" / "gff" / "sample1.annotation.gff"
    gff.parent.mkdir(parents = True)
    gff.write_text("##gff-version 3\n")

    sample = SampleReportData(
        sample_id = "sample1", metadata = {}, fastp_stats = None,
        assembly = None, blast = None, annotation = None, overall_verdict = None,
        aggregated_data_dir = aggregated_dir,
    )
    text = "\n".join(render_aggregated_data_section(sample, tmp_path / "sample1.report.md"))
    assert "Assembly FASTA" in text
    assert "Annotation GFF" in text
    # Files that don't exist (no annotation map/phylogeny for this sample) are simply omitted, not broken links.
    assert "Annotated genome map" not in text

def test_render_aggregated_data_section_links_real_archival_materials(tmp_path: Path) -> None:
    aggregated_dir = tmp_path / "aggregated_data"
    archival = aggregated_dir / "archival"
    (archival).mkdir(parents = True)
    (archival / "sra_metadata.tsv").write_text("sample_name\tstatus\nsample1\tREADY\n")
    (archival / "biosample_metadata.tsv").write_text("sample_name\tstatus\nsample1\tREADY\n")
    (archival / "submission_checklist.md").write_text("# checklist\n")
    genbank_dir = archival / "genbank" / "sample1"
    genbank_dir.mkdir(parents = True)
    (genbank_dir / "sample1.fasta").write_text(">sample1\nACGT\n")
    (genbank_dir / "sample1.gff").write_text("##gff-version 3\n")

    sample = SampleReportData(
        sample_id = "sample1", metadata = {}, fastp_stats = None,
        assembly = None, blast = None, annotation = None, overall_verdict = None,
        aggregated_data_dir = aggregated_dir,
    )
    text = "\n".join(render_aggregated_data_section(sample, tmp_path / "sample1.report.md"))
    assert "SRA metadata table" in text
    assert "BioSample metadata table" in text
    assert "GenBank assembly FASTA" in text
    assert "GenBank annotation GFF" in text
    assert "Submission checklist" in text
    # Files that don't exist (no readiness report/checksums produced in this fixture) are simply omitted.
    assert "Submission readiness report" not in text
    assert "MD5 checksums" not in text

def test_render_aggregated_data_section_handles_a_sample_with_nothing_aggregated(tmp_path: Path) -> None:
    sample = SampleReportData(
        sample_id = "sample_bare", metadata = {}, fastp_stats = None,
        assembly = None, blast = None, annotation = None, overall_verdict = None,
        aggregated_data_dir = tmp_path / "aggregated_data",
    )
    lines = render_aggregated_data_section(sample, tmp_path / "sample_bare.report.md")
    assert "_No aggregated files were available for this sample._" in lines

def test_render_sample_report_includes_methods_section_exactly_once(tmp_path: Path) -> None:
    text = render_sample_report(_empty_sample(), tmp_path / "sample1.report.md")
    assert text.count("## Methods") == 1
    assert "GetOrganelle" in text or "was not configured" in text

def test_render_sample_report_shows_title_block_fields_exactly_once(tmp_path: Path) -> None:
    text = render_sample_report(_empty_sample("sample1"), tmp_path / "sample1.report.md")
    assert text.count("Sample ID: sample1") == 1
    assert text.count("Generated ") == 1

def test_render_sample_report_annotation_breaks_down_gene_types(tmp_path: Path) -> None:
    sample = SampleReportData(
        sample_id = "sample1", metadata = {}, fastp_stats = None,
        assembly = AssemblyStatistics(sample_id = "sample1", status = "circular", source_attempt = 0,
                                      length_bp = 16596, contig_count = 1, gc_content_percent = 45.5),
        blast = None,
        annotation = AnnotationStatistics(sample_id = "sample1", gene_count = 13, core_gene_count = 13,
                                          trna_count = 22, rrna_count = 2),
        overall_verdict = None,
    )
    text = render_sample_report(sample, tmp_path / "sample1.report.md")
    assert "13 / 13" in text
    assert "22 / 22" in text
    assert "2 / 2" in text
    assert "Total gene features" not in text

def test_render_sample_report_blast_resolves_na_scientific_name(tmp_path: Path) -> None:
    sample = SampleReportData(
        sample_id = "sample1", metadata = {}, fastp_stats = None,
        assembly = AssemblyStatistics(sample_id = "sample1", status = "circular", source_attempt = 0,
                                      length_bp = 16596, contig_count = 1, gc_content_percent = 45.5),
        blast = BlastStatistics(sample_id = "sample1", top_accession = "ref|NC_002333.2|", top_scientific_name = "N/A",
                                top_identity_percent = 100.0, top_query_coverage_percent = 100.0, top_evalue = 0.0),
        annotation = None, overall_verdict = None,
    )
    text = render_sample_report(sample, tmp_path / "sample1.report.md")
    assert "Taxon unresolved" in text
    assert "| Scientific name | N/A |" not in text
