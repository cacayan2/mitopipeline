"""
test_methods.py

Unit tests for the live methods-section generator (Emil's explicit
request): real tool names/settings from resolved config JSON, and
rescue history reflected only when it actually happened.
"""

# Imports
from mitopipeline.utils.reporting.methods import render_methods_text
from mitopipeline.utils.reporting.models import (
    AnnotationConfigSummary,
    AssemblyConfigSummary,
    BlastConfigSummary,
    SampleReportData,
)
from mitopipeline.utils.statistics.models import AssemblyStatistics, BlastStatistics

def _sample(**overrides) -> SampleReportData:
    defaults = dict(
        sample_id = "sample1", metadata = {}, fastp_stats = None,
        assembly = None, blast = None, annotation = None, overall_verdict = None,
    )
    defaults.update(overrides)
    return SampleReportData(**defaults)

def test_methods_text_describes_seeded_assembly_without_rescue() -> None:
    sample = _sample(
        assembly = AssemblyStatistics(sample_id = "s", status = "circular", source_attempt = 0,
                                      length_bp = 16596, contig_count = 1, gc_content_percent = 45.0),
        assembly_config = AssemblyConfigSummary(organelle_type = "animal_mt", seeded = True, seed_file = "12S_rRNA.fasta"),
    )
    text = render_methods_text(sample)
    assert "GetOrganelle" in text
    assert "seeded with a reference locus (`12S_rRNA.fasta`)" in text
    assert "without requiring rescue" in text
    assert "rescue reassembly" not in text

def test_methods_text_describes_unseeded_assembly() -> None:
    sample = _sample(assembly_config = AssemblyConfigSummary(organelle_type = "animal_mt", seeded = False, seed_file = None))
    text = render_methods_text(sample)
    assert "run unseeded" in text

def test_methods_text_mentions_rescue_when_it_happened() -> None:
    sample = _sample(
        assembly = AssemblyStatistics(sample_id = "s", status = "fragmented", source_attempt = 1,
                                      length_bp = 16000, contig_count = 1, gc_content_percent = 44.0),
        assembly_config = AssemblyConfigSummary(organelle_type = "animal_mt", seeded = False, seed_file = None),
    )
    text = render_methods_text(sample)
    assert "rescue reassembly (attempt 1" in text
    assert "minimap2/samtools" in text

def test_methods_text_describes_blast_identification() -> None:
    sample = _sample(
        blast_config = BlastConfigSummary(database = "actinopterygii_mt_refs", mode = "blastn"),
        blast = BlastStatistics(sample_id = "s", top_accession = "NC_1", top_scientific_name = "Danio rerio",
                                top_identity_percent = 99.0, top_query_coverage_percent = 98.0, top_evalue = 0.0),
    )
    text = render_methods_text(sample)
    assert "BLAST" in text
    assert "Danio rerio" in text
    assert "99.00% identity" in text

def test_methods_text_describes_annotation() -> None:
    sample = _sample(annotation_config = AnnotationConfigSummary(genetic_code = 2, refseqver = "refseq89m"))
    text = render_methods_text(sample)
    assert "MITOS2" in text
    assert "genetic code 2" in text
    assert "refseq89m" in text

def test_methods_text_includes_validation_verdict_when_available() -> None:
    sample = _sample(overall_verdict = "PASS")
    text = render_methods_text(sample)
    assert "overall verdict of PASS" in text

def test_methods_text_handles_fully_missing_configuration() -> None:
    text = render_methods_text(_sample())
    assert "GetOrganelle" in text
    assert "not configured" in text
