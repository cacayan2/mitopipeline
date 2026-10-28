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
from mitopipeline.utils.statistics.models import AnnotationStatistics, AssemblyStatistics, BlastStatistics

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
    sample = _sample(
        annotation_config = AnnotationConfigSummary(genetic_code = 2, refseqver = "refseq89m"),
        annotation = AnnotationStatistics(sample_id = "s", gene_count = 13, core_gene_count = 13,
                                          trna_count = 22, rrna_count = 2),
    )
    text = render_methods_text(sample)
    assert "MITOS2" in text
    assert "genetic code 2" in text
    assert "refseq89m" in text
    assert "22 tRNA(s)" in text
    assert "2 rRNA(s)" in text

def test_methods_text_annotation_includes_mitofish_biomarker_alongside_genetic_code() -> None:
    sample = _sample(
        assembly_config = AssemblyConfigSummary(organelle_type = "animal_mt", seeded = True, seed_file = None,
                                                seed_gene = "COXI", seed_accession = "AC024175", seed_source = "mitofish"),
        annotation_config = AnnotationConfigSummary(genetic_code = 2, refseqver = "refseq89m"),
        annotation = AnnotationStatistics(sample_id = "s", gene_count = 13, core_gene_count = 13,
                                          trna_count = 22, rrna_count = 2),
    )
    text = render_methods_text(sample)
    assert "genetic code 2, MitoFish biomarker COXI" in text

def test_methods_text_annotation_omits_biomarker_when_unseeded() -> None:
    sample = _sample(
        assembly_config = AssemblyConfigSummary(organelle_type = "animal_mt", seeded = False, seed_file = None),
        annotation_config = AnnotationConfigSummary(genetic_code = 2, refseqver = "refseq89m"),
    )
    text = render_methods_text(sample)
    assert "genetic code 2, reference set" in text
    assert "MitoFish biomarker" not in text

def test_methods_text_includes_validation_verdict_when_available() -> None:
    sample = _sample(overall_verdict = "PASS")
    text = render_methods_text(sample)
    assert "plausibility verdict of Pass" in text

def test_methods_text_handles_fully_missing_configuration() -> None:
    text = render_methods_text(_sample())
    assert "GetOrganelle" in text
    assert "not configured" in text

def test_methods_text_distinguishes_skipped_from_not_configured() -> None:
    """Regression: a sample whose assembly genuinely ran (config was
    resolved) but produced nothing selected must read "skipped", not
    the misleading "not configured" (which implies the stage itself
    was turned off). fastp_stats is supplied so the trimming sentence
    (a genuinely separate, assembly-independent stage) doesn't itself
    say "not configured" and pollute this specific assertion."""
    sample = _sample(
        assembly_config = AssemblyConfigSummary(organelle_type = "animal_mt", seeded = False, seed_file = None),
        fastp_stats = {
            "reads_in": 1000, "reads_out": 900, "read_retention_rate": 0.9,
            "q20_rate_out": 0.98, "q30_rate_out": 0.95, "gc_content_out": 0.45, "low_quality_reads": 50,
        },
    )
    text = render_methods_text(sample)
    assert "was skipped: this sample had no selected assembly" in text
    assert "was not configured for this run" not in text

def test_methods_text_blast_rephrases_unresolved_taxon_instead_of_na() -> None:
    sample = _sample(
        blast_config = BlastConfigSummary(database = "actinopterygii_mt_refs", mode = "blastn"),
        blast = BlastStatistics(sample_id = "s", top_accession = "NC_1", top_scientific_name = "N/A",
                                top_identity_percent = 100.0, top_query_coverage_percent = 100.0, top_evalue = 0.0),
    )
    text = render_methods_text(sample)
    assert "The top match was N/A" not in text
    assert "`NC_1` (taxon unresolved)" in text

def test_methods_text_renders_blast_database_as_basename_not_full_path() -> None:
    """Regression: a real run's resolved BLAST database path is an
    absolute path under the execution host's own home directory
    (e.g. /home/eacac/.../resources/blast/actinopterygii_mt_refs) --
    confirmed via a real compiled report leaking it into the Methods
    prose. Only the basename should ever appear."""
    sample = _sample(
        blast_config = BlastConfigSummary(
            database = "/home/eacac/mitogenome-batch-processor/v0.2.0/resources/blast/actinopterygii_mt_refs",
            mode = "blastn",
        ),
    )
    text = render_methods_text(sample)
    assert "`actinopterygii_mt_refs`" in text
    assert "/home/" not in text

def test_methods_text_renders_seed_file_as_basename_not_full_path() -> None:
    sample = _sample(
        assembly_config = AssemblyConfigSummary(
            organelle_type = "animal_mt", seeded = True,
            seed_file = "/home/eacac/mitogenome-batch-processor/v0.2.0/resources/seeds/12S_rRNA.fasta",
        ),
    )
    text = render_methods_text(sample)
    assert "`12S_rRNA.fasta`" in text
    assert "/home/" not in text
