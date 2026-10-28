"""
test_methods.py

Unit tests for the live methods-section generator (Emil's explicit
request): real tool names/settings from resolved config JSON, and
rescue history reflected only when it actually happened.
"""

# Imports
from mitopipeline.utils.reporting.methods import render_methods_text
from mitopipeline.utils.reporting.compute_metrics import MachineInfo, StageComputeSummary
from mitopipeline.utils.reporting.models import (
    AnnotationConfigSummary,
    AssemblyConfigSummary,
    BlastConfigSummary,
    PhylogenySummary,
    RescueAttemptSummary,
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

def _machine_info(**overrides) -> MachineInfo:
    defaults = dict(
        os_name = "Linux", os_release = "6.18.54.1", architecture = "x86_64",
        cpu_model = "Intel(R) Xeon(R) Platinum 8259CL CPU @ 2.50GHz",
        cpu_physical_cores = 4, cpu_logical_threads = 8, cpu_clock_mhz = 2500.0,
        total_memory_gb = 16.0, memory_speed_mhz = 3200, memory_type = "DDR4",
        gpu_model = None, python_version = "3.12.3",
    )
    defaults.update(overrides)
    return MachineInfo(**defaults)

def test_methods_text_describes_the_real_compute_environment() -> None:
    sample = _sample(machine_info = _machine_info())
    text = render_methods_text(sample)
    assert "Linux x86_64" in text
    assert "Intel(R) Xeon(R) Platinum 8259CL CPU @ 2.50GHz" in text
    assert "4 physical core(s), 8 logical thread(s)" in text
    assert "2500 MHz" in text
    assert "16.0 GiB RAM" in text
    assert "DDR4" in text
    assert "3200 MHz" in text
    assert "Python 3.12.3" in text

def test_methods_text_mentions_gpu_with_not_utilized_disclaimer_when_present() -> None:
    sample = _sample(machine_info = _machine_info(gpu_model = "NVIDIA GeForce RTX 3080"))
    text = render_methods_text(sample)
    assert "NVIDIA GeForce RTX 3080" in text
    assert "not utilized by any part of this pipeline" in text

def test_methods_text_omits_gpu_clause_when_no_gpu_detected() -> None:
    sample = _sample(machine_info = _machine_info(gpu_model = None))
    text = render_methods_text(sample)
    assert "GPU" not in text

def test_methods_text_omits_environment_sentence_when_machine_info_unavailable() -> None:
    text = render_methods_text(_sample())
    assert "Computation was performed on" not in text

def test_methods_text_cites_real_fastp_and_fastqc_versions() -> None:
    sample = _sample(
        fastp_stats = {"reads_in": 1000, "reads_out": 900, "read_retention_rate": 0.9},
        compute_summary = (
            StageComputeSummary(stage = "trimming", tools = "fastp", invocation_count = 1,
                                total_runtime_seconds = 1.0, total_cpu_seconds = 1.0, peak_memory_mb = 10.0,
                                tool_versions = {"fastp": "fastp 1.4.0"}),
        ),
    )
    text = render_methods_text(sample)
    assert "fastp 1.4.0" in text

def test_methods_text_cites_real_getorganelle_version_for_the_primary_attempt() -> None:
    sample = _sample(
        assembly = AssemblyStatistics(sample_id = "s", status = "circular", source_attempt = 0,
                                      length_bp = 16596, contig_count = 1, gc_content_percent = 45.0),
        assembly_config = AssemblyConfigSummary(organelle_type = "animal_mt", seeded = False, seed_file = None),
        compute_summary = (
            StageComputeSummary(stage = "04d-assembly", tools = "getorganelle", invocation_count = 1,
                                total_runtime_seconds = 1.0, total_cpu_seconds = 1.0, peak_memory_mb = 10.0,
                                tool_versions = {"getorganelle": "GetOrganelle v1.7.7.1"}),
        ),
    )
    text = render_methods_text(sample)
    assert "GetOrganelle 1.7.7.1" in text

def test_methods_text_cites_distinct_real_versions_for_rescue_tools() -> None:
    """The nuanced case: rescue invokes minimap2, samtools, and a
    second GetOrganelle re-invocation -- each cited with its own real
    version from its own real 07-assembly_rescue metrics, separate
    from the primary attempt's own 04d-assembly GetOrganelle version."""
    sample = _sample(
        assembly = AssemblyStatistics(sample_id = "s", status = "fragmented", source_attempt = 1,
                                      length_bp = 16000, contig_count = 1, gc_content_percent = 44.0),
        assembly_config = AssemblyConfigSummary(organelle_type = "animal_mt", seeded = False, seed_file = None),
        rescue_attempts = (RescueAttemptSummary(attempt = 1, attempted = True, reason = "triggered",
                                                profile_name = "reduced_coverage", status = "fragmented"),),
        compute_summary = (
            StageComputeSummary(stage = "04d-assembly", tools = "getorganelle", invocation_count = 1,
                                total_runtime_seconds = 1.0, total_cpu_seconds = 1.0, peak_memory_mb = 10.0,
                                tool_versions = {"getorganelle": "GetOrganelle v1.7.7.0"}),
            StageComputeSummary(stage = "07-assembly_rescue", tools = "minimap2, samtools, getorganelle",
                                invocation_count = 3, total_runtime_seconds = 3.0, total_cpu_seconds = 3.0,
                                peak_memory_mb = 30.0, tool_versions = {
                                    "minimap2": "2.31-r1302", "samtools": "samtools 1.24",
                                    "getorganelle": "GetOrganelle v1.7.7.1",
                                }),
        ),
    )
    text = render_methods_text(sample)
    assert "GetOrganelle 1.7.7.0 (" in text  # the primary attempt's own version
    assert "minimap2 2.31-r1302" in text
    assert "samtools 1.24" in text
    assert "GetOrganelle 1.7.7.1):" in text  # the rescue re-invocation's own (different) version

def test_methods_text_cites_distinct_rescue_tool_versions_in_the_fallback_branch() -> None:
    """Same nuance, the other real code path: rescue_attempts is empty
    but source_attempt > 0 (the older/simpler fallback description)."""
    sample = _sample(
        assembly = AssemblyStatistics(sample_id = "s", status = "circular", source_attempt = 1,
                                      length_bp = 16596, contig_count = 1, gc_content_percent = 45.0),
        assembly_config = AssemblyConfigSummary(organelle_type = "animal_mt", seeded = False, seed_file = None),
        compute_summary = (
            StageComputeSummary(stage = "07-assembly_rescue", tools = "minimap2, samtools, getorganelle",
                                invocation_count = 3, total_runtime_seconds = 3.0, total_cpu_seconds = 3.0,
                                peak_memory_mb = 30.0, tool_versions = {
                                    "minimap2": "2.31-r1302", "samtools": "samtools 1.24",
                                    "getorganelle": "GetOrganelle v1.7.7.1",
                                }),
        ),
    )
    text = render_methods_text(sample)
    assert "minimap2 2.31-r1302" in text
    assert "samtools 1.24" in text
    assert "GetOrganelle 1.7.7.1)" in text

def test_methods_text_cites_real_blast_version() -> None:
    """Regression: BlastRunner's own real tool_name is "blast" (its
    generic BaseTool log-context label), not "blastn" (only the real
    build_command()'s own argv[0], the actual binary) -- the lookup
    key must match the real recorded tool_name, not the binary name."""
    sample = _sample(
        blast_config = BlastConfigSummary(database = "actinopterygii_mt_refs", mode = "blastn"),
        compute_summary = (
            StageComputeSummary(stage = "05b-blast", tools = "blast", invocation_count = 1,
                                total_runtime_seconds = 1.0, total_cpu_seconds = 1.0, peak_memory_mb = 10.0,
                                tool_versions = {"blast": "blastn: 2.16.0+"}),
        ),
    )
    text = render_methods_text(sample)
    assert "BLAST 2.16.0+ (" in text

def test_methods_text_cites_real_mitos2_version() -> None:
    sample = _sample(
        annotation_config = AnnotationConfigSummary(genetic_code = 2, refseqver = "refseq89m"),
        compute_summary = (
            StageComputeSummary(stage = "08b-annotation", tools = "mitos2", invocation_count = 1,
                                total_runtime_seconds = 1.0, total_cpu_seconds = 1.0, peak_memory_mb = 10.0,
                                tool_versions = {"mitos2": "runmitos.py 2.1.9"}),
        ),
    )
    text = render_methods_text(sample)
    assert "MITOS2 2.1.9 (" in text

def test_methods_text_mentions_mafft_and_cites_real_mafft_and_iqtree_versions() -> None:
    sample = _sample(
        phylogeny = PhylogenySummary(success = True, tip_count = 7, newick_path = None, reason = "ok", model = "MFP"),
        compute_summary = (
            StageComputeSummary(stage = "11-phylogeny", tools = "mafft, iqtree", invocation_count = 2,
                                total_runtime_seconds = 2.0, total_cpu_seconds = 2.0, peak_memory_mb = 20.0,
                                tool_versions = {"mafft": "v7.526 (2024/Apr/26)", "iqtree": "IQ-TREE version 3.1.2 for Linux"}),
        ),
    )
    text = render_methods_text(sample)
    assert "Sequences were aligned with MAFFT 7.526" in text
    assert "IQ-TREE 3.1.2" in text

def test_methods_text_omits_version_clauses_when_no_compute_summary_available() -> None:
    """A sample with no real compute_summary (e.g. metrics weren't
    passed through) must still render clean, valid prose -- no bare
    double-space or empty-parens artifact from a missing version."""
    sample = _sample(
        fastp_stats = {
            "reads_in": 1000, "reads_out": 900, "read_retention_rate": 0.9,
            "qualified_quality_phred": 20,
        },
        blast_config = BlastConfigSummary(database = "actinopterygii_mt_refs", mode = "blastn"),
    )
    text = render_methods_text(sample)
    assert "fastp (minimum base quality" in text
    assert "fastp  (" not in text
    assert "BLAST (`blastn`" in text
    assert "BLAST  (" not in text
