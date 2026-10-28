"""
models.py

Data required to render one sample's report and one run's report
(#358-361, #441-443) -- each field traces back to a real result JSON
an earlier stage already wrote, never re-derived or fabricated.
"""

# Imports
from dataclasses import dataclass, field
from pathlib import Path
from mitopipeline.utils.statistics.models import AnnotationStatistics, AssemblyStatistics, BlastStatistics

@dataclass(frozen = True, slots = True)
class AssemblyConfigSummary:
    """
    The real, deterministic GetOrganelle settings used for this
    sample's selected attempt -- read from #04a's own resolved
    config, the live source the methods section cites rather than
    hardcoding tool names/settings as boilerplate.

    Attributes:
        organelle_type (str): The GetOrganelle "-F" organelle type used.
        seeded (bool): Whether a seed reference was available for this run.
        seed_file (str | None): The seed file used, or None when unseeded.
        seed_gene (str | None): The mitochondrial locus the seed was selected for (e.g. "COXI"), from the real seed_selection.json, or None when unavailable.
        seed_accession (str | None): The real source accession for the selected seed, or the fallback marker when GetOrganelle's own default database was used instead.
        seed_source (str | None): Which reference source produced the seed ("mitofish" or "fallback_default_database"), or None when unavailable.
        seed_fallback_reason (str | None): Why a MitoFish seed couldn't be used, when seeded is False and a real reason was recorded (never fabricated).
    """
    organelle_type: str
    seeded: bool
    seed_file: str | None
    seed_gene: str | None = None
    seed_accession: str | None = None
    seed_source: str | None = None
    seed_fallback_reason: str | None = None

@dataclass(frozen = True, slots = True)
class RescueAttemptSummary:
    """
    One real rescue attempt's own outcome for this sample -- read
    from #07's own real rescue_mapping.json (whether the attempt even
    ran) and, if it did, the reassembly's own real assembly_status.json
    (what GetOrganelle actually produced). Recorded whether or not the
    attempt was ultimately selected, so the methods section can report
    every real attempt made, not only the one that succeeded.

    Attributes:
        attempt (int): The rescue attempt number (1-indexed).
        attempted (bool): Whether this attempt actually ran a reassembly, or was skipped because its own trigger condition wasn't met.
        reason (str): The real reason recorded for whether/why this attempt ran.
        profile_name (str | None): The GetOrganelle profile name used for this attempt's reassembly, or None if it never ran.
        status (str | None): The reassembly's own real outcome status (e.g. "circular", "fragmented", "missing"), or None if it never ran.
    """
    attempt: int
    attempted: bool
    reason: str
    profile_name: str | None = None
    status: str | None = None

@dataclass(frozen = True, slots = True)
class BlastConfigSummary:
    """The real, deterministic BLAST settings used for this sample."""
    database: str
    mode: str

@dataclass(frozen = True, slots = True)
class AnnotationConfigSummary:
    """The real, deterministic MITOS2 settings used for this sample."""
    genetic_code: int
    refseqver: str

@dataclass(frozen = True, slots = True)
class ArchivalSummary:
    """
    One sample's real NCBI archival preparation status (#367-376),
    read from the same SRA/BioSample/GenBank tables a human reviewer
    would open before uploading through the NCBI Submission Portal --
    this pipeline never submits directly.

    Attributes:
        sra_review_status (str | None): The SRA record's own review_status (READY/MANUAL_REVIEW_REQUIRED), or None if no record exists.
        sra_missing_fields (str): The SRA record's own missing_fields, semicolon-delimited.
        biosample_review_status (str | None): The BioSample record's own review_status, or None if no record exists.
        biosample_missing_fields (str): The BioSample record's own missing_fields, semicolon-delimited.
        genbank_review_status (str | None): The GenBank submission folder's own review_status, or None if no record exists.
        genbank_missing (str): The GenBank record's own missing_files + missing_metadata, semicolon-delimited.
        ready_for_review (bool): Whether every real archival source is READY for this sample.
    """
    sra_review_status: str | None = None
    sra_missing_fields: str = ""
    biosample_review_status: str | None = None
    biosample_missing_fields: str = ""
    genbank_review_status: str | None = None
    genbank_missing: str = ""

    @property
    def ready_for_review(self) -> bool:
        statuses = (self.sra_review_status, self.biosample_review_status, self.genbank_review_status)
        return all(s == "READY" for s in statuses)

@dataclass(frozen = True, slots = True)
class ReadQCSummary:
    """
    One read-set's own real FastQC-derived QC summary -- raw, trimmed,
    or recruited reads -- combining R1 and R2 into one display-ready
    record read straight from FastQC's own real zip output, never
    recomputed from the raw reads themselves.

    Attributes:
        total_sequences (int | None): R1 + R2 total sequence count.
        percent_gc (float | None): The average of R1's and R2's own real %GC.
        sequence_length (str | None): FastQC's own real sequence-length value (R1/R2 always match in practice).
        overall_verdict (str): The worse of R1's and R2's own overall FastQC verdict (PASS/WARN/FAIL).
        representative_image_png (Path | None): R1's own real per-base-quality chart, extracted from its FastQC zip, or None if unavailable.
    """
    total_sequences: int | None
    percent_gc: float | None
    sequence_length: str | None
    overall_verdict: str
    representative_image_png: Path | None = None

@dataclass(frozen = True, slots = True)
class SampleReportData:
    """
    Everything #358's per-sample report (and the live methods
    section) needs for one sample, assembled from real result JSON.

    Attributes:
        sample_id (str): The sample id.
        metadata (dict[str, str]): This sample's row from the runtime manifest.
        fastp_stats (dict | None): Parsed fastp summary stats, or None when unavailable.
        assembly (AssemblyStatistics | None): Assembly statistics, or None.
        blast (BlastStatistics | None): BLAST statistics, or None.
        annotation (AnnotationStatistics | None): Annotation statistics, or None.
        overall_verdict (str | None): The biological validation gate's overall verdict (PASS/REVIEW/FAIL) for the selected attempt, or None.
        validation_checks (tuple[dict, ...]): Every individual validation check's own result for the selected attempt.
        visualization_png (Path | None): The rendered circular genome map PNG, or None.
        assembly_config (AssemblyConfigSummary | None): The real GetOrganelle settings used, or None.
        blast_config (BlastConfigSummary | None): The real BLAST settings used, or None.
        annotation_config (AnnotationConfigSummary | None): The real MITOS2 settings used, or None.
        archival (ArchivalSummary | None): This sample's real NCBI archival preparation status, or None when the archival stages didn't run.
        phylogeny (PhylogenySummary | None): This sample's own real phylogenetic outcome (its assembled genome plus its own live-fetched BLAST references, restoring v0.1.0's per-sample design), or None when phylogeny didn't run for this sample.
        raw_read_qc (ReadQCSummary | None): FastQC's own real QC summary for the raw sequencing reads, or None when unavailable.
        trimmed_read_qc (ReadQCSummary | None): FastQC's own real QC summary for the trimmed reads, or None when unavailable.
        recruited_read_qc (ReadQCSummary | None): FastQC's own real QC summary for the reads GetOrganelle recruited as mitochondrial-like before assembly, or None when unavailable.
        aggregated_data_dir (Path | None): The run's own consolidated aggregated_data/ directory, or None when unavailable.
        rescue_attempts (tuple[RescueAttemptSummary, ...]): Every real rescue attempt made for this sample, in attempt order -- empty when rescue was never configured or never triggered.
        rescued (bool): Whether the selected attempt came from a rescue retry (source_attempt > 0), not the primary assembly.
    """
    sample_id: str
    metadata: dict[str, str]
    fastp_stats: dict | None
    assembly: AssemblyStatistics | None
    blast: BlastStatistics | None
    annotation: AnnotationStatistics | None
    overall_verdict: str | None
    validation_checks: tuple[dict, ...] = field(default_factory = tuple)
    visualization_png: Path | None = None
    assembly_config: AssemblyConfigSummary | None = None
    blast_config: BlastConfigSummary | None = None
    annotation_config: AnnotationConfigSummary | None = None
    archival: ArchivalSummary | None = None
    phylogeny: "PhylogenySummary | None" = None
    raw_read_qc: ReadQCSummary | None = None
    trimmed_read_qc: ReadQCSummary | None = None
    recruited_read_qc: ReadQCSummary | None = None
    aggregated_data_dir: Path | None = None
    rescue_attempts: tuple[RescueAttemptSummary, ...] = field(default_factory = tuple)

    @property
    def rescued(self) -> bool:
        return self.assembly is not None and self.assembly.source_attempt > 0

@dataclass(frozen = True, slots = True)
class PhylogenySummary:
    """
    One sample's own real phylogenetic outcome: its assembled genome
    plus its own live-fetched BLAST reference sequences, aligned and
    run through IQ-TREE (restoring v0.1.0's per-sample design -- a
    sample's own tree never depends on how many other samples are in
    the same run).

    Attributes:
        success (bool): Whether IQ-TREE produced a usable tree.
        tip_count (int | None): The tree's tip count, or None if no tree was built.
        newick_path (Path | None): The real Newick tree path, or None.
        reason (str): A short, human-readable explanation of the outcome.
        model (str | None): The substitution model IQ-TREE was run with, when known.
        tree_figure_png (Path | None): The rendered tree figure PNG, or None if rendering didn't run/succeed.
    """
    success: bool
    tip_count: int | None
    newick_path: Path | None
    reason: str
    model: str | None = None
    tree_figure_png: Path | None = None

@dataclass(frozen = True, slots = True)
class RunReportData:
    """
    Everything #360's run-level report needs -- every sample's own
    report data (each one carrying its own real phylogenetic outcome,
    restoring v0.1.0's per-sample design; there is no batch-level
    phylogeny artifact anymore).

    Attributes:
        job_id (str): The pipeline job id.
        generated_at (str): An ISO-8601 UTC timestamp for when this report was generated.
        samples (tuple[SampleReportData, ...]): Every sample's own report data, in manifest order.
    """
    job_id: str
    generated_at: str
    samples: tuple[SampleReportData, ...]

    @property
    def verdict_counts(self) -> dict[str, int]:
        """Counts samples by overall_verdict (PASS/REVIEW/FAIL/"not assessed")."""
        counts: dict[str, int] = {}
        for sample in self.samples:
            key = sample.overall_verdict or "not assessed"
            counts[key] = counts.get(key, 0) + 1
        return counts

    @property
    def archival_ready_count(self) -> int:
        """How many samples have every real archival source (SRA/BioSample/GenBank) at READY."""
        return sum(1 for sample in self.samples if sample.archival is not None and sample.archival.ready_for_review)
