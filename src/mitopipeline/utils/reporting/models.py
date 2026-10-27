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
    """
    organelle_type: str
    seeded: bool
    seed_file: str | None

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

    @property
    def rescued(self) -> bool:
        return self.assembly is not None and self.assembly.source_attempt > 0

@dataclass(frozen = True, slots = True)
class PhylogenySummary:
    """
    The batch-level phylogeny outcome (#344-348's own, job-level
    design -- not per-sample).

    Attributes:
        success (bool): Whether IQ-TREE produced a usable tree.
        tip_count (int | None): The tree's tip count, or None if no tree was built.
        newick_path (Path | None): The real Newick tree path, or None.
        reason (str): A short, human-readable explanation of the outcome.
        model (str | None): The substitution model IQ-TREE was run with, when known.
    """
    success: bool
    tip_count: int | None
    newick_path: Path | None
    reason: str
    model: str | None = None

@dataclass(frozen = True, slots = True)
class RunReportData:
    """
    Everything #360's run-level report needs -- every sample's own
    report data plus the batch-level phylogeny outcome.

    Attributes:
        job_id (str): The pipeline job id.
        generated_at (str): An ISO-8601 UTC timestamp for when this report was generated.
        samples (tuple[SampleReportData, ...]): Every sample's own report data, in manifest order.
        phylogeny (PhylogenySummary | None): The batch-level phylogeny outcome, or None if phylogeny didn't run.
    """
    job_id: str
    generated_at: str
    samples: tuple[SampleReportData, ...]
    phylogeny: PhylogenySummary | None = None

    @property
    def verdict_counts(self) -> dict[str, int]:
        """Counts samples by overall_verdict (PASS/REVIEW/FAIL/"not assessed")."""
        counts: dict[str, int] = {}
        for sample in self.samples:
            key = sample.overall_verdict or "not assessed"
            counts[key] = counts.get(key, 0) + 1
        return counts
