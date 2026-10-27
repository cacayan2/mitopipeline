"""
models.py

Typed, per-sample statistics records (#353-356) -- each extracted
from the same real result JSON a prior stage already wrote, never
re-derived from raw tool output a second time.
"""

# Imports
from dataclasses import asdict, dataclass, field

@dataclass(frozen = True, slots = True)
class AssemblyStatistics:
    """
    Statistics for one sample's selected assembly (#353).

    Attributes:
        sample_id (str): The sample id.
        status (str): "circular" or "fragmented" -- contig_selection's own recorded status.
        source_attempt (int): Which assembly attempt was selected (0 = primary, 1+ = rescue).
        length_bp (int): The selected FASTA's total length, in base pairs.
        contig_count (int): The number of sequence records in the selected FASTA.
        gc_content_percent (float): The selected FASTA's GC content, as a percentage.
        recruited_read_count (int | None): Total recruited reads (R1 + R2) for the selected attempt, or None when unavailable.
        runtime_seconds (float | None): The selected attempt's own assembly runtime, or None when unavailable.
    """
    sample_id: str
    status: str
    source_attempt: int
    length_bp: int
    contig_count: int
    gc_content_percent: float
    recruited_read_count: int | None = None
    runtime_seconds: float | None = None

    def to_dict(self) -> dict:
        return asdict(self)

@dataclass(frozen = True, slots = True)
class BlastStatistics:
    """
    Statistics for one sample's top BLAST hit (#354).

    Attributes:
        sample_id (str): The sample id.
        top_accession (str): The top hit's subject sequence id.
        top_scientific_name (str): The top hit's subject scientific name(s), as BLAST reports them.
        top_identity_percent (float): The top hit's percentage identity.
        top_query_coverage_percent (float): The top hit's query coverage, as a percentage.
        top_evalue (float): The top hit's expect value.
        decision_status (str | None): The biological validation gate's own "blast_identity" verdict (PASS/REVIEW/FAIL) for the selected attempt, or None when unavailable.
    """
    sample_id: str
    top_accession: str
    top_scientific_name: str
    top_identity_percent: float
    top_query_coverage_percent: float
    top_evalue: float
    decision_status: str | None = None

    def to_dict(self) -> dict:
        return asdict(self)

@dataclass(frozen = True, slots = True)
class AnnotationStatistics:
    """
    Statistics for one sample's MITOS2 annotation (#355).

    Attributes:
        sample_id (str): The sample id.
        gene_count (int): The total number of gene features MITOS2 found.
        core_gene_count (int): How many of the 13 core animal mitochondrial protein-coding genes were found.
        missing_core_genes (tuple[str, ...]): Which core genes were NOT found, sorted.
    """
    sample_id: str
    gene_count: int
    core_gene_count: int
    missing_core_genes: tuple[str, ...] = field(default_factory = tuple)

    def to_dict(self) -> dict:
        return {**asdict(self), "missing_core_genes": list(self.missing_core_genes)}

@dataclass(frozen = True, slots = True)
class SampleStatistics:
    """
    One sample's merged statistics record (#356), joined by sample_id.
    Any stage that didn't run, or found nothing usable, is None here
    -- never a fabricated zero standing in for "never measured."

    Attributes:
        sample_id (str): The sample id.
        assembly (AssemblyStatistics | None): Assembly statistics, or None if no usable assembly was selected.
        blast (BlastStatistics | None): BLAST statistics, or None if BLAST found no hits or didn't run.
        annotation (AnnotationStatistics | None): Annotation statistics, or None if MITOS2 didn't run or produced nothing.
    """
    sample_id: str
    assembly: AssemblyStatistics | None
    blast: BlastStatistics | None
    annotation: AnnotationStatistics | None

    def to_dict(self) -> dict:
        return {
            "sample_id": self.sample_id,
            "assembly": self.assembly.to_dict() if self.assembly else None,
            "blast": self.blast.to_dict() if self.blast else None,
            "annotation": self.annotation.to_dict() if self.annotation else None,
        }
