"""
models.py

Data models for mitochondrial seed-sequence selection.
"""

# Imports
from dataclasses import dataclass
from pathlib import Path

# The default mitochondrial locus to search for when a sample's config does
# not request a specific one. Configurable per call site (lookup(), seed
# selection orchestration), never hard-coded past this single constant.
DEFAULT_SEED_LOCUS = "12S_rRNA"

# The order select_seed() tries loci in, most densely documented
# first, so a sample with no usable 12S_rRNA annotation still gets a
# real seed instead of going straight to GetOrganelle's default
# database. Derived from REAL per-gene record counts in the local
# MitoFish database (tmp_playground/mitofish_tables/seq_annotation.parquet,
# measured 2026-10-19), not assumed -- COXI (the standard DNA-barcoding
# locus) turns out to be far better documented than 12S_rRNA, which
# used to be the only locus select_seed() ever tried:
#   COXI 431865, Cytb 249281, 16S_rRNA 77708, 12S_rRNA 68097,
#   ND2 46962, ATPase6 37921, ND4 26733, ND1 24582, ND5 24239,
#   ND6 21669, COXIII 21600, ND4L 20727, ND3 20716, COXII 20616.
LOCUS_FALLBACK_ORDER = (
    "COXI", "Cytb", "16S_rRNA", "12S_rRNA", "ND2", "ATPase6", "ND4",
    "ND1", "ND5", "ND6", "COXIII", "ND4L", "ND3", "COXII",
)

@dataclass(frozen = True, slots = True)
class SeedReferenceCandidate:
    """
    Candidate mitochondrial reference returned by a reference database.

    Attributes:
        accession (str): The source database accession for the mitochondrial sequence.
        scientific_name (str): The scientific name associated with the reference sequence.
        ncbi_taxon_id (int): The NCBI taxonomy identifier associated with the organism.
        gene (str): Mitochonddrial gene available from the reference.
        sequence_length (int): The length of the candidate gene sequence, if known.
    """
    accession: str
    scientific_name: str
    ncbi_taxon_id: int
    gene: str
    sequence_length: int | None = None

@dataclass(frozen = True, slots = True)
class SeedSelectionResult:
    """
    The outcome of selecting and writing a seed reference sequence for one
    sample -- what locus was used, where the sequence came from, and
    whether a fallback database or locus was needed to find it.

    Attributes:
        sample_id (str): The sample this seed selection was performed for.
        scientific_name (str): The scientific name the selection was queried against.
        gene (str): The mitochondrial locus actually used (may differ from the requested locus if a fallback was needed).
        accession (str): The source database accession for the selected reference sequence, or a fallback marker (e.g. "GETORGANELLE_DEFAULT") when no seed was selected.
        source (str): Which reference source produced this result (e.g. "mitofish" or "fallback_default_database").
        fasta_path (Path | None): The path to the written seed FASTA file, or None when no seed was written (fallback_used is True).
        fallback_used (bool): Whether the configured database/locus could not be used, and a fallback was substituted.
        fallback_reason (str | None): Why a fallback was needed (e.g. "no_mitofish_candidate", "ambiguous_candidates", "mitofish_lookup_failed"), or None when fallback_used is False.
        match_level (str | None): How confidently the selected candidate matched the requested organism/locus (e.g. "exact"), or None for a fallback result.
        checksum (str | None): A SHA-256 checksum of the selected sequence, for provenance/reproducibility -- None for a fallback result with no sequence.
        ncbi_taxon_id (int | None): The NCBI taxonomy identifier associated with the organism, if known -- a fallback source may not have one.
        sequence_length (int | None): The length of the selected sequence, if known.
    """
    sample_id: str
    scientific_name: str
    gene: str
    accession: str
    source: str
    fasta_path: Path | None = None
    fallback_used: bool = False
    fallback_reason: str | None = None
    match_level: str | None = None
    checksum: str | None = None
    ncbi_taxon_id: int | None = None
    sequence_length: int | None = None
