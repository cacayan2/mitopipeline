"""
models.py

Data models for mitochondrial seed-sequence selection.
"""

# Imports
from dataclasses import dataclass

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
