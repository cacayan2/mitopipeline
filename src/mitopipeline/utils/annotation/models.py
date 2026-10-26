"""
models.py

Data models for parsed MITOS2 gene features.
"""

# Imports
from dataclasses import dataclass, field

@dataclass(frozen = True, slots = True)
class GeneFeature:
    """
    One normalized gene feature from MITOS2's own GFF3 output
    (result.gff).

    Attributes:
        seqid (str): The sequence (contig/genome) this feature is on.
        gene_id (str): The gene's name/id (MITOS2's "Name" qualifier, e.g. "cox2").
        start (int): 1-based start coordinate.
        end (int): 1-based end coordinate.
        strand (str | None): "+" or "-", or None when MITOS2 reports "." (unknown).
        score (float | None): The feature's score, or None when MITOS2 reports "." (no score).
        qualifiers (dict[str, str]): The raw GFF attribute key/value pairs (e.g. ID, Name, gene_id), unmodified.
    """
    seqid: str
    gene_id: str
    start: int
    end: int
    strand: str | None
    score: float | None
    qualifiers: dict[str, str] = field(default_factory = dict)
