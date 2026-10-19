"""
models.py

Data models for parsed and ranked BLAST hits.
"""

# Imports
from dataclasses import dataclass

@dataclass(frozen = True, slots = True)
class BlastHit:
    """
    One row of tabular BLAST output (-outfmt 6, with the full explicit
    column list blast.BLAST_OUTPUT_COLUMNS requests).

    Attributes:
        qseqid (str): The query sequence ID.
        sseqid (str): The subject (database hit) sequence ID.
        pident (float): Percentage of identical matches.
        length (int): Alignment length.
        mismatch (int): Number of mismatches.
        gapopen (int): Number of gap openings.
        qstart (int): Start of alignment in the query.
        qend (int): End of alignment in the query.
        sstart (int): Start of alignment in the subject.
        send (int): End of alignment in the subject.
        evalue (float): Expect value -- lower is a stronger match.
        bitscore (float): Bit score -- higher is a stronger match.
        qcovs (float): Query coverage per subject, as a percentage.
        staxids (str): The subject's NCBI taxonomy ID(s), as BLAST reports them.
        sscinames (str): The subject's scientific name(s), as BLAST reports them.
        stitle (str): The subject's sequence title.
        rank (int | None): This hit's rank among the sample's selected top hits, or None before selection.
    """
    qseqid: str
    sseqid: str
    pident: float
    length: int
    mismatch: int
    gapopen: int
    qstart: int
    qend: int
    sstart: int
    send: int
    evalue: float
    bitscore: float
    qcovs: float
    staxids: str
    sscinames: str
    stitle: str
    rank: int | None = None
