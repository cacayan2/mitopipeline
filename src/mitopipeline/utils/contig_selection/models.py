"""
models.py

Data models for selecting a sample's final mitochondrial contig across
one or more GetOrganelle assembly attempts (the primary attempt, plus
any rescue attempts).
"""

# Imports
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path

@dataclass(frozen = True, slots = True)
class ContigCandidate:
    """
    One assembly attempt's result, as discovered by
    getorganelle.discover_assembly_artifacts(), labeled with which
    attempt it came from.

    Attributes:
        attempt (int): The attempt number this candidate came from (0 = primary, 1+ = rescue attempts, in the order they were tried).
        status (str): "circular", "fragmented", or "missing", as discover_assembly_artifacts() classified it.
        fasta_path (Path | None): The candidate's native FASTA, or None when status is "missing".
        gfa_path (Path | None): The candidate's native GFA, or None when status is "missing".
    """
    attempt: int
    status: str
    fasta_path: Path | None
    gfa_path: Path | None

@dataclass(frozen = True, slots = True)
class ContigSelection:
    """
    The single contig selected as a sample's final mitochondrial
    assembly, across every attempt considered.

    Attributes:
        status (str): The selected candidate's status ("circular" or "fragmented") -- never "missing", since a missing candidate is never selectable (see select_mitochondrial_contig).
        fasta_path (Path): The selected candidate's native FASTA.
        gfa_path (Path | None): The selected candidate's native GFA, if it has one.
        source_attempt (int): Which attempt the selection came from.
        length_bp (int): The selected sequence's length in base pairs, as read from fasta_path. The primary ranking criterion among same-status candidates.
        reason (str): A short, human-readable explanation of why this candidate won, for logging.
    """
    status: str
    fasta_path: Path
    gfa_path: Path | None
    source_attempt: int
    length_bp: int
    reason: str
