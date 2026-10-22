"""
selector.py

Selects a sample's single final mitochondrial contig across one or
more GetOrganelle assembly attempts.

Scope note: this ranks between attempts using only what
discover_assembly_artifacts() already knows (GetOrganelle's own
circular/fragmented classification) plus sequence length -- it does
NOT yet incorporate independent BLAST evidence the way v0.1.0's
assess_mitochondrial_assembly() did (taxon-aware length bounds,
per-contig BLAST support). That richer, evidence-based refinement
needs the BLAST reference-retrieval and exec stages (#401, #405a-c)
wired first, which aren't built yet in v0.2.0 -- tracked as a known
follow-on, not faked here.
"""

# Imports
from __future__ import annotations
from logging import Logger
from pathlib import Path
from Bio import SeqIO
from mitopipeline.utils.contig_selection.models import ContigCandidate, ContigSelection

# Lower is better -- "circular" always outranks "fragmented" regardless
# of length. "missing" never reaches this map: it's filtered out before
# ranking (see select_mitochondrial_contig).
_STATUS_PRIORITY = {"circular": 0, "fragmented": 1}

def _sequence_length(fasta_path: Path) -> int:
    """
    Returns the total base count across every record in a FASTA file.

    Args:
        fasta_path (Path): The FASTA file to measure.

    Returns:
        int: The total sequence length, in base pairs.
    """
    return sum(len(record.seq) for record in SeqIO.parse(fasta_path, "fasta"))

def select_mitochondrial_contig(candidates: list[ContigCandidate],
                                 logger: Logger | None = None) -> ContigSelection | None:
    """
    Selects the single best mitochondrial contig across every assembly
    attempt considered for one sample, using explicit, deterministic
    criteria: a circular/complete result always beats a fragmented one
    regardless of attempt or length; among equally-classified
    candidates, the longer sequence wins; a tie beyond that is broken
    by preferring the earlier attempt (the primary assembly over any
    later rescue attempt), so the result never depends on dict/set
    ordering or attempt-list order beyond what's documented here.

    A "missing" candidate (no native assembly output at all) is never
    selectable -- it carries no fasta_path to select, and candidates
    is expected to describe real attempts, not placeholders for ones
    that produced nothing.

    Args:
        candidates (list[ContigCandidate]): Every assembly attempt's result considered for this sample, in any order.
        logger (Logger | None): Logger to use for the decision. Defaults to None.

    Returns:
        ContigSelection | None: The selected contig, or None if every candidate was "missing" (nothing to select) or candidates was empty.
    """
    selectable = [c for c in candidates if c.status in _STATUS_PRIORITY]

    if not selectable:
        if logger is not None:
            logger.warning("No selectable (circular/fragmented) candidate "
                            f"among {len(candidates)} attempt(s); nothing to select.")
        return None

    lengths = {c.attempt: _sequence_length(c.fasta_path) for c in selectable}

    best = min(selectable, key = lambda c: (
        _STATUS_PRIORITY[c.status], -lengths[c.attempt], c.attempt))

    runner_up_note = ""
    if len(selectable) > 1:
        runner_up_note = f" ({len(selectable) - 1} other candidate(s) considered and not chosen)"

    reason = (f"attempt {best.attempt}: {best.status}, {lengths[best.attempt]}bp "
              f"-- best by status then length then earliest attempt{runner_up_note}")

    if logger is not None:
        logger.info(reason)

    return ContigSelection(
        status = best.status,
        fasta_path = best.fasta_path,
        gfa_path = best.gfa_path,
        source_attempt = best.attempt,
        length_bp = lengths[best.attempt],
        reason = reason,
    )
