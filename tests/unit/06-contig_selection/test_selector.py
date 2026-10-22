"""
test_selector.py

Tests for select_mitochondrial_contig() -- ranking a sample's final
mitochondrial contig across one or more GetOrganelle assembly attempts.
"""

# Imports
from pathlib import Path
import pytest
from mitopipeline.utils.contig_selection.models import ContigCandidate
from mitopipeline.utils.contig_selection.selector import select_mitochondrial_contig

def _write_fasta(path: Path, *, header: str = "seq", length: int) -> Path:
    """Writes a single-record FASTA of exactly `length` bases."""
    path.write_text(f">{header}\n{'A' * length}\n")
    return path

def test_no_candidates_returns_none():
    assert select_mitochondrial_contig([]) is None

def test_all_missing_returns_none(tmp_path: Path):
    candidates = [
        ContigCandidate(attempt = 0, status = "missing", fasta_path = None, gfa_path = None),
        ContigCandidate(attempt = 1, status = "missing", fasta_path = None, gfa_path = None),
    ]
    assert select_mitochondrial_contig(candidates) is None

def test_circular_always_beats_fragmented_regardless_of_length(tmp_path: Path):
    # The fragmented candidate is deliberately much longer -- status
    # priority must still win over raw length.
    circular_fasta = _write_fasta(tmp_path / "circular.fasta", length = 100)
    fragmented_fasta = _write_fasta(tmp_path / "fragmented.fasta", length = 50_000)
    candidates = [
        ContigCandidate(attempt = 0, status = "fragmented", fasta_path = fragmented_fasta, gfa_path = None),
        ContigCandidate(attempt = 1, status = "circular", fasta_path = circular_fasta, gfa_path = None),
    ]
    selected = select_mitochondrial_contig(candidates)
    assert selected is not None
    assert selected.status == "circular"
    assert selected.source_attempt == 1
    assert selected.length_bp == 100

def test_longer_sequence_wins_among_same_status(tmp_path: Path):
    shorter = _write_fasta(tmp_path / "shorter.fasta", length = 1_000)
    longer = _write_fasta(tmp_path / "longer.fasta", length = 15_000)
    candidates = [
        ContigCandidate(attempt = 0, status = "fragmented", fasta_path = shorter, gfa_path = None),
        ContigCandidate(attempt = 1, status = "fragmented", fasta_path = longer, gfa_path = None),
    ]
    selected = select_mitochondrial_contig(candidates)
    assert selected is not None
    assert selected.fasta_path == longer
    assert selected.length_bp == 15_000

def test_exact_tie_prefers_earlier_attempt(tmp_path: Path):
    primary = _write_fasta(tmp_path / "primary.fasta", length = 16_000)
    rescue = _write_fasta(tmp_path / "rescue.fasta", length = 16_000)
    candidates = [
        ContigCandidate(attempt = 2, status = "fragmented", fasta_path = rescue, gfa_path = None),
        ContigCandidate(attempt = 0, status = "fragmented", fasta_path = primary, gfa_path = None),
    ]
    selected = select_mitochondrial_contig(candidates)
    assert selected is not None
    assert selected.source_attempt == 0
    assert selected.fasta_path == primary

def test_missing_candidates_are_ignored_not_selected(tmp_path: Path):
    fragmented = _write_fasta(tmp_path / "fragmented.fasta", length = 12_000)
    candidates = [
        ContigCandidate(attempt = 0, status = "missing", fasta_path = None, gfa_path = None),
        ContigCandidate(attempt = 1, status = "fragmented", fasta_path = fragmented, gfa_path = None),
    ]
    selected = select_mitochondrial_contig(candidates)
    assert selected is not None
    assert selected.source_attempt == 1
    assert selected.status == "fragmented"

def test_single_circular_candidate_selected_directly(tmp_path: Path):
    fasta = _write_fasta(tmp_path / "only.fasta", length = 16_596)
    candidates = [ContigCandidate(attempt = 0, status = "circular", fasta_path = fasta, gfa_path = None)]
    selected = select_mitochondrial_contig(candidates)
    assert selected is not None
    assert selected.length_bp == 16_596
    assert selected.status == "circular"
    assert selected.source_attempt == 0

def test_gfa_path_carried_through_from_the_winning_candidate(tmp_path: Path):
    fasta = _write_fasta(tmp_path / "winner.fasta", length = 16_000)
    gfa = tmp_path / "winner.gfa"
    gfa.write_text("H\tVN:Z:1.0\n")
    candidates = [ContigCandidate(attempt = 0, status = "circular", fasta_path = fasta, gfa_path = gfa)]
    selected = select_mitochondrial_contig(candidates)
    assert selected is not None
    assert selected.gfa_path == gfa

def test_logger_receives_a_decision_message(tmp_path: Path, caplog):
    import logging
    logger = logging.getLogger("test-contig-selection")
    fasta = _write_fasta(tmp_path / "only.fasta", length = 500)
    candidates = [ContigCandidate(attempt = 0, status = "fragmented", fasta_path = fasta, gfa_path = None)]
    with caplog.at_level(logging.INFO, logger = "test-contig-selection"):
        select_mitochondrial_contig(candidates, logger = logger)
    assert any("attempt 0" in record.message for record in caplog.records)
