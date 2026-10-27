"""
test_circular_coverage.py

Unit tests for circular-genome BLAST coverage correction (#414):
doubling a circular query so an origin-spanning match reads as one
unbroken hit, then measuring real coverage against the true,
undoubled genome length.
"""

# Imports
from pathlib import Path
import pytest
from mitopipeline.utils.blast.circular_coverage import compute_circular_query_coverage, double_fasta_sequence
from mitopipeline.utils.blast.models import BlastHit

def _hit(**overrides) -> BlastHit:
    defaults = dict(
        qseqid = "q", sseqid = "ref1", pident = 100.0, length = 16596, mismatch = 0, gapopen = 0,
        qstart = 1, qend = 16596, sstart = 1, send = 16596, evalue = 0.0, bitscore = 30000.0,
        qcovs = 100.0, staxids = "1", sscinames = "Example species", stitle = "Example", rank = None,
    )
    defaults.update(overrides)
    return BlastHit(**defaults)

def test_double_fasta_sequence_doubles_the_real_sequence(tmp_path: Path) -> None:
    fasta = tmp_path / "data.fasta"
    fasta.write_text(">sample1 (circular)\nACGT\n")
    output_path = tmp_path / "doubled.fasta"

    genome_length = double_fasta_sequence(fasta, output_path)

    assert genome_length == 4
    text = output_path.read_text()
    assert text == ">sample1 (circular)\nACGTACGT\n"

def test_double_fasta_sequence_handles_multiline_sequence(tmp_path: Path) -> None:
    fasta = tmp_path / "data.fasta"
    fasta.write_text(">s\nACGT\nTTTT\n")
    output_path = tmp_path / "doubled.fasta"

    genome_length = double_fasta_sequence(fasta, output_path)

    assert genome_length == 8
    assert output_path.read_text() == ">s\nACGTTTTTACGTTTTT\n"

def test_double_fasta_sequence_missing_file_raises(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        double_fasta_sequence(tmp_path / "missing.fasta", tmp_path / "out.fasta")

def test_double_fasta_sequence_empty_file_raises(tmp_path: Path) -> None:
    fasta = tmp_path / "data.fasta"
    fasta.write_text("")
    with pytest.raises(ValueError):
        double_fasta_sequence(fasta, tmp_path / "out.fasta")

def test_compute_circular_query_coverage_full_length_match() -> None:
    # The real-world case this fixes: an origin-spanning match that
    # only shows as a single-HSP 53% hit on the undoubled query now
    # reads as one unbroken, full-length hit on the doubled query.
    hits = [_hit(length = 16596, qstart = 7740, qend = 24335)]
    result = compute_circular_query_coverage(hits, genome_length = 16596)
    assert result == (100.0, 16596)

def test_compute_circular_query_coverage_caps_at_100_percent() -> None:
    # A highly repetitive genome could in principle report an
    # alignment longer than the real genome itself once doubled --
    # coverage must never exceed 100%.
    hits = [_hit(length = 20000)]
    coverage, alignment_length = compute_circular_query_coverage(hits, genome_length = 16596)
    assert coverage == 100.0
    assert alignment_length == 20000

def test_compute_circular_query_coverage_partial_match() -> None:
    hits = [_hit(length = 8298)]
    coverage, _ = compute_circular_query_coverage(hits, genome_length = 16596)
    assert coverage == pytest.approx(50.0)

def test_compute_circular_query_coverage_picks_the_strongest_hit() -> None:
    weak = _hit(sseqid = "ref_weak", bitscore = 100.0, evalue = 1e-5, length = 100)
    strong = _hit(sseqid = "ref_strong", bitscore = 30000.0, evalue = 0.0, length = 16596)
    coverage, alignment_length = compute_circular_query_coverage([weak, strong], genome_length = 16596)
    assert alignment_length == 16596
    assert coverage == 100.0

def test_compute_circular_query_coverage_none_when_no_hits() -> None:
    assert compute_circular_query_coverage([], genome_length = 16596) is None

def test_compute_circular_query_coverage_rejects_invalid_genome_length() -> None:
    with pytest.raises(ValueError):
        compute_circular_query_coverage([_hit()], genome_length = 0)
