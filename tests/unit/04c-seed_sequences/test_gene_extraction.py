"""
test_gene_extraction.py

Unit tests for mitochondrial gene extraction from a reference record.
"""

# Imports
import logging
from pathlib import Path
import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord
from mitopipeline.utils.seed.gene_extraction import extract_gene, write_gene_fasta

# A synthetic 40bp reference. Positions 5-15 hold the (forward-strand)
# 12S locus; positions 25-35 hold the (reverse-strand) 16S locus.
_TWELVE_S_SEQ = "AAAACCCCGG"
_SIXTEEN_S_SEQ_FWD = "TTTTGGGGCC"

def _reference_record() -> SeqRecord:
    """
    Builds a synthetic reference record with one forward-strand and one
    reverse-strand gene feature, using real GenBank-style qualifiers.

    Returns:
        SeqRecord: The synthetic reference record.
    """
    full_seq = "NNNNN" + _TWELVE_S_SEQ + "NNNNNNNNNN" + _SIXTEEN_S_SEQ_FWD + "NNNNN"
    record = SeqRecord(Seq(full_seq), id = "REF_001")
    record.features = [
        SeqFeature(
            SimpleLocation(5, 15, strand = 1),
            type = "rRNA",
            qualifiers = {"product": ["12S ribosomal RNA"]},
        ),
        SeqFeature(
            SimpleLocation(25, 35, strand = -1),
            type = "rRNA",
            qualifiers = {"product": ["16S rRNA"]},
        ),
    ]
    return record

def test_extract_gene_forward_strand_returns_correct_sequence() -> None:
    """
    A forward-strand feature should be extracted exactly as stored.
    """
    result = extract_gene(_reference_record(), "12S_rRNA", logging.getLogger("test"))
    assert str(result.seq) == _TWELVE_S_SEQ

def test_extract_gene_matches_product_synonym() -> None:
    """
    A synonym for the canonical gene name ("12S ribosomal RNA") should
    resolve to the same locus as the canonical name itself.
    """
    canonical_result = extract_gene(_reference_record(), "12S_rRNA", logging.getLogger("test"))
    synonym_result = extract_gene(_reference_record(), "12s ribosomal rna", logging.getLogger("test"))
    assert str(canonical_result.seq) == str(synonym_result.seq)

def test_extract_gene_reverse_strand_returns_reverse_complement() -> None:
    """
    A minus-strand feature should be returned reverse-complemented, not
    as the raw forward-strand slice.
    """
    result = extract_gene(_reference_record(), "16S_rRNA", logging.getLogger("test"))
    expected = str(Seq(_SIXTEEN_S_SEQ_FWD).reverse_complement())
    assert str(result.seq) == expected
    assert str(result.seq) != _SIXTEEN_S_SEQ_FWD

def test_extract_gene_exits_on_unrecognized_gene_name() -> None:
    """
    A gene name with no known canonical mapping should exit rather than
    search the record at all.
    """
    with pytest.raises(SystemExit):
        extract_gene(_reference_record(), "not_a_real_gene", logging.getLogger("test"))

def test_extract_gene_exits_on_missing_locus() -> None:
    """
    A recognized gene name absent from this particular reference record
    (missing locus) should exit rather than return an empty/wrong result.
    """
    record = _reference_record()
    record.features = [record.features[0]]  # only the 12S feature remains
    with pytest.raises(SystemExit):
        extract_gene(record, "16S_rRNA", logging.getLogger("test"))

def test_write_gene_fasta_creates_readable_file(tmp_path: Path) -> None:
    """
    write_gene_fasta should create parent directories as needed and write
    a FASTA file whose content round-trips back to the same sequence.
    """
    from Bio import SeqIO
    seq_record = extract_gene(_reference_record(), "12S_rRNA", logging.getLogger("test"))
    output_path = tmp_path / "nested" / "REF_001_12S.fasta"
    written_path = write_gene_fasta(seq_record, output_path)
    assert written_path == output_path
    assert output_path.is_file()
    round_tripped = next(SeqIO.parse(output_path, "fasta"))
    assert str(round_tripped.seq) == _TWELVE_S_SEQ
