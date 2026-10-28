"""
test_sample_dataset.py

Unit tests for one sample's own real phylogenetic dataset builder
(restoring v0.1.0's per-sample design, ported here): the sample's own
assembled genome plus its own top-ranked BLAST reference hits,
fetched live from NCBI Entrez. The live Entrez call itself is never
exercised here -- build_sample_phylogeny_dataset() takes an
injectable retrieve_function, the same dependency-injection pattern
v0.1.0's own tests used, so these tests never touch the network.
"""

# Imports
from pathlib import Path
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from mitopipeline.utils.blast.blast_hits import write_top_blast_hits
from mitopipeline.utils.blast.models import BlastHit
from mitopipeline.utils.phylogeny.sample_dataset import (
    build_assembled_sequence_name,
    build_reference_sequence_name,
    build_sample_phylogeny_dataset,
    normalize_blast_accession,
    order_and_standardize_reference_sequences,
    read_assembled_genome,
    standardize_sequence_name,
    validate_sample_dataset,
)

def _hit(**overrides) -> BlastHit:
    defaults = dict(
        qseqid = "sample1", sseqid = "NC_000001.1", pident = 98.5, length = 16000,
        mismatch = 10, gapopen = 0, qstart = 1, qend = 16000, sstart = 1, send = 16000,
        evalue = 0.0, bitscore = 30000.0, qcovs = 99.0, staxids = "7998",
        sscinames = "Danio rerio", stitle = "Danio rerio mitochondrion, complete genome", rank = 1,
    )
    defaults.update(overrides)
    return BlastHit(**defaults)

def test_normalize_blast_accession_plain() -> None:
    assert normalize_blast_accession("NC_002333.2") == "NC_002333.2"

def test_normalize_blast_accession_ref_wrapper() -> None:
    assert normalize_blast_accession("ref|NC_002333.2|") == "NC_002333.2"

def test_normalize_blast_accession_empty_raises() -> None:
    with pytest.raises(ValueError):
        normalize_blast_accession("   ")

def test_standardize_sequence_name_replaces_special_characters() -> None:
    assert standardize_sequence_name("Danio rerio (mitochondrion)!") == "Danio_rerio_mitochondrion"

def test_build_assembled_sequence_name() -> None:
    assert build_assembled_sequence_name("sample 1") == "sample_1|assembled"

def test_build_reference_sequence_name() -> None:
    assert build_reference_sequence_name("NC_000001.1", "Danio rerio") == "Danio_rerio__NC_000001.1"

def test_read_assembled_genome(tmp_path: Path) -> None:
    fasta = tmp_path / "sample1.fasta"
    fasta.write_text(">contig_1\nATGCGT\n")
    record = read_assembled_genome(fasta, "sample1")
    assert record.id == "sample1|assembled"
    assert str(record.seq) == "ATGCGT"

def test_read_assembled_genome_missing_file_raises(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        read_assembled_genome(tmp_path / "does_not_exist.fasta", "sample1")

def test_read_assembled_genome_multiple_records_raises(tmp_path: Path) -> None:
    fasta = tmp_path / "sample1.fasta"
    fasta.write_text(">a\nATGC\n>b\nGGGG\n")
    with pytest.raises(ValueError):
        read_assembled_genome(fasta, "sample1")

def test_order_and_standardize_reference_sequences_orders_by_rank() -> None:
    hits = [_hit(sseqid = "NC_000002.1", rank = 2, sscinames = "Species two"),
            _hit(sseqid = "NC_000001.1", rank = 1, sscinames = "Species one")]
    retrieved = [SeqRecord(Seq("CCCC"), id = "NC_000002.1"), SeqRecord(Seq("AAAA"), id = "NC_000001.1")]

    ordered = order_and_standardize_reference_sequences(hits, retrieved)

    assert len(ordered) == 2
    assert ordered[0].id == "Species_two__NC_000002.1"
    assert ordered[1].id == "Species_one__NC_000001.1"

def test_order_and_standardize_reference_sequences_missing_accession_raises() -> None:
    hits = [_hit(sseqid = "NC_999999.1")]
    with pytest.raises(ValueError):
        order_and_standardize_reference_sequences(hits, [])

def test_validate_sample_dataset_requires_assembled_genome_first() -> None:
    records = [SeqRecord(Seq("AAAA"), id = "not_assembled"), SeqRecord(Seq("CCCC"), id = "ref__NC_1")]
    with pytest.raises(ValueError):
        validate_sample_dataset(records, expected_reference_count = 1)

def test_validate_sample_dataset_rejects_duplicate_identifiers() -> None:
    records = [SeqRecord(Seq("AAAA"), id = "s|assembled"), SeqRecord(Seq("CCCC"), id = "s|assembled")]
    with pytest.raises(ValueError):
        validate_sample_dataset(records, expected_reference_count = 1)

def test_build_sample_phylogeny_dataset(tmp_path: Path) -> None:
    """Confirm a complete combined FASTA (genome + references) is generated, without a real Entrez call."""
    assembly_fasta = tmp_path / "sample1.fasta"
    assembly_fasta.write_text(">contig_1\nATGCGT\n")

    top_hits_tsv = tmp_path / "sample1.top_hits.tsv"
    write_top_blast_hits(
        [_hit(sseqid = "NC_000001.1", rank = 1, sscinames = "Species one"),
         _hit(sseqid = "NC_000002.1", rank = 2, sscinames = "Species two")],
        top_hits_tsv, sample_id = "sample1",
    )

    output_fasta = tmp_path / "sample1.phylogeny_dataset.fasta"

    def fake_retrieve(accessions, entrez_email, entrez_api_key = None, logger = None):
        assert accessions == ["NC_000001.1", "NC_000002.1"]
        assert entrez_email == "test@example.com"
        return [SeqRecord(Seq("AAAA"), id = "NC_000001.1"), SeqRecord(Seq("CCCC"), id = "NC_000002.1")]

    reference_count = build_sample_phylogeny_dataset(
        sample_id = "sample1", assembly_fasta = assembly_fasta, top_hits_tsv = top_hits_tsv,
        output_fasta = output_fasta, entrez_email = "test@example.com", retrieve_function = fake_retrieve,
    )

    assert reference_count == 2
    from Bio import SeqIO
    records = list(SeqIO.parse(output_fasta, "fasta"))
    assert len(records) == 3
    assert records[0].id == "sample1|assembled"

def test_build_sample_phylogeny_dataset_no_hits_raises(tmp_path: Path) -> None:
    assembly_fasta = tmp_path / "sample1.fasta"
    assembly_fasta.write_text(">contig_1\nATGCGT\n")
    top_hits_tsv = tmp_path / "sample1.top_hits.tsv"
    write_top_blast_hits([], top_hits_tsv, sample_id = "sample1")

    with pytest.raises(ValueError):
        build_sample_phylogeny_dataset(
            sample_id = "sample1", assembly_fasta = assembly_fasta, top_hits_tsv = top_hits_tsv,
            output_fasta = tmp_path / "out.fasta", entrez_email = "test@example.com",
            retrieve_function = lambda **kwargs: [],
        )
