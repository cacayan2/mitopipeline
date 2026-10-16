"""
test_seed_selector.py

Unit tests for the MitoFish-aware seed selection orchestrator.
"""

# Imports
import gzip
import logging
from pathlib import Path
import pandas as pd
import pytest
from mitopipeline.utils.seed.models import DEFAULT_SEED_LOCUS
from mitopipeline.utils.seed.seed_selector import select_seed

_SEED_SEQUENCE = "ACGTACGTACGTACGTACGTACGTACGTACGT"

@pytest.fixture
def mitofish_database(tmp_path: Path) -> Path:
    """
    Creates a minimal synthetic MitoFish database (metadata tables plus
    the bgzipped sequence FASTA) for one known organism/accession.

    Represents: Semotilus atromaculatus, lineage_id 123, species_id 456,
    taxon_id 789, ACCESSION_001 (and a second accession for the
    ambiguous-match test).

    Args:
        tmp_path (Path): The pytest-provided temporary directory.

    Returns:
        Path: The path to the synthetic database directory.
    """
    db_dir = tmp_path / "mitofish_db"
    db_dir.mkdir()

    pd.DataFrame({
        "lineage_id": ["123"],
        "lineage_name": ["Semotilus atromaculatus"],
    }).to_parquet(db_dir / "taxonid_name.parquet")
    pd.DataFrame({
        "lineage_id": ["123"],
        "species_id": ["456"],
    }).to_parquet(db_dir / "speciesid_lineageid.parquet")
    pd.DataFrame({
        "species_id": ["456"],
        "taxon_id": ["789"],
    }).to_parquet(db_dir / "taxonid_speciesid.parquet")
    pd.DataFrame({
        "accession": ["ACCESSION_001"],
        "taxon_id": ["789"],
    }).to_parquet(db_dir / "seq_taxonid.parquet")
    pd.DataFrame({
        "accession": ["ACCESSION_001"],
        "gene": ["12S_rRNA"],
        "length": [len(_SEED_SEQUENCE)],
    }).to_parquet(db_dir / "seq_annotation.parquet")

    with gzip.open(db_dir / "mitofishdb.fa.gz", "wt") as handle:
        handle.write(f">ACCESSION_001 Semotilus atromaculatus 12S\n{_SEED_SEQUENCE}\n")

    return db_dir

def test_select_seed_writes_fasta_with_full_provenance(mitofish_database: Path, tmp_path: Path) -> None:
    """
    A known organism with exactly one matching candidate should produce
    a written FASTA and a result with the real accession, checksum,
    match level, and no fallback.
    """
    output_dir = tmp_path / "seeds"
    result = select_seed(
        sample_id = "SAMPLE_01",
        scientific_name = "Semotilus atromaculatus",
        database_dir = mitofish_database,
        output_dir = output_dir,
        logger = logging.getLogger("test"),
    )
    assert result.fallback_used is False
    assert result.fallback_reason is None
    assert result.accession == "ACCESSION_001"
    assert result.source == "mitofish"
    assert result.match_level == "exact"
    assert result.sequence_length == len(_SEED_SEQUENCE)
    assert result.checksum is not None
    assert result.fasta_path is not None
    assert result.fasta_path.is_file()
    assert _SEED_SEQUENCE in result.fasta_path.read_text()

def test_select_seed_falls_back_on_unknown_organism(mitofish_database: Path, tmp_path: Path) -> None:
    """
    An organism with no MitoFish candidate (missing locus for that
    organism, or the organism itself being absent) should fall back
    rather than raise, and write no FASTA.
    """
    result = select_seed(
        sample_id = "SAMPLE_02",
        scientific_name = "Definitely notreal species",
        database_dir = mitofish_database,
        output_dir = tmp_path / "seeds",
        logger = logging.getLogger("test"),
    )
    assert result.fallback_used is True
    assert result.fallback_reason == "no_mitofish_candidate"
    assert result.accession == "GETORGANELLE_DEFAULT"
    assert result.fasta_path is None

def test_select_seed_falls_back_on_ambiguous_candidates(mitofish_database: Path, tmp_path: Path) -> None:
    """
    More than one matching candidate for the same organism/gene is an
    ambiguous match that must fall back rather than silently guess one.
    """
    annotation_path = mitofish_database / "seq_annotation.parquet"
    pd.concat([
        pd.read_parquet(annotation_path),
        pd.DataFrame({"accession": ["ACCESSION_002"], "gene": ["12S_rRNA"], "length": [30]}),
    ]).to_parquet(annotation_path)
    pd.concat([
        pd.read_parquet(mitofish_database / "seq_taxonid.parquet"),
        pd.DataFrame({"accession": ["ACCESSION_002"], "taxon_id": ["789"]}),
    ]).to_parquet(mitofish_database / "seq_taxonid.parquet")

    result = select_seed(
        sample_id = "SAMPLE_03",
        scientific_name = "Semotilus atromaculatus",
        database_dir = mitofish_database,
        output_dir = tmp_path / "seeds",
        logger = logging.getLogger("test"),
    )
    assert result.fallback_used is True
    assert result.fallback_reason == "ambiguous_candidates"

def test_select_seed_falls_back_on_lookup_failure(tmp_path: Path) -> None:
    """
    A missing/corrupt/unreachable local database (standing in for a
    download or network failure) must fall back, never crash the pipeline
    over one sample.
    """
    result = select_seed(
        sample_id = "SAMPLE_04",
        scientific_name = "Semotilus atromaculatus",
        database_dir = tmp_path / "nonexistent_db",
        output_dir = tmp_path / "seeds",
        logger = logging.getLogger("test"),
    )
    assert result.fallback_used is True
    assert result.fallback_reason == "mitofish_lookup_failed"
    assert result.fasta_path is None

def test_select_seed_falls_back_when_accession_missing_from_sequence_file(mitofish_database: Path, tmp_path: Path) -> None:
    """
    A candidate accession present in the metadata tables but absent from
    the sequence FASTA (a partially-synced database) must fall back
    rather than crash.
    """
    with gzip.open(mitofish_database / "mitofishdb.fa.gz", "wt") as handle:
        handle.write(">SOME_OTHER_ACCESSION\nACGT\n")

    result = select_seed(
        sample_id = "SAMPLE_05",
        scientific_name = "Semotilus atromaculatus",
        database_dir = mitofish_database,
        output_dir = tmp_path / "seeds",
        logger = logging.getLogger("test"),
    )
    assert result.fallback_used is True
    assert result.fallback_reason == "sequence_fetch_failed"

def test_select_seed_uses_default_locus_when_unspecified(mitofish_database: Path, tmp_path: Path) -> None:
    """
    Not passing a gene should use DEFAULT_SEED_LOCUS, matching the
    "configurable seed locus" contract from #394.
    """
    result = select_seed(
        sample_id = "SAMPLE_06",
        scientific_name = "Semotilus atromaculatus",
        database_dir = mitofish_database,
        output_dir = tmp_path / "seeds",
        logger = logging.getLogger("test"),
    )
    assert result.gene == DEFAULT_SEED_LOCUS
