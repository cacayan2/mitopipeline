"""
test_reference_database.py

Unit tests for building and validating the local BLAST reference
database.
"""

# Imports
import logging
from pathlib import Path
from unittest.mock import patch
import pytest
from mitopipeline.utils.blast.reference_database import (
    DEFAULT_TAXON_SCOPE,
    build_blast_database,
    build_ncbi_query,
    count_fasta_records,
    download_reference_fasta,
    setup_reference_database,
    validate_database,
)

def test_build_ncbi_query_defaults_to_actinopterygii() -> None:
    """
    The default query is scoped to Actinopterygii -- broad enough to
    include real confound species (e.g. marine fish for a freshwater
    target), not narrowed to a single genus/species.
    """
    query = build_ncbi_query()
    assert "Actinopterygii[Organism]" in query
    assert "mitochondrion[filter]" in query
    assert '"complete genome"[Title]' in query

def test_build_ncbi_query_accepts_custom_scope() -> None:
    """
    A caller can scope the query to a different taxonomic group.
    """
    query = build_ncbi_query("Mammalia")
    assert "Mammalia[Organism]" in query

def test_count_fasta_records_counts_headers(tmp_path: Path) -> None:
    """
    count_fasta_records counts lines starting with ">", not sequence lines.
    """
    path = tmp_path / "refs.fasta"
    path.write_text(">seq1\nACGT\n>seq2\nGGCC\n")
    assert count_fasta_records(path) == 2

def test_count_fasta_records_empty_file(tmp_path: Path) -> None:
    """
    An empty file has zero records.
    """
    path = tmp_path / "refs.fasta"
    path.write_text("")
    assert count_fasta_records(path) == 0

def test_download_reference_fasta_skips_existing_file_by_default(tmp_path: Path) -> None:
    """
    An existing non-empty reference FASTA is reused, not
    re-downloaded, unless overwrite is set.
    """
    reference_fasta = tmp_path / "refs.fasta"
    reference_fasta.write_text(">existing\nACGT\n")
    with patch("mitopipeline.utils.blast.reference_database.ncbi_esearch") as mock_esearch:
        result = download_reference_fasta(
            reference_fasta = reference_fasta, taxon_scope = "Actinopterygii", email = "test@example.com",
            api_key = None, max_records = 100, batch_size = 50, overwrite = False,
            logger = logging.getLogger("test"),
        )
        mock_esearch.assert_not_called()
    assert result == 0

def test_download_reference_fasta_requires_email(tmp_path: Path) -> None:
    """
    Downloading without a contact email is rejected outright -- NCBI
    Entrez requires one.
    """
    result = download_reference_fasta(
        reference_fasta = tmp_path / "refs.fasta", taxon_scope = "Actinopterygii", email = "",
        api_key = None, max_records = 100, batch_size = 50, overwrite = False,
        logger = logging.getLogger("test"),
    )
    assert result == 1

def test_download_reference_fasta_writes_fetched_records(tmp_path: Path) -> None:
    """
    A successful search-then-fetch writes every returned record to
    the reference FASTA.
    """
    reference_fasta = tmp_path / "refs.fasta"
    with patch("mitopipeline.utils.blast.reference_database.ncbi_esearch", return_value = ["1", "2"]), \
         patch("mitopipeline.utils.blast.reference_database.ncbi_efetch_fasta",
               return_value = ">seq1 Semotilus atromaculatus\nACGT\n>seq2 Esox lucius\nGGCC\n"):
        result = download_reference_fasta(
            reference_fasta = reference_fasta, taxon_scope = "Actinopterygii", email = "test@example.com",
            api_key = None, max_records = 100, batch_size = 50, overwrite = False,
            logger = logging.getLogger("test"),
        )
    assert result == 0
    assert count_fasta_records(reference_fasta) == 2

def test_download_reference_fasta_fails_on_empty_search(tmp_path: Path) -> None:
    """
    A query matching zero records is a real failure, not silently
    accepted as an empty reference set.
    """
    with patch("mitopipeline.utils.blast.reference_database.ncbi_esearch", return_value = []):
        result = download_reference_fasta(
            reference_fasta = tmp_path / "refs.fasta", taxon_scope = "Actinopterygii", email = "test@example.com",
            api_key = None, max_records = 100, batch_size = 50, overwrite = False,
            logger = logging.getLogger("test"),
        )
    assert result == 1

def test_download_reference_fasta_overwrite_forces_redownload(tmp_path: Path) -> None:
    """
    overwrite = True re-downloads even when a non-empty reference
    FASTA already exists.
    """
    reference_fasta = tmp_path / "refs.fasta"
    reference_fasta.write_text(">stale\nACGT\n")
    with patch("mitopipeline.utils.blast.reference_database.ncbi_esearch", return_value = ["1"]), \
         patch("mitopipeline.utils.blast.reference_database.ncbi_efetch_fasta",
               return_value = ">fresh\nTTTT\n"):
        result = download_reference_fasta(
            reference_fasta = reference_fasta, taxon_scope = "Actinopterygii", email = "test@example.com",
            api_key = None, max_records = 100, batch_size = 50, overwrite = True,
            logger = logging.getLogger("test"),
        )
    assert result == 0
    assert "fresh" in reference_fasta.read_text()
    assert "stale" not in reference_fasta.read_text()

@pytest.mark.slow
def test_build_blast_database_creates_a_queryable_database(tmp_path: Path) -> None:
    """
    makeblastdb, run for real on a tiny synthetic FASTA, produces a
    database that validates successfully via blastdbcmd. Requires
    BLAST+ on PATH -- excluded from the default test run.
    """
    reference_fasta = tmp_path / "refs.fasta"
    reference_fasta.write_text(">seq1 Semotilus atromaculatus\n" + "ACGTACGTACGT" * 20 + "\n")
    database_prefix = tmp_path / "db" / "test_refs"
    result = build_blast_database(
        reference_fasta = reference_fasta, database_prefix = database_prefix,
        database_title = "Test reference database", parse_seqids = True, logger = logging.getLogger("test"),
    )
    assert result == 0
    assert validate_database(database_prefix, logging.getLogger("test"))

@pytest.mark.slow
def test_setup_reference_database_end_to_end(tmp_path: Path) -> None:
    """
    The full orchestration -- mocked download, real makeblastdb build
    and validation, metadata written -- for one taxonomic scope.
    """
    database_prefix = tmp_path / "db" / "test_refs"
    metadata_file = tmp_path / "metadata.json"
    reference_fasta = tmp_path / "refs.fasta"
    with patch("mitopipeline.utils.blast.reference_database.ncbi_esearch", return_value = ["1"]), \
         patch("mitopipeline.utils.blast.reference_database.ncbi_efetch_fasta",
               return_value = ">seq1 Semotilus atromaculatus\n" + "ACGTACGTACGT" * 20 + "\n"):
        result = setup_reference_database(
            reference_fasta = reference_fasta, database_prefix = database_prefix,
            database_title = "Test reference database", metadata_file = metadata_file,
            email = "test@example.com", logger = logging.getLogger("test"),
        )
    assert result == 0
    assert validate_database(database_prefix, logging.getLogger("test"))
    assert metadata_file.is_file()
    assert DEFAULT_TAXON_SCOPE in metadata_file.read_text()
