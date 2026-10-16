"""
test_mitofish.py

Unit tests for MitoFish mitochondrial reference lookup.
"""

# Imports
import logging
from pathlib import Path
import pandas as pd
import pytest
from mitopipeline.utils.seed.mitofish import MitoFishReferenceLookup

@pytest.fixture
def mitofish_database(tmp_path: Path) -> Path:
    """
    Creates a minimal synthetic MitoFish database for unit testing.

    The synthetic database represents:

        Semotilus atromaculatus
            - lineage_id 123
            - species_id 456
            - taxon_id 789
            - ACCESSION_001
            - ACCESSION_002

    Both accessions contain a 12S_rRNA annotation.

    Args:
        tmp_path (Path): The pytest-provided temporary directory.

    Returns:
        Path: The path to the synthetic database directory.
    """
    # Creating the MitoFish taxonomy-name table.
    pd.DataFrame({
        "lineage_id": ["123"],
        "lineage_name": ["Semotilus atromaculatus"],
    }).to_parquet(tmp_path / "taxonid_name.parquet")

    # Creating the MitoFish species-lineage table.
    pd.DataFrame({
        "lineage_id": ["123"],
        "species_id": ["456"],
    }).to_parquet(tmp_path / "speciesid_lineageid.parquet")

    # Creating the MitoFish taxon-species table.
    pd.DataFrame({
        "species_id": ["456"],
        "taxon_id": ["789"],
    }).to_parquet(tmp_path / "taxonid_speciesid.parquet")

    # Creating the MitoFish sequence-taxonomy table.
    pd.DataFrame({
        "accession": ["ACCESSION_001", "ACCESSION_002"],
        "taxon_id": ["789", "789"],
    }).to_parquet(tmp_path / "seq_taxonid.parquet")

    # Creating the MitoFish sequence-annotation table.
    pd.DataFrame({
        "accession": ["ACCESSION_001", "ACCESSION_002"],
        "gene": ["12S_rRNA", "12S_rRNA"],
        "length": [954, 951],
    }).to_parquet(tmp_path / "seq_annotation.parquet")

    return tmp_path

def _lookup(database_dir: Path) -> MitoFishReferenceLookup:
    """
    Builds a MitoFishReferenceLookup against the given database directory.

    Args:
        database_dir (Path): The path to the MitoFish database.

    Returns:
        MitoFishReferenceLookup: The constructed lookup.
    """
    return MitoFishReferenceLookup(database_dir = database_dir, logger = logging.getLogger("test"))

def test_lookup_unknown_species_returns_empty_list(mitofish_database: Path) -> None:
    """
    Unknown organisms should produce no reference candidates.
    """
    lookup = _lookup(mitofish_database)
    candidates = lookup.lookup(scientific_name = "Definitely notreal species", gene = "12S_rRNA")
    assert candidates == []

def test_lookup_known_species_returns_matching_candidates(mitofish_database: Path) -> None:
    """
    A known organism with matching gene annotations should return one
    candidate per matching accession.
    """
    lookup = _lookup(mitofish_database)
    candidates = lookup.lookup(scientific_name = "Semotilus atromaculatus", gene = "12S_rRNA")
    accessions = sorted(candidate.accession for candidate in candidates)
    assert accessions == ["ACCESSION_001", "ACCESSION_002"]
    assert all(candidate.ncbi_taxon_id == 789 for candidate in candidates)
    assert all(candidate.gene == "12S_rRNA" for candidate in candidates)

def test_lookup_known_species_wrong_gene_returns_empty_list(mitofish_database: Path) -> None:
    """
    A known organism without the requested gene should return no candidates.
    """
    lookup = _lookup(mitofish_database)
    candidates = lookup.lookup(scientific_name = "Semotilus atromaculatus", gene = "16S_rRNA")
    assert candidates == []

def test_lookup_missing_database_files_exits(tmp_path: Path) -> None:
    """
    A database directory missing the required taxonomy tables should exit
    rather than raise an unhandled pandas error.
    """
    lookup = _lookup(tmp_path)
    with pytest.raises(SystemExit):
        lookup.lookup(scientific_name = "Semotilus atromaculatus", gene = "12S_rRNA")

def test_lookup_empty_scientific_name_exits(mitofish_database: Path) -> None:
    """
    An empty scientific name should exit rather than query the database.
    """
    lookup = _lookup(mitofish_database)
    with pytest.raises(SystemExit):
        lookup.lookup(scientific_name = "   ", gene = "12S_rRNA")
