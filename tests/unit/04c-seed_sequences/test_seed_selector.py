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
    assert result.fallback_reason == "no_fetchable_mitofish_match"
    assert result.accession == "GETORGANELLE_DEFAULT"
    assert result.fasta_path is None

def test_select_seed_resolves_multiple_candidates_deterministically(mitofish_database: Path, tmp_path: Path) -> None:
    """
    Regression (Emil's direct ask, 2026-10-28): more than one matching
    candidate for the same organism/gene -- a real, common case for a
    well-documented organism (e.g. Danio rerio, 16 complete MitoFish
    genomes) -- must still resolve to a real seed, not fall back to
    unseeded just because more than one reference exists. The longer,
    more complete record (ACCESSION_001, 33bp) is selected over the
    shorter one (ACCESSION_002, 30bp) deterministically.
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
    assert result.fallback_used is False
    assert result.accession == "ACCESSION_001"
    assert result.match_level == "best_of_multiple"
    assert result.fasta_path is not None
    assert result.fasta_path.is_file()

def test_select_seed_tries_the_next_best_candidate_when_the_best_is_unfetchable(mitofish_database: Path, tmp_path: Path) -> None:
    """
    When the longest/best candidate at a locus has no sequence in the
    local FASTA (a partially-synced database), the next-best candidate
    at the SAME locus is tried before moving on to a different locus
    -- here ACCESSION_001 (33bp, "best") is unfetchable, so the real
    seed written should be the shorter ACCESSION_002 (30bp) instead of
    an unseeded fallback.
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
    # Replacing the FASTA so only ACCESSION_002's sequence is present --
    # ACCESSION_001 (the "best" candidate by length) is now unfetchable.
    with gzip.open(mitofish_database / "mitofishdb.fa.gz", "wt") as handle:
        handle.write(">ACCESSION_002 Semotilus atromaculatus 12S\nACGTACGTACGTACGTACGTACGTACGTAC\n")

    result = select_seed(
        sample_id = "SAMPLE_09",
        scientific_name = "Semotilus atromaculatus",
        database_dir = mitofish_database,
        output_dir = tmp_path / "seeds",
        logger = logging.getLogger("test"),
    )
    assert result.fallback_used is False
    assert result.accession == "ACCESSION_002"

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
    A candidate accession present in the metadata tables but absent
    from the sequence FASTA (a partially-synced database) must not
    crash -- it tries the next locus (a different locus points to a
    different accession, which may still be present), and since this
    fixture only has 12S_rRNA data, falls all the way through the
    rest of the chain to the final fallback.
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
    assert result.fallback_reason == "no_fetchable_mitofish_match"

@pytest.fixture
def mitofish_database_coxi(tmp_path: Path) -> Path:
    """
    Like mitofish_database, but the organism only has a COXI
    annotation -- not 12S_rRNA -- so a correct chain walk must try
    COXI (the most densely documented locus) and succeed there
    without ever needing 12S_rRNA.

    Args:
        tmp_path (Path): The pytest-provided temporary directory.

    Returns:
        Path: The path to the synthetic database directory.
    """
    db_dir = tmp_path / "mitofish_db_coxi"
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
        "accession": ["ACCESSION_COXI"],
        "taxon_id": ["789"],
    }).to_parquet(db_dir / "seq_taxonid.parquet")
    pd.DataFrame({
        "accession": ["ACCESSION_COXI"],
        "gene": ["COXI"],
        "length": [len(_SEED_SEQUENCE)],
    }).to_parquet(db_dir / "seq_annotation.parquet")

    with gzip.open(db_dir / "mitofishdb.fa.gz", "wt") as handle:
        handle.write(f">ACCESSION_COXI Semotilus atromaculatus COXI\n{_SEED_SEQUENCE}\n")

    return db_dir

def test_select_seed_tries_most_documented_locus_first(mitofish_database_coxi: Path, tmp_path: Path) -> None:
    """
    With no loci override, the default chain starts at COXI (the most
    densely documented locus in the real local MitoFish database) --
    an organism with only a COXI annotation should resolve there
    directly, never needing to reach 12S_rRNA at all.
    """
    result = select_seed(
        sample_id = "SAMPLE_06",
        scientific_name = "Semotilus atromaculatus",
        database_dir = mitofish_database_coxi,
        output_dir = tmp_path / "seeds",
        logger = logging.getLogger("test"),
    )
    assert result.fallback_used is False
    assert result.gene == "COXI"
    assert result.accession == "ACCESSION_COXI"

def test_select_seed_falls_through_empty_loci_to_a_deeper_match(mitofish_database: Path, tmp_path: Path) -> None:
    """
    The default mitofish_database fixture only has 12S_rRNA data --
    COXI, Cytb, and 16S_rRNA (all tried first, being more densely
    documented) are empty for this organism, so a correct chain walk
    falls through all three before succeeding at 12S_rRNA.
    """
    result = select_seed(
        sample_id = "SAMPLE_07",
        scientific_name = "Semotilus atromaculatus",
        database_dir = mitofish_database,
        output_dir = tmp_path / "seeds",
        logger = logging.getLogger("test"),
    )
    assert result.fallback_used is False
    assert result.gene == "12S_rRNA"

def test_select_seed_accepts_custom_loci_override(mitofish_database_coxi: Path, tmp_path: Path) -> None:
    """
    A caller can still restrict the chain to specific loci -- passing
    only 12S_rRNA should fall back here, since this fixture has no
    12S_rRNA annotation at all (only COXI).
    """
    result = select_seed(
        sample_id = "SAMPLE_08",
        scientific_name = "Semotilus atromaculatus",
        database_dir = mitofish_database_coxi,
        output_dir = tmp_path / "seeds",
        logger = logging.getLogger("test"),
        loci = ("12S_rRNA",),
    )
    assert result.fallback_used is True
    assert result.fallback_reason == "no_fetchable_mitofish_match"
