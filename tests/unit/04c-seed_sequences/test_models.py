"""
test_models.py

Unit tests for the seed-sequence selection data contract.
"""

# Imports
from pathlib import Path
import pytest
from mitopipeline.utils.seed.models import DEFAULT_SEED_LOCUS, SeedSelectionResult

def test_default_seed_locus_is_configurable_not_hardcoded() -> None:
    """
    The default locus should be a single named constant that callers can
    reference, not a literal repeated at every call site.
    """
    assert DEFAULT_SEED_LOCUS == "12S_rRNA"

def test_seed_selection_result_stores_primary_source_fields(tmp_path: Path) -> None:
    """
    A result from a primary (non-fallback) source should record the
    accession, fasta path, source, and leave fallback_used False.
    """
    fasta_path = tmp_path / "SAMPLE_01.seed.fasta"
    result = SeedSelectionResult(
        sample_id = "SAMPLE_01",
        scientific_name = "Semotilus atromaculatus",
        gene = DEFAULT_SEED_LOCUS,
        accession = "ACCESSION_001",
        fasta_path = fasta_path,
        source = "mitofish",
        ncbi_taxon_id = 789,
        sequence_length = 954,
    )
    assert result.fallback_used is False
    assert result.source == "mitofish"
    assert result.fasta_path == fasta_path
    assert result.ncbi_taxon_id == 789

def test_seed_selection_result_allows_missing_optional_metadata(tmp_path: Path) -> None:
    """
    A fallback result may not have a known NCBI taxonomy ID or sequence
    length -- the contract must allow constructing a result without them.
    """
    result = SeedSelectionResult(
        sample_id = "SAMPLE_02",
        scientific_name = "Unidentified species",
        gene = DEFAULT_SEED_LOCUS,
        accession = "FALLBACK_DEFAULT",
        fasta_path = tmp_path / "SAMPLE_02.seed.fasta",
        source = "fallback_default_database",
        fallback_used = True,
    )
    assert result.fallback_used is True
    assert result.ncbi_taxon_id is None
    assert result.sequence_length is None

def test_seed_selection_result_is_immutable(tmp_path: Path) -> None:
    """
    The result contract is frozen -- once recorded, a seed selection
    outcome must not be mutated in place.
    """
    result = SeedSelectionResult(
        sample_id = "SAMPLE_01",
        scientific_name = "Semotilus atromaculatus",
        gene = DEFAULT_SEED_LOCUS,
        accession = "ACCESSION_001",
        fasta_path = tmp_path / "SAMPLE_01.seed.fasta",
        source = "mitofish",
    )
    with pytest.raises(AttributeError):
        result.fallback_used = True
