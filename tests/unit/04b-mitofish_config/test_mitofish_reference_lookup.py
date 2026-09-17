"""
test_mitofish_reference_lookup.py

Unit tests for the MitoFishReferenceLookup class.
"""

# Imports
from pathlib import Path
from mitopipeline.utils.seed.mitofish import MitoFishReferenceLookup

def test_lookup_unknown_species_returns_empty_list(tmp_path: Path) -> None:
    """
    Unknown organisms should produce no reference candidates.
    """
    lookup = MitoFishReferenceLookup(database_dir = tmp_path)
    candidates = lookup.lookup(scientific_name="Definitely notreal species", gene = "12S")
    assert candidates == []