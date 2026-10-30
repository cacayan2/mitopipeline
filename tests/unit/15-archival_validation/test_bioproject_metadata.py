"""
test_bioproject_metadata.py

Unit tests for project-level BioProject metadata: real sample-scope
resolution from what a run actually contains, honest missing-field
tracking, and never a fabricated release date/accession.
"""

# Imports
import json
from pathlib import Path
from mitopipeline.utils.archival.bioproject_metadata import build_bioproject_metadata, write_bioproject_metadata

def test_single_organism_multiple_samples_is_multiisolate() -> None:
    record = build_bioproject_metadata({"s1": "Danio rerio", "s2": "Danio rerio", "s3": "Danio rerio"})
    assert record.sample_scope == "Multiisolate"
    assert record.organism == "Danio rerio"
    assert record.sample_ids == "s1;s2;s3"

def test_single_sample_single_organism_is_monoisolate() -> None:
    record = build_bioproject_metadata({"s1": "Danio rerio"})
    assert record.sample_scope == "Monoisolate"

def test_multiple_organisms_is_multispecies() -> None:
    record = build_bioproject_metadata({"s1": "Danio rerio", "s2": "Carassius gibelio"})
    assert record.sample_scope == "Multispecies"
    assert "Carassius gibelio" in record.organism
    assert "Danio rerio" in record.organism

def test_missing_publication_and_funding_tracked_not_fabricated() -> None:
    record = build_bioproject_metadata({"s1": "Danio rerio"})
    assert "publication" in record.missing_metadata
    assert "funding" in record.missing_metadata
    assert record.publication == ""
    assert record.funding == ""
    # Never a fabricated public release date.
    assert record.release_status == "HUP"

def test_real_publication_and_funding_clears_missing_metadata() -> None:
    record = build_bioproject_metadata({"s1": "Danio rerio"}, publication = "Smith et al. 2026", funding = "NSF DEB-1234567")
    assert "publication" not in record.missing_metadata
    assert "funding" not in record.missing_metadata

def test_no_samples_tracked_as_missing() -> None:
    record = build_bioproject_metadata({})
    assert "sample_ids" in record.missing_metadata
    assert "organism" in record.missing_metadata
    assert record.organism == "unknown organism (not yet resolved)"

def test_write_bioproject_metadata_writes_real_json(tmp_path: Path) -> None:
    record = build_bioproject_metadata({"s1": "Danio rerio"})
    output_path = tmp_path / "bioproject_metadata.json"

    write_bioproject_metadata(record, output_path)

    assert output_path.is_file()
    data = json.loads(output_path.read_text())
    assert data["sample_scope"] == "Monoisolate"
    assert data["organism"] == "Danio rerio"
