"""
test_biosample_metadata.py

Unit tests for BioSample metadata generation (#369/#372/#448):
required environmental/collection metadata is tracked when missing,
never fabricated, and manifest values win over configured defaults.
"""

# Imports
from pathlib import Path
import pytest
from mitopipeline.utils.archival.biosample_metadata import build_biosample_metadata, write_biosample_metadata

DEFAULTS = {
    "sample_type": "environmental sample", "description_template": "Biological sample for {sample_id}.",
}

def test_build_biosample_metadata_with_full_manifest_data(tmp_path: Path) -> None:
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text(
        "sample_id\tgenus\tspecies\tcollection_date\tgeo_loc_name\tisolation_source\n"
        "sample1\tDanio\trerio\t2026-06-01\tUSA: Illinois\tmuscle tissue\n"
    )
    records = build_biosample_metadata(manifest, defaults = DEFAULTS)

    assert len(records) == 1
    record = records[0]
    assert record["sample_name"] == "sample1"
    assert record["organism"] == "Danio rerio"
    assert record["collection_date"] == "2026-06-01"
    assert record["review_status"] == "READY"
    assert record["missing_fields"] == ""

def test_missing_required_environmental_metadata_is_tracked(tmp_path: Path) -> None:
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text("sample_id\tgenus\tspecies\nsample1\tDanio\trerio\n")
    records = build_biosample_metadata(manifest, defaults = {})

    record = records[0]
    assert record["review_status"] == "MANUAL_REVIEW_REQUIRED"
    assert "collection_date" in record["missing_fields"]
    assert "geo_loc_name" in record["missing_fields"]
    assert "isolation_source" in record["missing_fields"]

def test_manifest_values_win_over_defaults(tmp_path: Path) -> None:
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text("sample_id\tgeo_loc_name\nsample1\tUSA: Illinois\n")
    records = build_biosample_metadata(manifest, defaults = {"geo_loc_name": "not collected"})
    assert records[0]["geo_loc_name"] == "USA: Illinois"

def test_description_defaults_to_template_when_not_provided(tmp_path: Path) -> None:
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text("sample_id\n sample1\n".replace(" sample1", "sample1"))
    records = build_biosample_metadata(manifest, defaults = DEFAULTS)
    assert records[0]["description"] == "Biological sample for sample1."

def test_blank_sample_id_raises(tmp_path: Path) -> None:
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text("sample_id\tgenus\n\tDanio\n")
    with pytest.raises(ValueError):
        build_biosample_metadata(manifest)

def test_missing_manifest_raises(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        build_biosample_metadata(tmp_path / "missing.tsv")

def test_write_biosample_metadata_writes_a_real_tsv(tmp_path: Path) -> None:
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text("sample_id\tgenus\tspecies\nsample1\tDanio\trerio\n")
    output_path = tmp_path / "biosample_metadata.tsv"

    records = write_biosample_metadata(manifest, output_path, defaults = DEFAULTS)

    assert output_path.is_file()
    assert len(records) == 1
