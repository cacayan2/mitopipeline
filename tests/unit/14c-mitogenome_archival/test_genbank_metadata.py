"""
test_genbank_metadata.py

Unit tests for GenBank submission-folder preparation (#370/#373/
#449): a real selected assembly + MITOS2 GFF are collected into an
accession-independent folder; a missing selection or annotation is
tracked, not silently skipped or fabricated.
"""

# Imports
import json
from pathlib import Path
import pytest
from mitopipeline.utils.archival.genbank_metadata import build_genbank_submission, write_genbank_metadata

DEFAULTS = {
    "molecule_type": "DNA", "genome_location": "mitochondrion", "topology": "circular",
    "assembly_method": "GetOrganelle", "annotation_method": "MITOS2",
}

def _write_json(path: Path, data: dict) -> Path:
    path.parent.mkdir(parents = True, exist_ok = True)
    path.write_text(json.dumps(data))
    return path

def test_full_submission_folder_from_real_selection_and_annotation(tmp_path: Path) -> None:
    fasta = tmp_path / "data.fasta"
    fasta.write_text(">s\n" + "ACGT" * 4149 + "\n")
    gff = tmp_path / "result.gff"
    gff.write_text("##gff-version 3\n")

    selection_json = _write_json(tmp_path / "selection.json", {
        "selected": True, "status": "circular", "fasta_path": str(fasta), "length_bp": 16596,
    })
    annotation_json = _write_json(tmp_path / "annotation.json", {
        "success": True, "gff_path": str(gff),
    })
    output_dir = tmp_path / "genbank"

    record = build_genbank_submission("sample1", "Danio rerio", selection_json, annotation_json, output_dir, defaults = DEFAULTS)

    assert record["review_status"] == "READY"
    assert Path(record["assembly_fasta"]).is_file()
    assert Path(record["annotation_gff"]).is_file()
    assert record["sequence_length_bp"] == "16596"
    assert record["topology"] == "circular"
    assert "sample1_Danio_rerio" in record["submission_folder"]

def test_missing_selection_is_tracked_not_fabricated(tmp_path: Path) -> None:
    selection_json = _write_json(tmp_path / "selection.json", {"selected": False})
    annotation_json = tmp_path / "does_not_exist.json"
    output_dir = tmp_path / "genbank"

    record = build_genbank_submission("sample1", "Danio rerio", selection_json, annotation_json, output_dir, defaults = DEFAULTS)

    assert record["review_status"] == "MANUAL_REVIEW_REQUIRED"
    assert "assembly_fasta" in record["missing_files"]
    assert "annotation_gff" in record["missing_files"]
    assert record["assembly_fasta"] == ""

def test_missing_configured_defaults_are_tracked_in_metadata(tmp_path: Path) -> None:
    fasta = tmp_path / "data.fasta"
    fasta.write_text(">s\nACGT\n")
    selection_json = _write_json(tmp_path / "selection.json", {
        "selected": True, "status": "circular", "fasta_path": str(fasta), "length_bp": 4,
    })
    output_dir = tmp_path / "genbank"

    record = build_genbank_submission("sample1", "", selection_json, tmp_path / "missing.json", output_dir, defaults = {})

    assert "organism" in record["missing_metadata"]

def test_invalid_mode_raises(tmp_path: Path) -> None:
    selection_json = _write_json(tmp_path / "selection.json", {"selected": False})
    with pytest.raises(ValueError):
        build_genbank_submission("sample1", "Danio rerio", selection_json, tmp_path / "missing.json",
                                 tmp_path / "out", mode = "move")

def test_write_genbank_metadata_writes_a_real_tsv(tmp_path: Path) -> None:
    fasta = tmp_path / "data.fasta"
    fasta.write_text(">s\nACGT\n")
    selection_json = _write_json(tmp_path / "selection.json", {
        "selected": True, "status": "circular", "fasta_path": str(fasta), "length_bp": 4,
    })
    record = build_genbank_submission("sample1", "Danio rerio", selection_json, tmp_path / "missing.json",
                                      tmp_path / "genbank", defaults = DEFAULTS)
    output_path = tmp_path / "genbank_metadata.tsv"

    write_genbank_metadata([record], output_path)

    assert output_path.is_file()
    assert "sample_id" in output_path.read_text().splitlines()[0]
