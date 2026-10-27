"""
test_run_mitogenome_metadata.py

Fixture tests for the GenBank submission-folder exec (#370/#373).
"""

# Imports
import importlib
import json
import sys
from pathlib import Path
import pytest

run_mitogenome_metadata = importlib.import_module("mitopipeline.exec.14c-mitogenome_archival.run_mitogenome_metadata")

def test_builds_real_submission_folder_and_record(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text("sample_id\tgenus\tspecies\nsample1\tDanio\trerio\n")

    fasta = tmp_path / "data.fasta"
    fasta.write_text(">s\nACGT\n")
    selection_json = tmp_path / "selection.json"
    selection_json.write_text(json.dumps({"selected": True, "status": "circular", "fasta_path": str(fasta), "length_bp": 4}))
    annotation_json = tmp_path / "missing_annotation.json"
    output_dir = tmp_path / "output"

    monkeypatch.setattr(sys, "argv", [
        "run_mitogenome_metadata.py", "--sample-id", "sample1", "--runtime-manifest", str(manifest),
        "--contig-selection-json", str(selection_json), "--annotation-result-json", str(annotation_json),
        "--output-dir", str(output_dir), "--log-file", str(tmp_path / "logs/s.log"),
        "--global-log-file", str(tmp_path / "logs/g.log"),
    ])

    assert run_mitogenome_metadata.main() == 0
    record = json.loads((output_dir / "sample1.genbank_record.json").read_text())
    assert record["organism"] == "Danio rerio"
    assert record["review_status"] == "MANUAL_REVIEW_REQUIRED"
    assert Path(record["assembly_fasta"]).is_file()
