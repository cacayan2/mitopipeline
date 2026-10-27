"""
test_run_biosample_metadata.py

Fixture tests for the BioSample metadata exec (#369/#372).
"""

# Imports
import importlib
import sys
from pathlib import Path
import pytest

run_biosample_metadata = importlib.import_module("mitopipeline.exec.14b-biosample_archival.run_biosample_metadata")

def test_builds_biosample_metadata_table(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text("sample_id\tgenus\tspecies\tcollection_date\tgeo_loc_name\tisolation_source\n"
                        "sample1\tDanio\trerio\t2026-06-01\tUSA: Illinois\tmuscle\n")
    output_dir = tmp_path / "output"

    monkeypatch.setattr(sys, "argv", [
        "run_biosample_metadata.py", "--runtime-manifest", str(manifest),
        "--output-dir", str(output_dir), "--log-file", str(tmp_path / "logs/s.log"),
        "--global-log-file", str(tmp_path / "logs/g.log"),
    ])

    assert run_biosample_metadata.main() == 0
    text = (output_dir / "biosample_metadata.tsv").read_text()
    assert "sample1" in text
    assert "READY" in text

def test_missing_manifest_returns_nonzero(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(sys, "argv", [
        "run_biosample_metadata.py", "--runtime-manifest", str(tmp_path / "missing.tsv"),
        "--output-dir", str(tmp_path / "output"), "--log-file", str(tmp_path / "logs/s.log"),
        "--global-log-file", str(tmp_path / "logs/g.log"),
    ])
    assert run_biosample_metadata.main() == 1
