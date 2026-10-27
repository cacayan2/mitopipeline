"""
test_run_sra_metadata.py

Fixture tests for the SRA metadata exec (#368/#371).
"""

# Imports
import importlib
import json
import sys
from pathlib import Path
import pytest

run_sra_metadata = importlib.import_module("mitopipeline.exec.14a-sra_archival.run_sra_metadata")

def test_builds_sra_metadata_table(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    r1, r2 = tmp_path / "s1_R1.fastq.gz", tmp_path / "s1_R2.fastq.gz"
    r1.write_bytes(b""); r2.write_bytes(b"")
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text(f"sample_id\tr1\tr2\tgenus\tspecies\nsample1\t{r1}\t{r2}\tDanio\trerio\n")
    defaults_json = tmp_path / "defaults.json"
    defaults_json.write_text(json.dumps({"library_strategy": "WGS", "library_source": "GENOMIC",
                                         "library_selection": "RANDOM", "platform": "ILLUMINA",
                                         "instrument_model": "NovaSeq", "filetype": "fastq"}))
    output_dir = tmp_path / "output"

    monkeypatch.setattr(sys, "argv", [
        "run_sra_metadata.py", "--runtime-manifest", str(manifest), "--defaults-json", str(defaults_json),
        "--output-dir", str(output_dir), "--log-file", str(tmp_path / "logs/s.log"),
        "--global-log-file", str(tmp_path / "logs/g.log"),
    ])

    assert run_sra_metadata.main() == 0
    text = (output_dir / "sra_metadata.tsv").read_text()
    assert "sample1" in text
    assert "READY" in text

def test_missing_manifest_returns_nonzero(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(sys, "argv", [
        "run_sra_metadata.py", "--runtime-manifest", str(tmp_path / "missing.tsv"),
        "--output-dir", str(tmp_path / "output"), "--log-file", str(tmp_path / "logs/s.log"),
        "--global-log-file", str(tmp_path / "logs/g.log"),
    ])
    assert run_sra_metadata.main() == 1
