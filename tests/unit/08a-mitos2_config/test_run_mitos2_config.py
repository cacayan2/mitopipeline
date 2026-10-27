"""
test_run_mitos2_config.py

Fixture tests for the MITOS2 configuration exec (#431): deterministic
taxonomy/genetic-code/FASTA resolution, and bad-config rejection.
"""

# Imports
import importlib
import json
import sys
from pathlib import Path
import pytest

run_mitos2_config = importlib.import_module("mitopipeline.exec.08a-mitos2_config.run_mitos2_config")
resolve_config = run_mitos2_config.resolve_config

def _make_refdir(tmp_path: Path, refseqver: str = "refseq89m") -> Path:
    refdir = tmp_path / "resources" / "mitos"
    (refdir / refseqver).mkdir(parents = True)
    return refdir

# ---- resolve_config() -- the deterministic core ----

def test_resolve_config_is_deterministic_across_repeat_calls(tmp_path: Path) -> None:
    fasta = tmp_path / "selected.fasta"
    fasta.write_text(">s\nACGT\n")
    refdir = _make_refdir(tmp_path)
    first = resolve_config("SAMPLE_01", fasta, 2, refdir, "refseq89m", True, {}, tmp_path / "annotation")
    second = resolve_config("SAMPLE_01", fasta, 2, refdir, "refseq89m", True, {}, tmp_path / "annotation")
    assert first == second

def test_resolve_config_records_tool_options_without_mutating_caller_dict(tmp_path: Path) -> None:
    fasta = tmp_path / "selected.fasta"
    fasta.write_text(">s\nACGT\n")
    refdir = _make_refdir(tmp_path)
    caller_options = {"noplots": True}
    config = resolve_config("SAMPLE_01", fasta, 2, refdir, "refseq89m", True, caller_options, tmp_path / "annotation")
    config["tool_options"]["best"] = True
    assert "best" not in caller_options

def test_resolve_config_rejects_missing_fasta(tmp_path: Path) -> None:
    refdir = _make_refdir(tmp_path)
    with pytest.raises(FileNotFoundError):
        resolve_config("SAMPLE_01", tmp_path / "does_not_exist.fasta", 2, refdir, "refseq89m", True, {}, tmp_path / "annotation")

def test_resolve_config_rejects_invalid_genetic_code(tmp_path: Path) -> None:
    fasta = tmp_path / "selected.fasta"
    fasta.write_text(">s\nACGT\n")
    refdir = _make_refdir(tmp_path)
    with pytest.raises(ValueError, match = "Invalid genetic code"):
        resolve_config("SAMPLE_01", fasta, 0, refdir, "refseq89m", True, {}, tmp_path / "annotation")

def test_resolve_config_rejects_missing_taxonomy_reference_dir(tmp_path: Path) -> None:
    fasta = tmp_path / "selected.fasta"
    fasta.write_text(">s\nACGT\n")
    refdir = _make_refdir(tmp_path)
    with pytest.raises(ValueError, match = "No such reference data directory"):
        resolve_config("SAMPLE_01", fasta, 2, refdir, "does_not_exist", True, {}, tmp_path / "annotation")

def test_resolve_config_is_json_serializable(tmp_path: Path) -> None:
    fasta = tmp_path / "selected.fasta"
    fasta.write_text(">s\nACGT\n")
    refdir = _make_refdir(tmp_path)
    config = resolve_config("SAMPLE_01", fasta, 2, refdir, "refseq89m", True, {"noplots": True}, tmp_path / "annotation")
    json.dumps(config)  # must not raise

# ---- main() -- the exec-script wiring ----

def test_main_writes_a_resolved_config_json(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    fasta = tmp_path / "selected.fasta"
    fasta.write_text(">s\nACGT\n")
    refdir = _make_refdir(tmp_path)
    output_dir = tmp_path / "config_out"
    monkeypatch.setattr(sys, "argv", [
        "run_mitos2_config.py",
        "--sample-id", "SAMPLE_01",
        "--input-fasta", str(fasta),
        "--genetic-code", "2",
        "--refdir", str(refdir),
        "--refseqver", "refseq89m",
        "--annotation-output-dir", str(tmp_path / "annotation"),
        "--output-dir", str(output_dir),
        "--log-file", str(tmp_path / "logs" / "config.log"),
        "--global-log-file", str(tmp_path / "logs" / "global.log"),
    ])
    assert run_mitos2_config.main() == 0
    config = json.loads((output_dir / "SAMPLE_01.mitos2_config.json").read_text())
    assert config["genetic_code"] == 2
    assert config["circular"] is True

def test_main_linear_flag_is_recorded(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    fasta = tmp_path / "selected.fasta"
    fasta.write_text(">s\nACGT\n")
    refdir = _make_refdir(tmp_path)
    output_dir = tmp_path / "config_out"
    monkeypatch.setattr(sys, "argv", [
        "run_mitos2_config.py",
        "--sample-id", "SAMPLE_02",
        "--input-fasta", str(fasta),
        "--genetic-code", "2",
        "--refdir", str(refdir),
        "--refseqver", "refseq89m",
        "--linear",
        "--annotation-output-dir", str(tmp_path / "annotation"),
        "--output-dir", str(output_dir),
        "--log-file", str(tmp_path / "logs" / "config.log"),
        "--global-log-file", str(tmp_path / "logs" / "global.log"),
    ])
    assert run_mitos2_config.main() == 0
    config = json.loads((output_dir / "SAMPLE_02.mitos2_config.json").read_text())
    assert config["circular"] is False

def test_main_returns_nonzero_on_invalid_genetic_code(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    fasta = tmp_path / "selected.fasta"
    fasta.write_text(">s\nACGT\n")
    refdir = _make_refdir(tmp_path)
    monkeypatch.setattr(sys, "argv", [
        "run_mitos2_config.py",
        "--sample-id", "SAMPLE_03",
        "--input-fasta", str(fasta),
        "--genetic-code", "0",
        "--refdir", str(refdir),
        "--refseqver", "refseq89m",
        "--annotation-output-dir", str(tmp_path / "annotation"),
        "--output-dir", str(tmp_path / "config_out"),
        "--log-file", str(tmp_path / "logs" / "config.log"),
        "--global-log-file", str(tmp_path / "logs" / "global.log"),
    ])
    assert run_mitos2_config.main() == 1
