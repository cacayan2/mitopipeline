"""
test_run_getorganelle_config.py

Fixture tests for the GetOrganelle configuration exec (#420):
deterministic assembly parameters, seed/no-seed inputs, output paths.
"""

# Imports
import importlib
import json
import sys
from pathlib import Path
import pytest

run_getorganelle_config = importlib.import_module("mitopipeline.exec.04a-getorganelle_config.run_getorganelle_config")
resolve_config = run_getorganelle_config.resolve_config

# ---- resolve_config() -- the deterministic core ----

def test_resolve_config_is_deterministic_across_repeat_calls(tmp_path: Path) -> None:
    seed = tmp_path / "seed.fasta"
    seed.write_text(">seed\nACGT\n")
    first = resolve_config("SAMPLE_01", "animal_mt", seed, tmp_path / "assembly")
    second = resolve_config("SAMPLE_01", "animal_mt", seed, tmp_path / "assembly")
    assert first == second

def test_resolve_config_seeded_input_is_recorded(tmp_path: Path) -> None:
    seed = tmp_path / "seed.fasta"
    seed.write_text(">seed\nACGT\n")
    config = resolve_config("SAMPLE_01", "animal_mt", seed, tmp_path / "assembly")
    assert config["seeded"] is True
    assert config["seed_file"] == str(seed)

def test_resolve_config_no_seed_input_is_recorded(tmp_path: Path) -> None:
    config = resolve_config("SAMPLE_01", "animal_mt", None, tmp_path / "assembly")
    assert config["seeded"] is False
    assert config["seed_file"] is None

def test_resolve_config_rejects_unsupported_organelle_type(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match = "Unsupported organelle type"):
        resolve_config("SAMPLE_01", "not_a_real_type", None, tmp_path / "assembly")

def test_resolve_config_records_the_assembly_output_path(tmp_path: Path) -> None:
    assembly_dir = tmp_path / "assembly" / "SAMPLE_01"
    config = resolve_config("SAMPLE_01", "animal_mt", None, assembly_dir)
    assert config["assembly_output_dir"] == str(assembly_dir)

def test_resolve_config_is_json_serializable(tmp_path: Path) -> None:
    config = resolve_config("SAMPLE_01", "animal_mt", tmp_path / "seed.fasta", tmp_path / "assembly")
    json.dumps(config)  # must not raise

# ---- main() -- the exec-script wiring ----

def test_main_writes_a_resolved_config_json(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    seed = tmp_path / "seed.fasta"
    seed.write_text(">seed\nACGT\n")
    output_dir = tmp_path / "config_out"
    monkeypatch.setattr(sys, "argv", [
        "run_getorganelle_config.py",
        "--sample-id", "SAMPLE_01",
        "--organelle-type", "animal_mt",
        "--seed-file", str(seed),
        "--assembly-output-dir", str(tmp_path / "assembly"),
        "--output-dir", str(output_dir),
        "--log-file", str(tmp_path / "logs" / "config.log"),
        "--global-log-file", str(tmp_path / "logs" / "global.log"),
    ])
    assert run_getorganelle_config.main() == 0
    config_path = output_dir / "SAMPLE_01.getorganelle_config.json"
    assert config_path.is_file()
    config = json.loads(config_path.read_text())
    assert config["seeded"] is True

def test_main_without_seed_file_flag_resolves_an_unseeded_config(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    output_dir = tmp_path / "config_out"
    monkeypatch.setattr(sys, "argv", [
        "run_getorganelle_config.py",
        "--sample-id", "SAMPLE_02",
        "--organelle-type", "animal_mt",
        "--assembly-output-dir", str(tmp_path / "assembly"),
        "--output-dir", str(output_dir),
        "--log-file", str(tmp_path / "logs" / "config.log"),
        "--global-log-file", str(tmp_path / "logs" / "global.log"),
    ])
    assert run_getorganelle_config.main() == 0
    config = json.loads((output_dir / "SAMPLE_02.getorganelle_config.json").read_text())
    assert config["seeded"] is False

def test_main_returns_nonzero_on_unsupported_organelle_type(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(sys, "argv", [
        "run_getorganelle_config.py",
        "--sample-id", "SAMPLE_03",
        "--organelle-type", "not_a_real_type",
        "--assembly-output-dir", str(tmp_path / "assembly"),
        "--output-dir", str(tmp_path / "config_out"),
        "--log-file", str(tmp_path / "logs" / "config.log"),
        "--global-log-file", str(tmp_path / "logs" / "global.log"),
    ])
    assert run_getorganelle_config.main() == 1
