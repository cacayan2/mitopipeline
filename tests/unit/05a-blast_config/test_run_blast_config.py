"""
test_run_blast_config.py

Fixture tests for the BLAST configuration exec (#428): deterministic
query/database/output-format resolution, and bad-config rejection.
"""

# Imports
import importlib
import json
import sys
from pathlib import Path
import pytest

run_blast_config = importlib.import_module("mitopipeline.exec.05a-blast_config.run_blast_config")
resolve_config = run_blast_config.resolve_config

# ---- resolve_config() -- the deterministic core ----

def test_resolve_config_is_deterministic_across_repeat_calls(tmp_path: Path) -> None:
    query = tmp_path / "query.fasta"
    query.write_text(">q\nACGT\n")
    first = resolve_config("SAMPLE_01", query, "actinopterygii_mt_refs", "local", tmp_path / "blast")
    second = resolve_config("SAMPLE_01", query, "actinopterygii_mt_refs", "local", tmp_path / "blast")
    assert first == second

def test_resolve_config_records_the_output_path(tmp_path: Path) -> None:
    query = tmp_path / "query.fasta"
    query.write_text(">q\nACGT\n")
    blast_dir = tmp_path / "blast" / "SAMPLE_01"
    config = resolve_config("SAMPLE_01", query, "actinopterygii_mt_refs", "local", blast_dir)
    assert config["output_file"] == str(blast_dir / "SAMPLE_01.blast.tsv")

def test_resolve_config_requests_the_full_parser_column_list(tmp_path: Path) -> None:
    query = tmp_path / "query.fasta"
    query.write_text(">q\nACGT\n")
    config = resolve_config("SAMPLE_01", query, "actinopterygii_mt_refs", "local", tmp_path / "blast")
    assert "qcovs" in config["output_format"]
    assert "staxids" in config["output_format"]
    assert len(config["output_columns"]) == 16

def test_resolve_config_rejects_unsupported_mode(tmp_path: Path) -> None:
    query = tmp_path / "query.fasta"
    query.write_text(">q\nACGT\n")
    with pytest.raises(ValueError, match = "Unsupported BLAST mode"):
        resolve_config("SAMPLE_01", query, "actinopterygii_mt_refs", "not_a_real_mode", tmp_path / "blast")

def test_resolve_config_rejects_missing_query_fasta(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        resolve_config("SAMPLE_01", tmp_path / "does_not_exist.fasta", "actinopterygii_mt_refs", "local", tmp_path / "blast")

def test_resolve_config_rejects_empty_database(tmp_path: Path) -> None:
    query = tmp_path / "query.fasta"
    query.write_text(">q\nACGT\n")
    with pytest.raises(ValueError, match = "No BLAST database"):
        resolve_config("SAMPLE_01", query, "", "local", tmp_path / "blast")

def test_resolve_config_is_json_serializable(tmp_path: Path) -> None:
    query = tmp_path / "query.fasta"
    query.write_text(">q\nACGT\n")
    config = resolve_config("SAMPLE_01", query, "actinopterygii_mt_refs", "local", tmp_path / "blast")
    json.dumps(config)  # must not raise

# ---- main() -- the exec-script wiring ----

def test_main_writes_a_resolved_config_json(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    query = tmp_path / "query.fasta"
    query.write_text(">q\nACGT\n")
    output_dir = tmp_path / "config_out"
    monkeypatch.setattr(sys, "argv", [
        "run_blast_config.py",
        "--sample-id", "SAMPLE_01",
        "--query-fasta", str(query),
        "--database", "actinopterygii_mt_refs",
        "--mode", "local",
        "--blast-output-dir", str(tmp_path / "blast"),
        "--output-dir", str(output_dir),
        "--log-file", str(tmp_path / "logs" / "config.log"),
        "--global-log-file", str(tmp_path / "logs" / "global.log"),
    ])
    assert run_blast_config.main() == 0
    config = json.loads((output_dir / "SAMPLE_01.blast_config.json").read_text())
    assert config["database"] == "actinopterygii_mt_refs"

def test_main_returns_nonzero_on_missing_query_fasta(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(sys, "argv", [
        "run_blast_config.py",
        "--sample-id", "SAMPLE_02",
        "--query-fasta", str(tmp_path / "does_not_exist.fasta"),
        "--database", "actinopterygii_mt_refs",
        "--blast-output-dir", str(tmp_path / "blast"),
        "--output-dir", str(tmp_path / "config_out"),
        "--log-file", str(tmp_path / "logs" / "config.log"),
        "--global-log-file", str(tmp_path / "logs" / "global.log"),
    ])
    assert run_blast_config.main() == 1
