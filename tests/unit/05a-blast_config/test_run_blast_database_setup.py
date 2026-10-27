"""
test_run_blast_database_setup.py

Fixture tests for the job-level BLAST reference database setup exec:
wires setup_reference_database() to a CLI and records the outcome -- a
download/build failure must never fail the whole pipeline run, since
the biological validation gate's own "blast_identity" check already
treats missing BLAST evidence as a normal REVIEW outcome, not a crash.
"""

# Imports
import importlib
import json
import sys
from pathlib import Path
import pytest

run_blast_database_setup = importlib.import_module("mitopipeline.exec.05a-blast_config.run_blast_database_setup")

def _argv(tmp_path: Path) -> list[str]:
    return [
        "run_blast_database_setup.py",
        "--reference-fasta", str(tmp_path / "reference.fasta"),
        "--database-prefix", str(tmp_path / "db" / "actinopterygii_mt_refs"),
        "--database-title", "Test reference database",
        "--taxon-scope", "Actinopterygii",
        "--entrez-email", "test@example.com",
        "--max-records", "10000",
        "--batch-size", "50",
        "--parse-seqids",
        "--output-dir", str(tmp_path / "out"),
        "--log-file", str(tmp_path / "logs" / "setup.log"),
        "--global-log-file", str(tmp_path / "logs" / "global.log"),
    ]

def test_main_records_success_when_setup_succeeds(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(run_blast_database_setup, "setup_reference_database", lambda **kwargs: 0)
    monkeypatch.setattr(sys, "argv", _argv(tmp_path))

    assert run_blast_database_setup.main() == 0
    result = json.loads((tmp_path / "out" / "blast_database.json").read_text())
    assert result["success"] is True
    assert result["database_prefix"] == str(tmp_path / "db" / "actinopterygii_mt_refs")

def test_main_records_failure_without_a_nonzero_exit_on_a_nonzero_return_code(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(run_blast_database_setup, "setup_reference_database", lambda **kwargs: 1)
    monkeypatch.setattr(sys, "argv", _argv(tmp_path))

    assert run_blast_database_setup.main() == 0   # non-fatal, see module docstring
    result = json.loads((tmp_path / "out" / "blast_database.json").read_text())
    assert result["success"] is False

def test_main_records_failure_without_a_nonzero_exit_on_an_exception(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """A real NCBI network failure (DNS, timeout, malformed XML) must
    be equally non-fatal to the whole pipeline run."""
    def _raise(**kwargs):
        raise OSError("name resolution failed")
    monkeypatch.setattr(run_blast_database_setup, "setup_reference_database", _raise)
    monkeypatch.setattr(sys, "argv", _argv(tmp_path))

    assert run_blast_database_setup.main() == 0
    result = json.loads((tmp_path / "out" / "blast_database.json").read_text())
    assert result["success"] is False

def test_main_passes_through_the_real_cli_arguments(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    seen_kwargs = {}

    def _stub(**kwargs):
        seen_kwargs.update(kwargs)
        return 0

    monkeypatch.setattr(run_blast_database_setup, "setup_reference_database", _stub)
    monkeypatch.setattr(sys, "argv", _argv(tmp_path) + ["--entrez-api-key", "real-key", "--overwrite"])

    assert run_blast_database_setup.main() == 0
    assert seen_kwargs["taxon_scope"] == "Actinopterygii"
    assert seen_kwargs["email"] == "test@example.com"
    assert seen_kwargs["api_key"] == "real-key"
    assert seen_kwargs["max_records"] == 10000
    assert seen_kwargs["batch_size"] == 50
    assert seen_kwargs["parse_seqids"] is True
    assert seen_kwargs["overwrite"] is True
