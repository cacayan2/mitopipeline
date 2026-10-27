"""
test_run_mitofish_setup.py

Fixture tests for the job-level MitoFish database setup exec: resolves
a release (or uses a pinned one), prepares the database, and records
the outcome -- a download/network failure must never fail the whole
pipeline run, since seed selection already falls back to unseeded
when the database is unavailable.
"""

# Imports
import importlib
import json
import sys
from pathlib import Path
import pytest
from mitopipeline.utils.seed.mitofish_database import MitoFishDatabase

run_mitofish_setup = importlib.import_module("mitopipeline.exec.04c-seed_sequences.run_mitofish_setup")

def _argv(tmp_path: Path, version: str | None = None) -> list[str]:
    argv = [
        "run_mitofish_setup.py",
        "--database-root", str(tmp_path / "mitofish"),
        "--output-dir", str(tmp_path / "out"),
        "--log-file", str(tmp_path / "logs" / "setup.log"),
        "--global-log-file", str(tmp_path / "logs" / "global.log"),
    ]
    if version:
        argv += ["--version", version]
    return argv

def test_main_resolves_the_latest_version_and_prepares_the_database(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(run_mitofish_setup, "resolve_latest_mitofish_version", lambda logger: "2026.09")

    def _stub_prepare(self):
        return self.database_dir

    monkeypatch.setattr(MitoFishDatabase, "prepare", _stub_prepare)
    monkeypatch.setattr(sys, "argv", _argv(tmp_path))

    assert run_mitofish_setup.main() == 0
    result = json.loads((tmp_path / "out" / "mitofish_database.json").read_text())
    assert result["success"] is True
    assert result["version"] == "2026.09"
    assert result["database_dir"] == str(tmp_path / "mitofish" / "2026.09")

def test_main_uses_a_pinned_version_without_resolving_the_latest(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    def _fail_if_called(logger):
        raise AssertionError("must not resolve the latest version when one is pinned")
    monkeypatch.setattr(run_mitofish_setup, "resolve_latest_mitofish_version", _fail_if_called)
    monkeypatch.setattr(MitoFishDatabase, "prepare", lambda self: self.database_dir)
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, version = "2025.01"))

    assert run_mitofish_setup.main() == 0
    result = json.loads((tmp_path / "out" / "mitofish_database.json").read_text())
    assert result["version"] == "2025.01"

def test_main_records_failure_without_a_nonzero_exit_when_preparation_fails(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """A real network failure or a corrupt archive must not fail the
    whole pipeline run -- seed selection treats a missing database as
    a normal fallback-to-unseeded outcome."""
    monkeypatch.setattr(run_mitofish_setup, "resolve_latest_mitofish_version", lambda logger: "2026.09")

    def _stub_prepare(self):
        raise SystemExit(1)

    monkeypatch.setattr(MitoFishDatabase, "prepare", _stub_prepare)
    monkeypatch.setattr(sys, "argv", _argv(tmp_path))

    assert run_mitofish_setup.main() == 0
    result = json.loads((tmp_path / "out" / "mitofish_database.json").read_text())
    assert result["success"] is False

def test_main_records_failure_when_version_resolution_itself_fails(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """A real network failure reaching MitoFish's own stat API (not
    just a failed download) must be equally non-fatal."""
    def _fail(logger):
        raise OSError("name resolution failed")
    monkeypatch.setattr(run_mitofish_setup, "resolve_latest_mitofish_version", _fail)
    monkeypatch.setattr(sys, "argv", _argv(tmp_path))

    assert run_mitofish_setup.main() == 0
    result = json.loads((tmp_path / "out" / "mitofish_database.json").read_text())
    assert result["success"] is False
