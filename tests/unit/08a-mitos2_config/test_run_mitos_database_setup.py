"""
test_run_mitos_database_setup.py

Fixture tests for the job-level MITOS2 reference data setup exec:
wires MitosReferenceData.prepare() to a CLI and records the outcome.
"""

# Imports
import importlib
import json
import sys
from pathlib import Path
import pytest
from mitopipeline.utils.annotation.mitos_database import MitosReferenceData

run_mitos_database_setup = importlib.import_module("mitopipeline.exec.08a-mitos2_config.run_mitos_database_setup")

def _argv(tmp_path: Path) -> list[str]:
    return [
        "run_mitos_database_setup.py",
        "--refdir", str(tmp_path / "mitos"),
        "--refseqver", "refseq89m",
        "--zenodo-record", "4284483",
        "--output-dir", str(tmp_path / "out"),
        "--log-file", str(tmp_path / "logs" / "setup.log"),
        "--global-log-file", str(tmp_path / "logs" / "global.log"),
    ]

def test_main_records_success_when_preparation_succeeds(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    def _stub_prepare(self):
        return self.version_dir
    monkeypatch.setattr(MitosReferenceData, "prepare", _stub_prepare)
    monkeypatch.setattr(sys, "argv", _argv(tmp_path))

    assert run_mitos_database_setup.main() == 0
    result = json.loads((tmp_path / "out" / "mitos_database.json").read_text())
    assert result["success"] is True
    assert result["refdir"] == str(tmp_path / "mitos" / "refseq89m")

def test_main_records_failure_without_a_nonzero_exit_on_a_systemexit(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    def _stub_prepare(self):
        raise SystemExit(1)
    monkeypatch.setattr(MitosReferenceData, "prepare", _stub_prepare)
    monkeypatch.setattr(sys, "argv", _argv(tmp_path))

    assert run_mitos_database_setup.main() == 0
    result = json.loads((tmp_path / "out" / "mitos_database.json").read_text())
    assert result["success"] is False

def test_main_records_failure_without_a_nonzero_exit_on_an_exception(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """A real Zenodo network failure must not fail the whole run."""
    def _stub_prepare(self):
        raise OSError("name resolution failed")
    monkeypatch.setattr(MitosReferenceData, "prepare", _stub_prepare)
    monkeypatch.setattr(sys, "argv", _argv(tmp_path))

    assert run_mitos_database_setup.main() == 0
    result = json.loads((tmp_path / "out" / "mitos_database.json").read_text())
    assert result["success"] is False
