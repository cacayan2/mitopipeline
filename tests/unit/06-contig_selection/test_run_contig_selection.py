"""
test_run_contig_selection.py

Fixture tests for the contig-selection exec (#402): reads each
assembly attempt's recorded status and standardized output, selects
the winner, and records it.
"""

# Imports
import importlib
import json
import sys
from pathlib import Path
import pytest

run_contig_selection = importlib.import_module("mitopipeline.exec.06-contig_selection.run_contig_selection")

def _argv(tmp_path: Path, sample_id: str, attempt_dirs: list[Path]) -> list[str]:
    argv = ["run_contig_selection.py", "--sample-id", sample_id]
    for d in attempt_dirs:
        argv += ["--attempt-dir", str(d)]
    argv += [
        "--output-dir", str(tmp_path / "selection"),
        "--log-file", str(tmp_path / "logs" / "contig_selection.log"),
        "--global-log-file", str(tmp_path / "logs" / "global.log"),
    ]
    return argv

def _write_attempt(attempt_dir: Path, sample_id: str, status: str, length: int = 100) -> None:
    attempt_dir.mkdir(parents = True, exist_ok = True)
    (attempt_dir / f"{sample_id}.assembly_status.json").write_text(
        json.dumps({"sample_id": sample_id, "status": status, "success": status != "missing"}))
    if status in ("circular", "fragmented"):
        (attempt_dir / "data.fasta").write_text(f">seq\n{'A' * length}\n")
        (attempt_dir / "graph.gfa").write_text("H\tVN:Z:1.0\n")

def test_single_circular_attempt_is_selected(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    sample_id = "SAMPLE_01"
    attempt0 = tmp_path / "attempt0"
    _write_attempt(attempt0, sample_id, "circular", length = 16_596)

    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, [attempt0]))
    assert run_contig_selection.main() == 0

    result = json.loads((tmp_path / "selection" / f"{sample_id}.contig_selection.json").read_text())
    assert result["selected"] is True
    assert result["status"] == "circular"
    assert result["length_bp"] == 16_596
    assert result["source_attempt"] == 0

def test_rescue_attempt_wins_when_better(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    sample_id = "SAMPLE_02"
    attempt0 = tmp_path / "attempt0"
    attempt1 = tmp_path / "attempt1"
    _write_attempt(attempt0, sample_id, "fragmented", length = 5_000)
    _write_attempt(attempt1, sample_id, "circular", length = 16_000)

    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, [attempt0, attempt1]))
    assert run_contig_selection.main() == 0

    result = json.loads((tmp_path / "selection" / f"{sample_id}.contig_selection.json").read_text())
    assert result["status"] == "circular"
    assert result["source_attempt"] == 1

def test_all_missing_is_recorded_as_unselected_and_fails(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    sample_id = "SAMPLE_03"
    attempt0 = tmp_path / "attempt0"
    _write_attempt(attempt0, sample_id, "missing")

    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, [attempt0]))
    assert run_contig_selection.main() == 1

    result = json.loads((tmp_path / "selection" / f"{sample_id}.contig_selection.json").read_text())
    assert result["selected"] is False
    assert result["attempts_considered"] == 1

def test_missing_status_file_is_treated_as_missing_not_a_crash(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    sample_id = "SAMPLE_04"
    attempt0 = tmp_path / "attempt0"
    attempt0.mkdir(parents = True)   # status file never written -- e.g. the assembly step crashed first

    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, [attempt0]))
    assert run_contig_selection.main() == 1

    result = json.loads((tmp_path / "selection" / f"{sample_id}.contig_selection.json").read_text())
    assert result["selected"] is False

def test_status_says_circular_but_fasta_missing_is_treated_as_missing(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    sample_id = "SAMPLE_05"
    attempt0 = tmp_path / "attempt0"
    attempt0.mkdir(parents = True)
    (attempt0 / f"{sample_id}.assembly_status.json").write_text(
        json.dumps({"sample_id": sample_id, "status": "circular", "success": True}))
    # data.fasta deliberately never written -- an inconsistent/corrupt attempt dir

    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, [attempt0]))
    assert run_contig_selection.main() == 1
