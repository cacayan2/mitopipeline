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

def test_multi_scaffold_candidate_is_resolved_to_its_longest_record(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    # Confirmed live (etheostoma_tippecanoe_run_may2026_115): GetOrganelle
    # can leave multiple disjoint scaffold records in one data.fasta when
    # it can't link them into a single contig. Every downstream single-
    # sequence consumer (phylogeny dataset construction, MITOS2 annotation)
    # rejects a multi-record FASTA outright, so this step must resolve it
    # down to one real record -- the longest -- rather than pass the raw
    # multi-record file through.
    sample_id = "SAMPLE_06"
    attempt0 = tmp_path / "attempt0"
    attempt0.mkdir(parents = True)
    (attempt0 / f"{sample_id}.assembly_status.json").write_text(
        json.dumps({"sample_id": sample_id, "status": "fragmented", "success": True}))
    (attempt0 / "data.fasta").write_text(
        ">scaffold_1\n" + ("A" * 374) + "\n"
        ">scaffold_2\n" + ("C" * 414) + "\n"
        ">scaffold_5\n" + ("T" * 854) + "\n"
    )
    (attempt0 / "graph.gfa").write_text("H\tVN:Z:1.0\n")

    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, [attempt0]))
    assert run_contig_selection.main() == 0

    result = json.loads((tmp_path / "selection" / f"{sample_id}.contig_selection.json").read_text())
    assert result["length_bp"] == 854

    resolved_fasta = Path(result["fasta_path"])
    assert resolved_fasta != attempt0 / "data.fasta"
    resolved_text = resolved_fasta.read_text()
    assert resolved_text.count(">") == 1
    assert "scaffold_5" in resolved_text
    assert "A" * 374 not in resolved_text

def test_single_record_candidate_round_trips_unchanged(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    sample_id = "SAMPLE_07"
    attempt0 = tmp_path / "attempt0"
    _write_attempt(attempt0, sample_id, "circular", length = 16_596)

    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, [attempt0]))
    assert run_contig_selection.main() == 0

    result = json.loads((tmp_path / "selection" / f"{sample_id}.contig_selection.json").read_text())
    resolved_text = Path(result["fasta_path"]).read_text()
    assert resolved_text.count(">") == 1
    assert resolved_text.count("A") == 16_596

def test_status_says_circular_but_fasta_missing_is_treated_as_missing(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    sample_id = "SAMPLE_05"
    attempt0 = tmp_path / "attempt0"
    attempt0.mkdir(parents = True)
    (attempt0 / f"{sample_id}.assembly_status.json").write_text(
        json.dumps({"sample_id": sample_id, "status": "circular", "success": True}))
    # data.fasta deliberately never written -- an inconsistent/corrupt attempt dir

    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, [attempt0]))
    assert run_contig_selection.main() == 1
