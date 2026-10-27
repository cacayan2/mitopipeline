"""
test_run_circular_coverage.py

Fixture tests for the circular-coverage correction exec (#414): the
real BLAST binary isn't available in this environment, so BlastRunner
itself is stubbed (same approach used throughout this project for
external tools) -- these tests verify this script's own wiring
(skip logic for non-circular samples, result recording), not BLAST's
own behavior.
"""

# Imports
import importlib
import json
import sys
from datetime import datetime
from pathlib import Path
import pytest
from mitopipeline.api.blast import BlastRunner
from mitopipeline.models.command_result.command_result import CommandResult

run_circular_coverage = importlib.import_module("mitopipeline.exec.05c-blast_hits.run_circular_coverage")

def _fake_result(success: bool = True, return_code: int = 0) -> CommandResult:
    now = datetime.now()
    return CommandResult(
        command = [], return_code = return_code, stdout = "", stderr = "", runtime_seconds = 1.0,
        success = success, tool_name = "blast", started_at = now, ended_at = now,
        peak_memory_mb = 0.0, user_cpu_seconds = 0.0, system_cpu_seconds = 0.0,
    )

def _write_selection(tmp_path: Path, sample_id: str, status: str, fasta_text: str = ">s\nACGT\n") -> Path:
    fasta = tmp_path / f"{sample_id}.fasta"
    fasta.write_text(fasta_text)
    selection_path = tmp_path / f"{sample_id}.contig_selection.json"
    selection_path.write_text(json.dumps({
        "sample_id": sample_id, "selected": status is not None, "status": status,
        "fasta_path": str(fasta), "length_bp": len(fasta_text.splitlines()[1]) if status else None,
    }))
    return selection_path

def _argv(tmp_path: Path, sample_id: str, selection_path: Path, output_dir: Path) -> list[str]:
    return [
        "run_circular_coverage.py", "--sample-id", sample_id,
        "--contig-selection-json", str(selection_path), "--database", "fake_db",
        "--working-dir", str(tmp_path), "--output-dir", str(output_dir),
        "--log-file", str(tmp_path / "logs" / "cc.log"), "--global-log-file", str(tmp_path / "logs" / "global.log"),
    ]

def test_skips_cleanly_for_a_non_circular_sample(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    sample_id = "sample_fragmented"
    selection_path = _write_selection(tmp_path, sample_id, status = "fragmented")
    output_dir = tmp_path / "output"

    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, selection_path, output_dir))

    assert run_circular_coverage.main() == 0
    result = json.loads((output_dir / f"{sample_id}.circular_coverage_result.json").read_text())
    assert result["success"] is True
    assert result["corrected_query_coverage"] is None
    assert "not circular" in result["reason"]

def test_skips_cleanly_when_nothing_selected(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    sample_id = "sample_missing"
    selection_path = tmp_path / f"{sample_id}.contig_selection.json"
    selection_path.write_text(json.dumps({"sample_id": sample_id, "selected": False}))
    output_dir = tmp_path / "output"

    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, selection_path, output_dir))

    assert run_circular_coverage.main() == 0
    result = json.loads((output_dir / f"{sample_id}.circular_coverage_result.json").read_text())
    assert result["corrected_query_coverage"] is None

def test_computes_real_corrected_coverage_for_a_circular_sample(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    sample_id = "sample_circular"
    selection_path = _write_selection(tmp_path, sample_id, status = "circular", fasta_text = ">s\n" + ("ACGT" * 4149) + "\n")

    def _stub_run(self: BlastRunner) -> CommandResult:
        self.output_file.parent.mkdir(parents = True, exist_ok = True)
        self.output_file.write_text(
            "q\tref1\t100.0\t16596\t0\t0\t1\t16596\t1\t16596\t0.0\t30000.0\t100.0\t1\tExample species\tExample\n"
        )
        return _fake_result()

    monkeypatch.setattr(BlastRunner, "run", _stub_run)
    output_dir = tmp_path / "output"
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, selection_path, output_dir))

    assert run_circular_coverage.main() == 0
    result = json.loads((output_dir / f"{sample_id}.circular_coverage_result.json").read_text())
    assert result["success"] is True
    assert result["corrected_query_coverage"] == pytest.approx(100.0)
    assert result["genome_length"] == 16596

def test_tool_failure_is_recorded_not_fabricated(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    sample_id = "sample_circular"
    selection_path = _write_selection(tmp_path, sample_id, status = "circular", fasta_text = ">s\nACGT\n")

    monkeypatch.setattr(BlastRunner, "run", lambda self: _fake_result(success = False, return_code = 1))
    output_dir = tmp_path / "output"
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, selection_path, output_dir))

    assert run_circular_coverage.main() == 1
    result = json.loads((output_dir / f"{sample_id}.circular_coverage_result.json").read_text())
    assert result["success"] is False
    assert result["corrected_query_coverage"] is None
