"""
test_run_mafft_alignment.py

Fixture tests for the MAFFT alignment exec (#348): runs through
BaseTool and validates the aligned output against the real input.
"""

# Imports
import importlib
import json
import sys
from datetime import datetime
from pathlib import Path
import pytest
from mitopipeline.api.mafft import MafftRunner
from mitopipeline.models.command_result.command_result import CommandResult

run_mafft_alignment = importlib.import_module("mitopipeline.exec.11-phylogeny.run_mafft_alignment")

def _fake_result(success: bool = True, return_code: int = 0) -> CommandResult:
    now = datetime.now()
    return CommandResult(command = [], return_code = return_code, stdout = "", stderr = "", runtime_seconds = 1.0,
                         success = success, tool_name = "mafft", started_at = now, ended_at = now,
                         peak_memory_mb = 0.0, user_cpu_seconds = 0.0, system_cpu_seconds = 0.0)

def _argv(tmp_path: Path, input_fasta: Path, output_fasta: Path) -> list[str]:
    return [
        "run_mafft_alignment.py", "--run-id", "test_run",
        "--input-fasta", str(input_fasta), "--output-fasta", str(output_fasta),
        "--threads", "2", "--working-dir", str(tmp_path), "--output-dir", str(tmp_path / "out"),
        "--log-file", str(tmp_path / "logs" / "mafft.log"), "--global-log-file", str(tmp_path / "logs" / "global.log"),
    ]

def test_successful_alignment_preserving_all_ids(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    input_fasta = tmp_path / "combined.fasta"
    input_fasta.write_text(">a\nACGT\n>b\nACGT\n")
    output_fasta = tmp_path / "aligned.fasta"

    def _stub_run(self: MafftRunner) -> CommandResult:
        self.output_fasta.parent.mkdir(parents = True, exist_ok = True)
        self.output_fasta.write_text(">a\nACGT\n>b\nACGT\n")
        return _fake_result()

    monkeypatch.setattr(MafftRunner, "run", _stub_run)
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, input_fasta, output_fasta))

    assert run_mafft_alignment.main() == 0
    result = json.loads((tmp_path / "out" / "mafft_alignment_result.json").read_text())
    assert result["success"] is True
    assert result["sequence_count"] == 2

def test_aligned_output_missing_a_sequence_is_rejected(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    input_fasta = tmp_path / "combined.fasta"
    input_fasta.write_text(">a\nACGT\n>b\nACGT\n")
    output_fasta = tmp_path / "aligned.fasta"

    def _stub_run(self: MafftRunner) -> CommandResult:
        self.output_fasta.parent.mkdir(parents = True, exist_ok = True)
        self.output_fasta.write_text(">a\nACGT\n")   # lost "b"
        return _fake_result()

    monkeypatch.setattr(MafftRunner, "run", _stub_run)
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, input_fasta, output_fasta))

    assert run_mafft_alignment.main() == 1
    result = json.loads((tmp_path / "out" / "mafft_alignment_result.json").read_text())
    assert result["success"] is False
    assert "doesn't match" in result["reason"]

def test_inconsistent_alignment_lengths_is_rejected(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    input_fasta = tmp_path / "combined.fasta"
    input_fasta.write_text(">a\nACGT\n>b\nACGT\n")
    output_fasta = tmp_path / "aligned.fasta"

    def _stub_run(self: MafftRunner) -> CommandResult:
        self.output_fasta.parent.mkdir(parents = True, exist_ok = True)
        self.output_fasta.write_text(">a\nACGT\n>b\nAC\n")   # not actually aligned
        return _fake_result()

    monkeypatch.setattr(MafftRunner, "run", _stub_run)
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, input_fasta, output_fasta))

    assert run_mafft_alignment.main() == 1
    result = json.loads((tmp_path / "out" / "mafft_alignment_result.json").read_text())
    assert "inconsistent lengths" in result["reason"]

def test_tool_failure_is_recorded_not_fabricated(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    input_fasta = tmp_path / "combined.fasta"
    input_fasta.write_text(">a\nACGT\n>b\nACGT\n")
    monkeypatch.setattr(MafftRunner, "run", lambda self: _fake_result(success = False, return_code = 1))
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, input_fasta, tmp_path / "aligned.fasta"))

    assert run_mafft_alignment.main() == 1
    result = json.loads((tmp_path / "out" / "mafft_alignment_result.json").read_text())
    assert result["success"] is False

def test_missing_input_fasta_is_recorded_not_a_crash(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, tmp_path / "does_not_exist.fasta", tmp_path / "aligned.fasta"))
    assert run_mafft_alignment.main() == 1
