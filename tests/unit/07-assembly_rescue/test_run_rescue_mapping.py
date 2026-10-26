"""
test_run_rescue_mapping.py

Fixture tests for the rescue-mapping exec (#426): connects the rescue
trigger decision to the real minimap2/samtools chain, and records the
outcome for #406 to read next.
"""

# Imports
import importlib
import json
import sys
from pathlib import Path
import pytest
from mitopipeline.api.minimap2 import Minimap2Runner
from mitopipeline.api.samtools import SamtoolsRunner
from mitopipeline.models.command_result.command_result import CommandResult

run_rescue_mapping = importlib.import_module("mitopipeline.exec.07-assembly_rescue.run_rescue_mapping")

def _fake_result() -> CommandResult:
    from datetime import datetime
    now = datetime.now()
    return CommandResult(
        command = [], return_code = 0, stdout = "", stderr = "", runtime_seconds = 0.0,
        success = True, tool_name = "stub", started_at = now, ended_at = now,
        peak_memory_mb = 0.0, user_cpu_seconds = 0.0, system_cpu_seconds = 0.0,
    )

def _stub_minimap2_run(self: Minimap2Runner) -> CommandResult:
    self.output_dir.mkdir(parents = True, exist_ok = True)
    self.output_sam.write_text("fake sam")
    return _fake_result()

def _stub_samtools_run(self: SamtoolsRunner) -> CommandResult:
    self.output_dir.mkdir(parents = True, exist_ok = True)
    if self.operation == "sort":
        self.sorted_bam.write_bytes(b"fake bam")
    elif self.operation == "index":
        self.bam_index.write_bytes(b"fake index")
    else:
        self.rescued_r1.write_text("fake r1")
        self.rescued_r2.write_text("fake r2")
    return _fake_result()

def _argv(tmp_path: Path, sample_id: str, status: str, attempt: int = 0, max_attempts: int = 2) -> list[str]:
    r1 = tmp_path / f"{sample_id}_R1.fastq.gz"
    r2 = tmp_path / f"{sample_id}_R2.fastq.gz"
    r1.touch(); r2.touch()
    reference = tmp_path / "reference.fasta"
    reference.write_text(">ref\nACGT\n")
    return [
        "run_rescue_mapping.py",
        "--sample-id", sample_id,
        "--r1", str(r1), "--r2", str(r2),
        "--reference-fasta", str(reference),
        "--status", status,
        "--attempt", str(attempt),
        "--max-attempts", str(max_attempts),
        "--output-dir", str(tmp_path / "rescue"),
        "--working-dir", str(tmp_path),
        "--log-file", str(tmp_path / "logs" / "rescue.log"),
        "--global-log-file", str(tmp_path / "logs" / "global.log"),
    ]

def test_circular_status_is_never_rescued(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    sample_id = "SAMPLE_01"
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, "circular"))
    assert run_rescue_mapping.main() == 0
    result = json.loads((tmp_path / "rescue" / f"{sample_id}.rescue_mapping.json").read_text())
    assert result["attempted"] is False
    assert result["rescued_r1"] is None

def test_fragmented_status_triggers_the_real_chain(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(Minimap2Runner, "run", _stub_minimap2_run)
    monkeypatch.setattr(SamtoolsRunner, "run", _stub_samtools_run)
    sample_id = "SAMPLE_02"
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, "fragmented"))

    assert run_rescue_mapping.main() == 0

    result = json.loads((tmp_path / "rescue" / f"{sample_id}.rescue_mapping.json").read_text())
    assert result["attempted"] is True
    assert result["rescued_r1"].endswith(f"{sample_id}.rescue.R1.fastq.gz")
    assert result["rescued_r2"].endswith(f"{sample_id}.rescue.R2.fastq.gz")
    assert Path(result["rescued_r1"]).is_file()

def test_attempts_exhausted_is_not_rescued(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    sample_id = "SAMPLE_03"
    # attempt 2 of max_attempts 2 -- every profile already tried.
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, "fragmented", attempt = 2, max_attempts = 2))
    assert run_rescue_mapping.main() == 0
    result = json.loads((tmp_path / "rescue" / f"{sample_id}.rescue_mapping.json").read_text())
    assert result["attempted"] is False
    assert "exhausted" in result["reason"]

def test_chain_failure_is_recorded_not_an_uncaught_crash(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    def _failing_run(self):
        raise SystemExit(1)
    monkeypatch.setattr(Minimap2Runner, "run", _failing_run)
    sample_id = "SAMPLE_04"
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, "partial"))

    assert run_rescue_mapping.main() == 1

    result = json.loads((tmp_path / "rescue" / f"{sample_id}.rescue_mapping.json").read_text())
    assert result["attempted"] is True
    assert result["rescued_r1"] is None

def test_status_not_in_trigger_list_is_not_rescued(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    sample_id = "SAMPLE_05"
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, "missing"))
    assert run_rescue_mapping.main() == 0
    result = json.loads((tmp_path / "rescue" / f"{sample_id}.rescue_mapping.json").read_text())
    assert result["attempted"] is False
