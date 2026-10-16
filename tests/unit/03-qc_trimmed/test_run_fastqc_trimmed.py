"""
test_run_fastqc_trimmed.py

Unit tests for the run_fastqc_trimmed exec-stage script.
"""

# Imports
import importlib
import sys
from datetime import datetime
from pathlib import Path
import pytest
from mitopipeline.api.fastqc import FastQCRunner
from mitopipeline.models.command_result.command_result import CommandResult

# The exec/ stage directories are numbered ("03-qc_trimmed") to match the
# Snakemake rule files, which is not a valid dotted-import identifier --
# importlib.import_module takes a plain string, so it resolves fine even
# though `from mitopipeline.exec.03-qc_trimmed import ...` could never
# be written as literal Python syntax.
run_fastqc_trimmed = importlib.import_module(
    "mitopipeline.exec.03-qc_trimmed.run_fastqc_trimmed"
)

def _fake_success_result() -> CommandResult:
    """
    Builds a trivially successful CommandResult, used to stand in for a
    real FastQC invocation that this unit test does not run.

    Returns:
        CommandResult: A successful, zero-cost CommandResult.
    """
    now = datetime.now()
    return CommandResult(
        command = ["fastqc"],
        return_code = 0,
        stdout = "",
        stderr = "",
        runtime_seconds = 0.0,
        success = True,
        tool_name = "fastqc",
        started_at = now,
        ended_at = now,
        peak_memory_mb = 0.0,
        user_cpu_seconds = 0.0,
        system_cpu_seconds = 0.0,
    )

def test_main_isolates_trimmed_output_from_raw(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """
    main() should build the FastQCRunner with exactly the trimmed-stage
    output directory passed on the command line, never a raw-stage path
    -- the two must never collide for the same sample.
    """
    r1 = tmp_path / "SAMPLE_01.R1.trimmed.fastq.gz"
    r2 = tmp_path / "SAMPLE_01.R2.trimmed.fastq.gz"
    r1.touch()
    r2.touch()
    raw_dir = tmp_path / "qc" / "raw" / "SAMPLE_01"
    trimmed_dir = tmp_path / "qc" / "trimmed" / "SAMPLE_01"
    raw_dir.mkdir(parents = True)

    captured: dict = {}
    real_init = FastQCRunner.__init__

    def _capture_init(self, **kwargs):
        captured.update(kwargs)
        real_init(self, **kwargs)

    monkeypatch.setattr(FastQCRunner, "__init__", _capture_init)
    monkeypatch.setattr(FastQCRunner, "run", lambda self: _fake_success_result())
    monkeypatch.setattr(sys, "argv", [
        "run_fastqc_trimmed.py",
        "--sample-id", "SAMPLE_01",
        "--r1", str(r1),
        "--r2", str(r2),
        "--output-dir", str(trimmed_dir),
        "--working-dir", str(tmp_path),
        "--log-file", str(tmp_path / "logs" / "qc_trimmed.log"),
        "--global-log-file", str(tmp_path / "logs" / "global.log"),
        "--threads", "2",
    ])

    exit_code = run_fastqc_trimmed.main()

    assert exit_code == 0
    assert captured["output_dir"] == trimmed_dir
    assert captured["output_dir"] != raw_dir
    assert captured["sample"].sample_id == "SAMPLE_01"
    assert captured["sample"].r1 == r1.resolve()
    assert captured["sample"].r2 == r2.resolve()

def test_main_returns_failure_code_on_unsuccessful_run(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """
    main() should propagate a non-zero return code when the FastQC run
    itself fails, rather than always returning 0.
    """
    r1 = tmp_path / "SAMPLE_01.R1.trimmed.fastq.gz"
    r2 = tmp_path / "SAMPLE_01.R2.trimmed.fastq.gz"
    r1.touch()
    r2.touch()

    def _failed_run(self) -> CommandResult:
        result = _fake_success_result()
        result.success = False
        result.return_code = 2
        return result

    monkeypatch.setattr(FastQCRunner, "run", _failed_run)
    monkeypatch.setattr(sys, "argv", [
        "run_fastqc_trimmed.py",
        "--sample-id", "SAMPLE_01",
        "--r1", str(r1),
        "--r2", str(r2),
        "--output-dir", str(tmp_path / "qc" / "trimmed" / "SAMPLE_01"),
        "--working-dir", str(tmp_path),
        "--log-file", str(tmp_path / "logs" / "qc_trimmed.log"),
        "--global-log-file", str(tmp_path / "logs" / "global.log"),
        "--threads", "2",
    ])

    assert run_fastqc_trimmed.main() == 2
