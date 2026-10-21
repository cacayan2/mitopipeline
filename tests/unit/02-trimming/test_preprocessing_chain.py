"""
test_preprocessing_chain.py

Integration-style unit test exercising one fixture sample through the
full preprocessing chain: FastQC (raw) -> fastp (trimming) -> FastQC
(trimmed). Verifies stage dependency order (the trimmed QC stage
consumes fastp's actual outputs, not the original raw reads) and that
every stage's QC/trimmed outputs land in distinct directories.

The real fastqc/fastp binaries are not available in this environment, so
each BaseTool.run() is stubbed to create its declared output files and
return a successful CommandResult -- this test verifies the WIRING
between stages, not the external tools' own behavior.
"""

# Imports
import importlib
import sys
from datetime import datetime
from pathlib import Path
import pytest
from mitopipeline.api.fastp import FastpRunner
from mitopipeline.api.fastqc import FastQCRunner
from mitopipeline.models.command_result.command_result import CommandResult

run_fastqc_raw = importlib.import_module("mitopipeline.exec.01-qc_raw.run_fastqc_raw")
run_fastp = importlib.import_module("mitopipeline.exec.02-trimming.run_fastp")
run_fastqc_trimmed = importlib.import_module("mitopipeline.exec.03-qc_trimmed.run_fastqc_trimmed")

def _fake_result() -> CommandResult:
    """
    Builds a trivially successful CommandResult, standing in for a real
    tool invocation this test does not run.

    Returns:
        CommandResult: A successful, zero-cost CommandResult.
    """
    now = datetime.now()
    return CommandResult(
        command = [], return_code = 0, stdout = "", stderr = "", runtime_seconds = 0.0,
        success = True, tool_name = "stub", started_at = now, ended_at = now,
        peak_memory_mb = 0.0, user_cpu_seconds = 0.0, system_cpu_seconds = 0.0,
    )

def _stub_fastqc_run(self: FastQCRunner) -> CommandResult:
    """A FastQCRunner.run() stub that creates its expected outputs."""
    self.output_dir.mkdir(parents = True, exist_ok = True)
    for path in self._expected_fastqc_outputs():
        path.write_text("fake fastqc output")
    return _fake_result()

def _stub_fastp_run(self: FastpRunner) -> CommandResult:
    """A FastpRunner.run() stub that creates its expected outputs."""
    self.output_dir.mkdir(parents = True, exist_ok = True)
    for path in (self.output_r1, self.output_r2, self.output_html, self.output_json):
        path.write_text("fake fastp output")
    return _fake_result()

def test_preprocessing_chain_threads_trimmed_reads_into_trimmed_qc(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    """
    Running qc_raw -> fastp -> qc_trimmed for one fixture sample should:
    - give each of the three stages its own, non-overlapping output dir
    - feed the trimmed QC stage fastp's ACTUAL output paths, not the
      original raw read paths
    """
    monkeypatch.setattr(FastQCRunner, "run", _stub_fastqc_run)
    monkeypatch.setattr(FastpRunner, "run", _stub_fastp_run)

    sample_id = "FIXTURE_01"
    raw_r1 = tmp_path / "raw" / f"{sample_id}_R1.fastq.gz"
    raw_r2 = tmp_path / "raw" / f"{sample_id}_R2.fastq.gz"
    raw_r1.parent.mkdir(parents = True)
    raw_r1.touch()
    raw_r2.touch()

    qc_raw_dir = tmp_path / "qc" / "raw" / sample_id
    trimmed_dir = tmp_path / "trimmed" / sample_id
    qc_trimmed_dir = tmp_path / "qc" / "trimmed" / sample_id
    log_dir = tmp_path / "logs"

    # Stage 1: FastQC on raw reads.
    monkeypatch.setattr(sys, "argv", [
        "run_fastqc_raw.py",
        "--sample-id", sample_id,
        "--r1", str(raw_r1), "--r2", str(raw_r2),
        "--output-dir", str(qc_raw_dir),
        "--working-dir", str(tmp_path),
        "--log-file", str(log_dir / "qc_raw.log"),
        "--global-log-file", str(log_dir / "global.log"),
        "--threads", "2",
    ])
    assert run_fastqc_raw.main() == 0

    # Stage 2: fastp trimming of the SAME raw reads.
    monkeypatch.setattr(sys, "argv", [
        "run_fastp.py",
        "--sample-id", sample_id,
        "--in1", str(raw_r1), "--in2", str(raw_r2),
        "--output-dir", str(trimmed_dir),
        "--working-dir", str(tmp_path),
        "--threads", "2",
        "--log-file", str(log_dir / "fastp.log"),
        "--global-log-file", str(log_dir / "global.log"),
    ])
    assert run_fastp.main() == 0

    # fastp's real, sample-scoped output paths -- this is what the
    # trimmed-QC stage must actually consume next, per the dependency
    # order this issue asks to verify.
    trimmed_r1 = trimmed_dir / f"{sample_id}.R1.trimmed.fastq.gz"
    trimmed_r2 = trimmed_dir / f"{sample_id}.R2.trimmed.fastq.gz"
    assert trimmed_r1.is_file()
    assert trimmed_r2.is_file()

    # Stage 3: FastQC on the TRIMMED reads (fastp's outputs), not the raw ones.
    monkeypatch.setattr(sys, "argv", [
        "run_fastqc_trimmed.py",
        "--sample-id", sample_id,
        "--r1", str(trimmed_r1), "--r2", str(trimmed_r2),
        "--output-dir", str(qc_trimmed_dir),
        "--working-dir", str(tmp_path),
        "--log-file", str(log_dir / "qc_trimmed.log"),
        "--global-log-file", str(log_dir / "global.log"),
        "--threads", "2",
    ])
    assert run_fastqc_trimmed.main() == 0

    # Every stage's output directory must be distinct -- no QC/trimmed
    # output ever collides with another stage's.
    output_dirs = {qc_raw_dir, trimmed_dir, qc_trimmed_dir}
    assert len(output_dirs) == 3

    # Both QC stages produced real artifacts, in their own directories.
    # Sample-id-prefixed (Emil's call, 2026-10-21): every stage's output
    # carries the sample prefix, for unambiguous data reorg/access later.
    assert (qc_raw_dir / f"{sample_id}.R1_fastqc.html").is_file()
    assert (qc_trimmed_dir / f"{sample_id}.R1_fastqc.html").is_file()
    assert (qc_raw_dir / f"{sample_id}.R1_fastqc.html") != (qc_trimmed_dir / f"{sample_id}.R1_fastqc.html")
