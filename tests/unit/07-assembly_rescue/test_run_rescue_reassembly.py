"""
test_run_rescue_reassembly.py

Fixture tests for the rescue-reassembly exec (#406): re-runs
GetOrganelle on rescued reads with one getorganelle_profiles entry's
tuning, through the same GetOrganelleRunner the primary attempt uses.
"""

# Imports
import importlib
import json
import sys
from datetime import datetime
from pathlib import Path
import pytest
from mitopipeline.api.getorganelle import GetOrganelleRunner
from mitopipeline.models.command_result.command_result import CommandResult

run_rescue_reassembly = importlib.import_module("mitopipeline.exec.07-assembly_rescue.run_rescue_reassembly")

def _fake_result(success: bool = True, return_code: int = 0) -> CommandResult:
    now = datetime.now()
    return CommandResult(
        command = [], return_code = return_code, stdout = "", stderr = "", runtime_seconds = 1.5,
        success = success, tool_name = "getorganelle", started_at = now, ended_at = now,
        peak_memory_mb = 0.0, user_cpu_seconds = 0.0, system_cpu_seconds = 0.0,
    )

def _argv(tmp_path: Path, sample_id: str, output_dir: Path, **extra) -> list[str]:
    r1 = tmp_path / f"{sample_id}.rescue.R1.fastq"
    r2 = tmp_path / f"{sample_id}.rescue.R2.fastq"
    r1.touch(); r2.touch()
    argv = [
        "run_rescue_reassembly.py",
        "--sample-id", sample_id,
        "--r1", str(r1), "--r2", str(r2),
        "--profile-name", extra.get("profile_name", "reduced_coverage"),
        "--output-dir", str(output_dir),
        "--working-dir", str(tmp_path),
        "--threads", "2",
        "--log-file", str(tmp_path / "logs" / "rescue_reassembly.log"),
        "--global-log-file", str(tmp_path / "logs" / "global.log"),
    ]
    if "kmer" in extra:
        argv += ["--kmer", extra["kmer"]]
    if "target_coverage" in extra:
        argv += ["--target-coverage", str(extra["target_coverage"])]
    return argv

def test_successful_rescue_writes_status_with_profile_name(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    sample_id = "SAMPLE_01"
    output_dir = tmp_path / "rescue_attempt_1"

    def _stub_run(self: GetOrganelleRunner) -> CommandResult:
        self.output_dir.mkdir(parents = True, exist_ok = True)
        (self.output_dir / "k85.complete.path_sequence.fasta").write_text(">seq\nACGT\n")
        (self.output_dir / "k85.complete.selected_graph.gfa").write_text("H\n")
        from mitopipeline.api.getorganelle import discover_assembly_artifacts
        self.assembly_status = discover_assembly_artifacts(self.output_dir, self.logger)
        return _fake_result()

    monkeypatch.setattr(GetOrganelleRunner, "run", _stub_run)
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, output_dir,
                                            kmer = "21,45,65,85,105", target_coverage = 100))

    assert run_rescue_reassembly.main() == 0
    status = json.loads((output_dir / f"{sample_id}.assembly_status.json").read_text())
    assert status == {"sample_id": sample_id, "status": "circular", "success": True,
                       "profile_name": "reduced_coverage"}

def test_tool_options_carry_the_profile_tuning_through(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    sample_id = "SAMPLE_02"
    output_dir = tmp_path / "rescue_attempt_1"
    seen = {}

    def _stub_run(self: GetOrganelleRunner) -> CommandResult:
        seen.update(self.tool_options)
        self.output_dir.mkdir(parents = True, exist_ok = True)
        (self.output_dir / "k85.path_sequence.fasta").write_text(">seq\nACGT\n")
        (self.output_dir / "k85.selected_graph.gfa").write_text("H\n")
        from mitopipeline.api.getorganelle import discover_assembly_artifacts
        self.assembly_status = discover_assembly_artifacts(self.output_dir, self.logger)
        return _fake_result()

    monkeypatch.setattr(GetOrganelleRunner, "run", _stub_run)
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, output_dir,
                                            kmer = "21,33,55,85,115", target_coverage = 100))

    assert run_rescue_reassembly.main() == 0
    assert seen.get("kmer") == "21,33,55,85,115"
    assert seen.get("target_coverage") == 100.0

def test_no_native_candidate_is_classified_missing_and_fails(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    sample_id = "SAMPLE_03"
    output_dir = tmp_path / "rescue_attempt_1"

    def _stub_run(self: GetOrganelleRunner) -> CommandResult:
        self.output_dir.mkdir(parents = True, exist_ok = True)
        raise SystemExit(1)

    monkeypatch.setattr(GetOrganelleRunner, "run", _stub_run)
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, output_dir))

    assert run_rescue_reassembly.main() == 1
    status = json.loads((output_dir / f"{sample_id}.assembly_status.json").read_text())
    assert status["status"] == "missing"
    assert status["success"] is False

def test_two_rescue_attempts_never_collide(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Different attempts must write to different output dirs -- the
    caller's job (a Snakemake rule/orchestrator), but this confirms
    the script itself never assumes or hardcodes a shared path."""
    sample_id = "SAMPLE_04"
    attempt1_dir = tmp_path / "rescue_attempt_1"
    attempt2_dir = tmp_path / "rescue_attempt_2"

    def _stub_run(self: GetOrganelleRunner) -> CommandResult:
        self.output_dir.mkdir(parents = True, exist_ok = True)
        (self.output_dir / "k85.path_sequence.fasta").write_text(">seq\nACGT\n")
        (self.output_dir / "k85.selected_graph.gfa").write_text("H\n")
        from mitopipeline.api.getorganelle import discover_assembly_artifacts
        self.assembly_status = discover_assembly_artifacts(self.output_dir, self.logger)
        return _fake_result()

    monkeypatch.setattr(GetOrganelleRunner, "run", _stub_run)

    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, attempt1_dir, profile_name = "reduced_coverage"))
    assert run_rescue_reassembly.main() == 0
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, attempt2_dir, profile_name = "alternate_kmers"))
    assert run_rescue_reassembly.main() == 0

    status1 = json.loads((attempt1_dir / f"{sample_id}.assembly_status.json").read_text())
    status2 = json.loads((attempt2_dir / f"{sample_id}.assembly_status.json").read_text())
    assert status1["profile_name"] == "reduced_coverage"
    assert status2["profile_name"] == "alternate_kmers"
