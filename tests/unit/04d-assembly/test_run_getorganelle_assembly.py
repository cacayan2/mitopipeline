"""
test_run_getorganelle_assembly.py

Fixture tests for the GetOrganelle assembly exec (#421): runs through
BaseTool, classifies native output, and records runtime/failures.

The real get_organelle_from_reads.py binary is not available in this
environment, so GetOrganelleRunner.run() is stubbed (same approach as
test_preprocessing_chain.py) -- these tests verify the exec script's
own wiring and status recording, not GetOrganelle's own behavior.
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

run_getorganelle_assembly = importlib.import_module("mitopipeline.exec.04d-assembly.run_getorganelle_assembly")

def _fake_result(success: bool = True, return_code: int = 0) -> CommandResult:
    now = datetime.now()
    return CommandResult(
        command = [], return_code = return_code, stdout = "", stderr = "", runtime_seconds = 1.5,
        success = success, tool_name = "getorganelle", started_at = now, ended_at = now,
        peak_memory_mb = 0.0, user_cpu_seconds = 0.0, system_cpu_seconds = 0.0,
    )

def _write_config(tmp_path: Path, sample_id: str, assembly_output_dir: Path, seeded: bool = False) -> Path:
    config_path = tmp_path / f"{sample_id}.getorganelle_config.json"
    config_path.write_text(json.dumps({
        "sample_id": sample_id, "organelle_type": "animal_mt",
        "seeded": seeded, "seed_file": str(tmp_path / "seed.fasta") if seeded else None,
        "assembly_output_dir": str(assembly_output_dir),
    }))
    return config_path

def _argv(tmp_path: Path, sample_id: str, config_json: Path) -> list[str]:
    return [
        "run_getorganelle_assembly.py",
        "--sample-id", sample_id,
        "--r1", str(tmp_path / f"{sample_id}_R1.fastq.gz"),
        "--r2", str(tmp_path / f"{sample_id}_R2.fastq.gz"),
        "--config-json", str(config_json),
        "--working-dir", str(tmp_path),
        "--log-file", str(tmp_path / "logs" / "assembly.log"),
        "--global-log-file", str(tmp_path / "logs" / "global.log"),
        "--threads", "2",
    ]

def test_successful_circular_assembly_writes_status_json(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    sample_id = "SAMPLE_01"
    assembly_dir = tmp_path / "assembly" / sample_id
    config_json = _write_config(tmp_path, sample_id, assembly_dir)

    def _stub_run(self: GetOrganelleRunner) -> CommandResult:
        self.output_dir.mkdir(parents = True, exist_ok = True)
        (self.output_dir / "k75.complete.path_sequence.fasta").write_text(">seq\nACGT\n")
        (self.output_dir / "k75.complete.selected_graph.gfa").write_text("H\n")
        return _fake_result()

    monkeypatch.setattr(GetOrganelleRunner, "run", _stub_run)
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, config_json))

    assert run_getorganelle_assembly.main() == 0
    status = json.loads((assembly_dir / f"{sample_id}.assembly_status.json").read_text())
    assert status == {"sample_id": sample_id, "status": "circular", "success": True}

def test_successful_but_fragmented_assembly_is_still_exit_zero(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """GetOrganelle completing without error but without a complete
    candidate is a normal outcome, not a tool failure -- the rescue
    trigger (#405), not this script, decides what happens next."""
    sample_id = "SAMPLE_02"
    assembly_dir = tmp_path / "assembly" / sample_id
    config_json = _write_config(tmp_path, sample_id, assembly_dir)

    def _stub_run(self: GetOrganelleRunner) -> CommandResult:
        self.output_dir.mkdir(parents = True, exist_ok = True)
        (self.output_dir / "k75.path_sequence.fasta").write_text(">seq\nACGT\n")
        (self.output_dir / "k75.selected_graph.gfa").write_text("H\n")
        return _fake_result()

    monkeypatch.setattr(GetOrganelleRunner, "run", _stub_run)
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, config_json))

    assert run_getorganelle_assembly.main() == 0
    status = json.loads((assembly_dir / f"{sample_id}.assembly_status.json").read_text())
    assert status["status"] == "fragmented"
    assert status["success"] is True

def test_assembly_with_no_native_candidate_is_classified_missing_and_fails(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """GetOrganelleRunner.run() itself sys.exit(1)s when postprocess_outputs
    finds zero native candidates -- this must be caught, classified, and
    recorded as a real failure, not an uncaught crash."""
    sample_id = "SAMPLE_03"
    assembly_dir = tmp_path / "assembly" / sample_id
    config_json = _write_config(tmp_path, sample_id, assembly_dir)

    def _stub_run(self: GetOrganelleRunner) -> CommandResult:
        self.output_dir.mkdir(parents = True, exist_ok = True)
        raise SystemExit(1)  # no native candidate -> postprocess_outputs would exit

    monkeypatch.setattr(GetOrganelleRunner, "run", _stub_run)
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, config_json))

    assert run_getorganelle_assembly.main() == 1
    status = json.loads((assembly_dir / f"{sample_id}.assembly_status.json").read_text())
    assert status == {"sample_id": sample_id, "status": "missing", "success": False}

def test_nonzero_return_code_is_recorded_as_failure(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    sample_id = "SAMPLE_04"
    assembly_dir = tmp_path / "assembly" / sample_id
    config_json = _write_config(tmp_path, sample_id, assembly_dir)

    def _stub_run(self: GetOrganelleRunner) -> CommandResult:
        self.output_dir.mkdir(parents = True, exist_ok = True)
        return _fake_result(success = False, return_code = 2)

    monkeypatch.setattr(GetOrganelleRunner, "run", _stub_run)
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, config_json))

    assert run_getorganelle_assembly.main() == 2
    status = json.loads((assembly_dir / f"{sample_id}.assembly_status.json").read_text())
    assert status["success"] is False

def test_missing_config_json_fails_cleanly(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    sample_id = "SAMPLE_05"
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, tmp_path / "does_not_exist.json"))
    assert run_getorganelle_assembly.main() == 1

def test_seeded_config_passes_seed_file_through_to_the_runner(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    sample_id = "SAMPLE_06"
    assembly_dir = tmp_path / "assembly" / sample_id
    seed = tmp_path / "seed.fasta"
    seed.write_text(">seed\nACGT\n")
    config_json = _write_config(tmp_path, sample_id, assembly_dir, seeded = True)

    seen_tool_options = {}

    def _stub_run(self: GetOrganelleRunner) -> CommandResult:
        seen_tool_options.update(self.tool_options)
        self.output_dir.mkdir(parents = True, exist_ok = True)
        (self.output_dir / "k75.complete.path_sequence.fasta").write_text(">seq\nACGT\n")
        (self.output_dir / "k75.complete.selected_graph.gfa").write_text("H\n")
        return _fake_result()

    monkeypatch.setattr(GetOrganelleRunner, "run", _stub_run)
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, config_json))

    assert run_getorganelle_assembly.main() == 0
    assert seen_tool_options.get("seed_file") == str(seed)
