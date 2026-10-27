"""
test_run_blast.py

Fixture tests for the BLAST exec stage (#429): runs through BaseTool
and records a result distinguishing zero-hits success from tool
failure.

The real blastn binary is not available in this environment, so
BlastRunner.run() is stubbed (same approach used throughout this
project for external tools) -- these tests verify the exec script's
own wiring and result recording, not blastn's own behavior.
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

run_blast = importlib.import_module("mitopipeline.exec.05b-blast.run_blast")

def _fake_result(success: bool = True, return_code: int = 0) -> CommandResult:
    now = datetime.now()
    return CommandResult(
        command = [], return_code = return_code, stdout = "", stderr = "", runtime_seconds = 1.2,
        success = success, tool_name = "blast", started_at = now, ended_at = now,
        peak_memory_mb = 0.0, user_cpu_seconds = 0.0, system_cpu_seconds = 0.0,
    )

_TSV_COLUMNS = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend",
                "sstart", "send", "evalue", "bitscore", "qcovs", "staxids", "sscinames", "stitle"]

def _write_config(tmp_path: Path, sample_id: str, output_file: Path) -> Path:
    query = tmp_path / f"{sample_id}.fasta"
    query.write_text(">q\nACGT\n")
    config_path = tmp_path / f"{sample_id}.blast_config.json"
    config_path.write_text(json.dumps({
        "sample_id": sample_id, "query_fasta": str(query), "database": "actinopterygii_mt_refs",
        "mode": "local", "output_format": "6 " + " ".join(_TSV_COLUMNS),
        "output_columns": _TSV_COLUMNS, "output_file": str(output_file),
    }))
    return config_path

def _argv(tmp_path: Path, sample_id: str, config_json: Path) -> list[str]:
    return [
        "run_blast.py",
        "--sample-id", sample_id,
        "--config-json", str(config_json),
        "--output-dir", str(tmp_path / "out"),
        "--working-dir", str(tmp_path),
        "--threads", "2",
        "--log-file", str(tmp_path / "logs" / "blast.log"),
        "--global-log-file", str(tmp_path / "logs" / "global.log"),
    ]

def _write_hit_row(output_file: Path, sseqid: str = "NC_000001.1") -> None:
    row = ["q1", sseqid, "98.5", "500", "2", "0", "1", "501", "1", "501", "1e-50", "200", "95", "9606", "Homo sapiens", "title"]
    output_file.parent.mkdir(parents = True, exist_ok = True)
    output_file.write_text("\t".join(row) + "\n")

def test_successful_run_with_hits_records_hit_count(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    sample_id = "SAMPLE_01"
    output_file = tmp_path / "blast" / f"{sample_id}.blast.tsv"
    config_json = _write_config(tmp_path, sample_id, output_file)

    def _stub_run(self: BlastRunner) -> CommandResult:
        _write_hit_row(self.output_file)
        return _fake_result()

    monkeypatch.setattr(BlastRunner, "run", _stub_run)
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, config_json))

    assert run_blast.main() == 0
    result = json.loads((tmp_path / "out" / f"{sample_id}.blast_result.json").read_text())
    assert result["success"] is True
    assert result["hit_count"] == 1
    assert result["output_file"] == str(output_file)

def test_successful_run_with_zero_hits_is_not_a_failure(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    sample_id = "SAMPLE_02"
    output_file = tmp_path / "blast" / f"{sample_id}.blast.tsv"
    config_json = _write_config(tmp_path, sample_id, output_file)

    def _stub_run(self: BlastRunner) -> CommandResult:
        self.output_file.parent.mkdir(parents = True, exist_ok = True)
        self.output_file.write_text("")   # zero hits -- a valid, non-failure outcome
        return _fake_result()

    monkeypatch.setattr(BlastRunner, "run", _stub_run)
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, config_json))

    assert run_blast.main() == 0
    result = json.loads((tmp_path / "out" / f"{sample_id}.blast_result.json").read_text())
    assert result["success"] is True
    assert result["hit_count"] == 0

def test_tool_failure_is_recorded_not_fabricated(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    sample_id = "SAMPLE_03"
    output_file = tmp_path / "blast" / f"{sample_id}.blast.tsv"
    config_json = _write_config(tmp_path, sample_id, output_file)

    monkeypatch.setattr(BlastRunner, "run", lambda self: _fake_result(success = False, return_code = 1))
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, config_json))

    assert run_blast.main() == 1
    result = json.loads((tmp_path / "out" / f"{sample_id}.blast_result.json").read_text())
    assert result["success"] is False
    assert result["hit_count"] is None
    assert result["output_file"] is None

def test_validation_failure_systemexit_is_caught_not_raised(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    sample_id = "SAMPLE_04"
    output_file = tmp_path / "blast" / f"{sample_id}.blast.tsv"
    config_json = _write_config(tmp_path, sample_id, output_file)

    def _raise(self):
        raise SystemExit(1)
    monkeypatch.setattr(BlastRunner, "run", _raise)
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, config_json))

    assert run_blast.main() == 1   # must not raise out of main()

def test_missing_config_json_is_recorded_not_a_crash(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, "SAMPLE_05", tmp_path / "does_not_exist.json"))
    assert run_blast.main() == 1
