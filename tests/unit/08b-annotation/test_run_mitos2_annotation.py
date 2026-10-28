"""
test_run_mitos2_annotation.py

Fixture tests for the MITOS2 annotation exec (#432): runs through
BaseTool, parses real gene features out of its own GFF output, and
records runtime/failures/artifact checks.

The real runmitos.py binary is not available in this environment, so
Mitos2Runner.run() is stubbed (same approach used throughout this
project for external tools) -- these tests verify the exec script's
own wiring and result recording, not MITOS2's own behavior.
"""

# Imports
import importlib
import json
import sys
from datetime import datetime
from pathlib import Path
import pytest
from mitopipeline.api.mitos2 import Mitos2Runner
from mitopipeline.models.command_result.command_result import CommandResult

run_mitos2_annotation = importlib.import_module("mitopipeline.exec.08b-annotation.run_mitos2_annotation")

def _fake_result(success: bool = True, return_code: int = 0) -> CommandResult:
    now = datetime.now()
    return CommandResult(
        command = [], return_code = return_code, stdout = "", stderr = "", runtime_seconds = 2.0,
        success = success, tool_name = "mitos2", started_at = now, ended_at = now,
        peak_memory_mb = 0.0, user_cpu_seconds = 0.0, system_cpu_seconds = 0.0,
    )

def _write_config(tmp_path: Path, sample_id: str, annotation_output_dir: Path) -> Path:
    fasta = tmp_path / f"{sample_id}.fasta"
    fasta.write_text(">s\nACGT\n")
    refdir = tmp_path / "resources" / "mitos"
    (refdir / "refseq89m").mkdir(parents = True, exist_ok = True)
    config_path = tmp_path / f"{sample_id}.mitos2_config.json"
    config_path.write_text(json.dumps({
        "sample_id": sample_id, "input_fasta": str(fasta), "genetic_code": 2,
        "refdir": str(refdir), "refseqver": "refseq89m", "circular": True,
        "tool_options": {}, "annotation_output_dir": str(annotation_output_dir),
    }))
    return config_path

def _argv(tmp_path: Path, sample_id: str, config_json: Path) -> list[str]:
    return [
        "run_mitos2_annotation.py",
        "--sample-id", sample_id,
        "--config-json", str(config_json),
        "--working-dir", str(tmp_path),
        "--log-file", str(tmp_path / "logs" / "mitos2.log"),
        "--global-log-file", str(tmp_path / "logs" / "global.log"),
    ]

def _write_real_gff(result_gff: Path) -> None:
    result_gff.parent.mkdir(parents = True, exist_ok = True)
    result_gff.write_text(
        "##gff-version 3\n"
        "test\tmitos\tgene\t1\t600\t.\t+\t.\tID=gene_cox2;Name=cox2\n"
        "test\tmitos\tgene\t650\t1200\t.\t-\t.\tID=gene_nad1;Name=nad1\n"
    )

def test_successful_annotation_records_feature_count_and_genes(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    sample_id = "SAMPLE_01"
    output_dir = tmp_path / "annotation" / sample_id
    config_json = _write_config(tmp_path, sample_id, output_dir)

    def _stub_run(self: Mitos2Runner) -> CommandResult:
        _write_real_gff(self.result_gff)
        self.result_fas.write_text(">cox2\nACGT\n")
        return _fake_result()

    monkeypatch.setattr(Mitos2Runner, "run", _stub_run)
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, config_json))

    assert run_mitos2_annotation.main() == 0
    result = json.loads((output_dir / f"{sample_id}.mitos2_annotation_result.json").read_text())
    assert result["success"] is True
    assert result["feature_count"] == 2
    assert result["gene_names"] == ["cox2", "nad1"]
    assert result["trna_count"] == 0
    assert result["rrna_count"] == 0

def test_successful_annotation_records_trna_and_rrna_counts_separately(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Regression: feature_count only ever counted GFF type "gene"
    (protein-coding), silently mislabeled as "Total gene features" in
    the report -- tRNA/rRNA rows need their own real counts."""
    sample_id = "SAMPLE_03"
    output_dir = tmp_path / "annotation" / sample_id
    config_json = _write_config(tmp_path, sample_id, output_dir)

    def _stub_run(self: Mitos2Runner) -> CommandResult:
        self.result_gff.parent.mkdir(parents = True, exist_ok = True)
        self.result_gff.write_text(
            "##gff-version 3\n"
            "test\tmitos\tgene\t1\t600\t.\t+\t.\tID=gene_cox2;Name=cox2\n"
            "test\tmitfi\ttRNA\t650\t700\t.\t+\t.\tID=trnF;Name=trnF\n"
            "test\tmitfi\ttRNA\t710\t760\t.\t+\t.\tID=trnV;Name=trnV\n"
            "test\tmitfi\trRNA\t800\t1700\t.\t+\t.\tID=rrnS;Name=rrnS\n"
        )
        self.result_fas.write_text(">cox2\nACGT\n")
        return _fake_result()

    monkeypatch.setattr(Mitos2Runner, "run", _stub_run)
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, config_json))

    assert run_mitos2_annotation.main() == 0
    result = json.loads((output_dir / f"{sample_id}.mitos2_annotation_result.json").read_text())
    assert result["feature_count"] == 1
    assert result["trna_count"] == 2
    assert result["rrna_count"] == 1

def test_zero_features_is_recorded_not_a_crash(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    sample_id = "SAMPLE_02"
    output_dir = tmp_path / "annotation" / sample_id
    config_json = _write_config(tmp_path, sample_id, output_dir)

    def _stub_run(self: Mitos2Runner) -> CommandResult:
        self.result_gff.parent.mkdir(parents = True, exist_ok = True)
        self.result_gff.write_text("##gff-version 3\ntest\tmitos\tregion\t1\t1200\t.\t+\t.\tID=region1\n")
        self.result_fas.write_text(">x\nACGT\n")
        return _fake_result()

    monkeypatch.setattr(Mitos2Runner, "run", _stub_run)
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, config_json))

    assert run_mitos2_annotation.main() == 0
    result = json.loads((output_dir / f"{sample_id}.mitos2_annotation_result.json").read_text())
    assert result["feature_count"] == 0
    assert result["gene_names"] == []

def test_tool_failure_is_recorded_not_fabricated(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    sample_id = "SAMPLE_03"
    output_dir = tmp_path / "annotation" / sample_id
    config_json = _write_config(tmp_path, sample_id, output_dir)

    monkeypatch.setattr(Mitos2Runner, "run", lambda self: _fake_result(success = False, return_code = 1))
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, config_json))

    assert run_mitos2_annotation.main() == 1
    result = json.loads((output_dir / f"{sample_id}.mitos2_annotation_result.json").read_text())
    assert result["success"] is False
    assert result["feature_count"] is None

def test_validation_failure_systemexit_is_caught_not_raised(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    sample_id = "SAMPLE_04"
    output_dir = tmp_path / "annotation" / sample_id
    config_json = _write_config(tmp_path, sample_id, output_dir)

    def _raise(self):
        raise SystemExit(1)
    monkeypatch.setattr(Mitos2Runner, "run", _raise)
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, config_json))

    assert run_mitos2_annotation.main() == 1   # must not raise

def test_malformed_gff_after_success_is_recorded_not_a_crash(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    # MITOS2 itself exits 0, but its GFF is somehow malformed (e.g. a
    # gene row missing its Name qualifier) -- must not crash main(),
    # must be recorded as a real failure, not silently accepted.
    sample_id = "SAMPLE_05"
    output_dir = tmp_path / "annotation" / sample_id
    config_json = _write_config(tmp_path, sample_id, output_dir)

    def _stub_run(self: Mitos2Runner) -> CommandResult:
        self.result_gff.parent.mkdir(parents = True, exist_ok = True)
        self.result_gff.write_text("test\tmitos\tgene\t1\t600\t.\t+\t.\tID=gene_cox2\n")   # no Name/gene_id
        self.result_fas.write_text(">x\nACGT\n")
        return _fake_result()

    monkeypatch.setattr(Mitos2Runner, "run", _stub_run)
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, config_json))

    assert run_mitos2_annotation.main() == 1
    result = json.loads((output_dir / f"{sample_id}.mitos2_annotation_result.json").read_text())
    assert result["success"] is False
    assert "gff parse error" in result["reason"]

def test_missing_config_json_is_recorded_not_a_crash(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, "SAMPLE_06", tmp_path / "does_not_exist.json"))
    assert run_mitos2_annotation.main() == 1
