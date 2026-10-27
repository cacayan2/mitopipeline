"""
test_run_build_combined_fasta.py

Fixture tests for the combined-FASTA exec (#344): reads every
sample's own contig_selection/biological_validation records and
builds the batch-level dataset.
"""

# Imports
import importlib
import json
import sys
from pathlib import Path
import pytest

run_build_combined_fasta = importlib.import_module("mitopipeline.exec.10-alignment_dataset.run_build_combined_fasta")

def _write_sample(contig_dir: Path, sample_id: str, verdict: str = "PASS", selected: bool = True) -> None:
    sample_dir = contig_dir / sample_id
    sample_dir.mkdir(parents = True, exist_ok = True)
    fasta_path = sample_dir / "data.fasta"
    fasta_path.write_text(f">{sample_id}_native\nACGT\n")
    (sample_dir / f"{sample_id}.contig_selection.json").write_text(json.dumps({
        "sample_id": sample_id, "selected": selected, "fasta_path": str(fasta_path),
    }))
    (sample_dir / f"{sample_id}.biological_validation.json").write_text(json.dumps({
        "sample_id": sample_id, "selected_attempt": 0,
        "results": [{"attempt": 0, "overall_verdict": verdict}],
    }))

def _argv(tmp_path: Path, sample_ids: list[str], contig_dir: Path) -> list[str]:
    argv = ["run_build_combined_fasta.py"]
    for sid in sample_ids:
        argv += ["--sample-id", sid]
    argv += [
        "--contig-selection-dir", str(contig_dir),
        "--output-dir", str(tmp_path / "out"),
        "--log-file", str(tmp_path / "logs" / "combined.log"),
        "--global-log-file", str(tmp_path / "logs" / "global.log"),
    ]
    return argv

def test_accepted_samples_are_combined(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    contig_dir = tmp_path / "06-contig_selection"
    _write_sample(contig_dir, "sample_a", verdict = "PASS")
    _write_sample(contig_dir, "sample_b", verdict = "REVIEW")
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, ["sample_a", "sample_b"], contig_dir))

    assert run_build_combined_fasta.main() == 0
    result = json.loads((tmp_path / "out" / "combined_fasta_result.json").read_text())
    assert result["success"] is True
    assert set(result["included_sample_ids"]) == {"sample_a", "sample_b"}

def test_failed_validation_excludes_the_sample(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    contig_dir = tmp_path / "06-contig_selection"
    _write_sample(contig_dir, "sample_a", verdict = "PASS")
    _write_sample(contig_dir, "sample_bad", verdict = "FAIL")
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, ["sample_a", "sample_bad"], contig_dir))

    assert run_build_combined_fasta.main() == 0
    result = json.loads((tmp_path / "out" / "combined_fasta_result.json").read_text())
    assert result["included_sample_ids"] == ["sample_a"]
    assert "sample_bad" in result["excluded_sample_ids"]

def test_unselected_sample_is_excluded(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    contig_dir = tmp_path / "06-contig_selection"
    _write_sample(contig_dir, "sample_a", verdict = "PASS")
    _write_sample(contig_dir, "sample_none", selected = False)
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, ["sample_a", "sample_none"], contig_dir))

    assert run_build_combined_fasta.main() == 0
    result = json.loads((tmp_path / "out" / "combined_fasta_result.json").read_text())
    assert "sample_none" in result["excluded_sample_ids"]

def test_missing_sample_records_is_recorded_not_a_crash(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    contig_dir = tmp_path / "06-contig_selection"
    _write_sample(contig_dir, "sample_a", verdict = "PASS")
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, ["sample_a", "sample_never_ran"], contig_dir))

    assert run_build_combined_fasta.main() == 0
    result = json.loads((tmp_path / "out" / "combined_fasta_result.json").read_text())
    assert "sample_never_ran" in result["excluded_sample_ids"]

def test_zero_accepted_samples_is_recorded_not_a_crash(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    contig_dir = tmp_path / "06-contig_selection"
    _write_sample(contig_dir, "sample_bad", verdict = "FAIL")
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, ["sample_bad"], contig_dir))

    assert run_build_combined_fasta.main() == 1
    result = json.loads((tmp_path / "out" / "combined_fasta_result.json").read_text())
    assert result["success"] is False
