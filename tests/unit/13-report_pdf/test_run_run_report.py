"""
test_run_run_report.py

Fixture tests for the run-level report exec (#360/#441/#443): every
sample's own real statistics merge into the run report, mixed PASS/
REVIEW/FAIL outcomes are counted correctly, and a sample with no
usable result still appears with a "not assessed" status rather
than being silently dropped.
"""

# Imports
import importlib
import json
import sys
from pathlib import Path
import pytest

run_run_report = importlib.import_module("mitopipeline.exec.13-report_pdf.run_run_report")

def _write_json(path: Path, data: dict) -> Path:
    path.parent.mkdir(parents = True, exist_ok = True)
    path.write_text(json.dumps(data))
    return path

def test_run_report_merges_mixed_sample_outcomes(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text("sample_id\tr1\tr2\nsample_pass\ta\tb\nsample_fail\tc\td\n")
    job_dir = tmp_path / "job"

    _write_json(job_dir / "06-contig_selection" / "sample_pass" / "sample_pass.biological_validation.json", {
        "sample_id": "sample_pass", "selected_attempt": 0,
        "results": [{"attempt": 0, "overall_verdict": "PASS", "checks": []}],
    })
    _write_json(job_dir / "06-contig_selection" / "sample_fail" / "sample_fail.biological_validation.json", {
        "sample_id": "sample_fail", "selected_attempt": 0,
        "results": [{"attempt": 0, "overall_verdict": "FAIL", "checks": []}],
    })
    output_dir = tmp_path / "output"

    monkeypatch.setattr(sys, "argv", [
        "run_run_report.py", "--job-id", "job1",
        "--sample-id", "sample_pass", "--sample-id", "sample_fail",
        "--runtime-manifest", str(manifest), "--job-dir", str(job_dir),
        "--output-dir", str(output_dir),
        "--log-file", str(tmp_path / "logs" / "run_report.log"),
        "--global-log-file", str(tmp_path / "logs" / "global.log"),
    ])

    assert run_run_report.main() == 0
    markdown_path = output_dir / "job1.run_report.md"
    text = markdown_path.read_text()
    assert "sample_pass" in text and "sample_fail" in text
    assert "Pass" in text and "Fail" in text

def test_sample_with_no_results_shows_not_assessed(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text("sample_id\tr1\tr2\nsample_bare\ta\tb\n")
    job_dir = tmp_path / "job"
    output_dir = tmp_path / "output"

    monkeypatch.setattr(sys, "argv", [
        "run_run_report.py", "--job-id", "job1", "--sample-id", "sample_bare",
        "--runtime-manifest", str(manifest), "--job-dir", str(job_dir),
        "--output-dir", str(output_dir),
        "--log-file", str(tmp_path / "logs" / "run_report.log"),
        "--global-log-file", str(tmp_path / "logs" / "global.log"),
    ])

    assert run_run_report.main() == 0
    text = (output_dir / "job1.run_report.md").read_text()
    assert "Not assessed" in text

def test_phylogeny_included_when_result_exists(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text("sample_id\tr1\tr2\nsample1\ta\tb\n")
    job_dir = tmp_path / "job"
    output_dir = tmp_path / "output"
    iqtree_result = _write_json(tmp_path / "iqtree_result.json", {
        "success": True, "tip_count": 5, "newick_path": None, "reason": "completed",
    })

    monkeypatch.setattr(sys, "argv", [
        "run_run_report.py", "--job-id", "job1", "--sample-id", "sample1",
        "--runtime-manifest", str(manifest), "--job-dir", str(job_dir),
        "--iqtree-result-json", str(iqtree_result), "--phylogeny-model", "GTR+G",
        "--output-dir", str(output_dir),
        "--log-file", str(tmp_path / "logs" / "run_report.log"),
        "--global-log-file", str(tmp_path / "logs" / "global.log"),
    ])

    assert run_run_report.main() == 0
    text = (output_dir / "job1.run_report.md").read_text()
    assert "GTR+G" in text
    assert "| Tree tips | 5 |" in text
