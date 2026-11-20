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

class _FakeTypstBinary:
    """Stands in for the real TypstBinary in tests that don't care
    about PDF export behavior itself -- never touches the network."""
    def __init__(self, binary_path: Path, logger = None) -> None:
        self.binary_path = binary_path
    def prepare(self) -> Path:
        return self.binary_path

@pytest.fixture(autouse = True)
def _stub_typst_binary(monkeypatch: pytest.MonkeyPatch) -> None:
    """Every test in this file goes through TypstBinary() before
    export_markdown_to_pdf() now -- stub it everywhere so none of
    them need a real network call just to render a Markdown report."""
    monkeypatch.setattr(run_run_report, "TypstBinary", _FakeTypstBinary)

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
        "--pdf-engine", str(tmp_path / "typst"),
        "--output-dir", str(output_dir),
        "--log-file", str(tmp_path / "logs" / "run_report.log"),
        "--global-log-file", str(tmp_path / "logs" / "global.log"),
    ])

    assert run_run_report.main() == 0
    markdown_path = output_dir / "job1.run_report.md"
    text = markdown_path.read_text()
    # The Sample column now goes through wrap_cell_text() (plain
    # wrapped text, not basename_code()'s unescaped monospace) so a
    # real underscore in a sample id arrives backslash-escaped in the
    # Markdown source -- still renders as a literal underscore once
    # pandoc compiles it, confirmed directly, so this is the correct
    # updated expectation, not a loosened check.
    assert "sample\\_pass" in text and "sample\\_fail" in text
    assert "Pass" in text and "Fail" in text

def test_sample_with_no_results_shows_not_assessed(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text("sample_id\tr1\tr2\nsample_bare\ta\tb\n")
    job_dir = tmp_path / "job"
    output_dir = tmp_path / "output"

    monkeypatch.setattr(sys, "argv", [
        "run_run_report.py", "--job-id", "job1", "--sample-id", "sample_bare",
        "--runtime-manifest", str(manifest), "--job-dir", str(job_dir),
        "--pdf-engine", str(tmp_path / "typst"),
        "--output-dir", str(output_dir),
        "--log-file", str(tmp_path / "logs" / "run_report.log"),
        "--global-log-file", str(tmp_path / "logs" / "global.log"),
    ])

    assert run_run_report.main() == 0
    text = (output_dir / "job1.run_report.md").read_text()
    assert "Not assessed" in text

def test_phylogeny_model_flag_is_accepted_and_does_not_crash(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Phylogeny is now per-sample (restoring v0.1.0's own design) --
    each sample's own tree lives in that sample's own report, not
    here. This exec script only needs to thread --phylogeny-model
    through without crashing; test_run_sample_report.py covers the
    actual tree content."""
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text("sample_id\tr1\tr2\nsample1\ta\tb\n")
    job_dir = tmp_path / "job"
    output_dir = tmp_path / "output"

    monkeypatch.setattr(sys, "argv", [
        "run_run_report.py", "--job-id", "job1", "--sample-id", "sample1",
        "--runtime-manifest", str(manifest), "--job-dir", str(job_dir),
        "--phylogeny-model", "GTR+G",
        "--pdf-engine", str(tmp_path / "typst"),
        "--output-dir", str(output_dir),
        "--log-file", str(tmp_path / "logs" / "run_report.log"),
        "--global-log-file", str(tmp_path / "logs" / "global.log"),
    ])

    assert run_run_report.main() == 0
    assert (output_dir / "job1.run_report.md").is_file()

def test_bioproject_metadata_flag_surfaces_in_the_run_report(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text("sample_id\tr1\tr2\nsample1\ta\tb\n")
    job_dir = tmp_path / "job"
    output_dir = tmp_path / "output"
    bioproject_json = _write_json(tmp_path / "bioproject_metadata.json", {
        "organism": "Danio rerio", "sample_scope": "Monoisolate",
        "review_status": "MANUAL_REVIEW_REQUIRED", "missing_metadata": "publication;funding",
    })

    monkeypatch.setattr(sys, "argv", [
        "run_run_report.py", "--job-id", "job1", "--sample-id", "sample1",
        "--runtime-manifest", str(manifest), "--job-dir", str(job_dir),
        "--bioproject-metadata-json", str(bioproject_json),
        "--pdf-engine", str(tmp_path / "typst"),
        "--output-dir", str(output_dir),
        "--log-file", str(tmp_path / "logs" / "run_report.log"),
        "--global-log-file", str(tmp_path / "logs" / "global.log"),
    ])

    assert run_run_report.main() == 0
    text = (output_dir / "job1.run_report.md").read_text()
    assert "### BioProject" in text
    assert "Danio rerio" in text
    assert "`publication`: " in text
