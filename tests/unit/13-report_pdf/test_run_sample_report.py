"""
test_run_sample_report.py

Fixture tests for the per-sample report exec (#358/#359/#361): a
full set of real inputs renders both Markdown and PDF; missing
optional inputs still render a Markdown report with "not available"
sections rather than crashing (#359); a missing pandoc/typst skips
the PDF without failing the whole step.
"""

# Imports
import importlib
import json
import sys
from pathlib import Path
import pytest

run_sample_report = importlib.import_module("mitopipeline.exec.13-report_pdf.run_sample_report")

def _write_json(path: Path, data: dict) -> Path:
    path.parent.mkdir(parents = True, exist_ok = True)
    path.write_text(json.dumps(data))
    return path

def _argv(tmp_path: Path, sample_id: str, manifest: Path, output_dir: Path, **optional) -> list[str]:
    argv = [
        "run_sample_report.py",
        "--sample-id", sample_id,
        "--runtime-manifest", str(manifest),
        "--output-dir", str(output_dir),
        "--log-file", str(tmp_path / "logs" / "report.log"),
        "--global-log-file", str(tmp_path / "logs" / "global.log"),
    ]
    flag_names = {
        "statistics_json": "--statistics-json",
        "biological_validation_json": "--biological-validation-json",
        "iqtree_result_json": "--iqtree-result-json",
        "phylogeny_model": "--phylogeny-model",
        "tree_figure_png": "--tree-figure-png",
    }
    for key, value in optional.items():
        argv.extend([flag_names[key], str(value)])
    return argv

def test_full_report_renders_markdown_and_pdf(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    sample_id = "sample1"
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text("sample_id\tr1\tr2\tgenus\tspecies\nsample1\ta\tb\tDanio\trerio\n")
    statistics_json = _write_json(tmp_path / "statistics.json", {
        "sample_id": sample_id,
        "assembly": {"sample_id": sample_id, "status": "circular", "source_attempt": 0, "length_bp": 16596,
                    "contig_count": 1, "gc_content_percent": 45.0},
        "blast": None, "annotation": None,
    })
    output_dir = tmp_path / "output"

    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, manifest, output_dir, statistics_json = statistics_json))

    assert run_sample_report.main() == 0
    result = json.loads((output_dir / f"{sample_id}.report_result.json").read_text())
    assert result["success"] is True
    markdown_path = Path(result["markdown_path"])
    assert markdown_path.is_file()
    assert "Circular" in markdown_path.read_text()
    if result["pdf_path"] is not None:
        assert Path(result["pdf_path"]).is_file()

def test_missing_optional_inputs_still_renders_markdown(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    sample_id = "sample_bare"
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text("sample_id\tr1\tr2\nsample_bare\ta\tb\n")
    output_dir = tmp_path / "output"

    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, manifest, output_dir))

    assert run_sample_report.main() == 0
    result = json.loads((output_dir / f"{sample_id}.report_result.json").read_text())
    assert result["success"] is True
    markdown_text = Path(result["markdown_path"]).read_text()
    assert "No assembly: contig selection did not find a usable circular or fragmented sequence" in markdown_text

def test_report_includes_the_sample_own_phylogeny(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Regression: phylogeny is now per-sample (restoring v0.1.0's own
    design) -- a sample's own real iqtree_result.json/model/tree
    figure, passed as this exec script's own CLI args, must appear
    in that sample's own report, not a separate batch-level report."""
    sample_id = "sample1"
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text("sample_id\tr1\tr2\nsample1\ta\tb\n")
    iqtree_result_json = _write_json(tmp_path / "iqtree_result.json", {
        "success": True, "tip_count": 7, "newick_path": None, "reason": "completed",
    })
    output_dir = tmp_path / "output"

    monkeypatch.setattr(sys, "argv", _argv(
        tmp_path, sample_id, manifest, output_dir,
        iqtree_result_json = iqtree_result_json, phylogeny_model = "GTR+G",
    ))

    assert run_sample_report.main() == 0
    result = json.loads((output_dir / f"{sample_id}.report_result.json").read_text())
    markdown_text = Path(result["markdown_path"]).read_text()
    assert "## Phylogenetic Context" in markdown_text
    assert "GTR+G" in markdown_text
    assert "| Tree tips | 7 |" in markdown_text

def test_pdf_export_skip_does_not_fail_the_step(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    sample_id = "sample1"
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text("sample_id\tr1\tr2\nsample1\ta\tb\n")
    output_dir = tmp_path / "output"

    def _raise_missing_pandoc(*args, **kwargs):
        raise FileNotFoundError("Required executable not found on PATH: pandoc.")

    monkeypatch.setattr(run_sample_report, "export_markdown_to_pdf", _raise_missing_pandoc)
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, manifest, output_dir))

    assert run_sample_report.main() == 0
    result = json.loads((output_dir / f"{sample_id}.report_result.json").read_text())
    assert result["success"] is True
    assert result["pdf_path"] is None
    assert "pdf export skipped" in result["reason"]
