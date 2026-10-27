"""
test_run_tree_figure.py

Fixture tests for the job-level tree figure rendering exec: reads
run_iqtree.py's own result JSON and renders a PNG when a usable tree
exists -- a missing/unsuccessful IQ-TREE result is a normal outcome
this step must still complete for, not a failure of its own.
"""

# Imports
import importlib
import json
import sys
from pathlib import Path
import pytest

run_tree_figure = importlib.import_module("mitopipeline.exec.11-phylogeny.run_tree_figure")

def _argv(tmp_path: Path, iqtree_result_json: Path) -> list[str]:
    return [
        "run_tree_figure.py",
        "--iqtree-result-json", str(iqtree_result_json),
        "--dpi", "600",
        "--midpoint-root",
        "--output-dir", str(tmp_path / "out"),
        "--log-file", str(tmp_path / "logs" / "tree_figure.log"),
        "--global-log-file", str(tmp_path / "logs" / "global.log"),
    ]

def test_main_renders_a_figure_when_a_usable_tree_exists(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    newick_path = tmp_path / "tree.nwk"
    newick_path.write_text("((A:0.1,B:0.2):0.05,C:0.3);")
    result_json = tmp_path / "iqtree_result.json"
    result_json.write_text(json.dumps({"success": True, "tip_count": 3, "newick_path": str(newick_path), "reason": "completed successfully"}))

    monkeypatch.setattr(sys, "argv", _argv(tmp_path, result_json))

    assert run_tree_figure.main() == 0
    result = json.loads((tmp_path / "out" / "tree_figure_result.json").read_text())
    assert result["success"] is True
    assert Path(result["tree_figure_png"]).is_file()

def test_main_skips_rendering_when_iqtree_did_not_succeed(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    result_json = tmp_path / "iqtree_result.json"
    result_json.write_text(json.dumps({"success": False, "tip_count": None, "newick_path": None, "reason": "not enough taxa"}))

    monkeypatch.setattr(sys, "argv", _argv(tmp_path, result_json))

    assert run_tree_figure.main() == 0
    result = json.loads((tmp_path / "out" / "tree_figure_result.json").read_text())
    assert result["success"] is False
    assert result["tree_figure_png"] is None

def test_main_skips_rendering_when_the_result_json_is_missing(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, tmp_path / "does_not_exist.json"))

    assert run_tree_figure.main() == 0
    result = json.loads((tmp_path / "out" / "tree_figure_result.json").read_text())
    assert result["success"] is False
