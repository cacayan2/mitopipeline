"""
test_run_iqtree.py

Fixture tests for the IQ-TREE exec (#436): runs through BaseTool and
validates the real tree, including missing and malformed tree
artifacts (#437).
"""

# Imports
import importlib
import json
import sys
from datetime import datetime
from pathlib import Path
import pytest
from mitopipeline.api.iqtree import IqtreeRunner
from mitopipeline.models.command_result.command_result import CommandResult

run_iqtree = importlib.import_module("mitopipeline.exec.11-phylogeny.run_iqtree")

def _fake_result(success: bool = True, return_code: int = 0) -> CommandResult:
    now = datetime.now()
    return CommandResult(command = [], return_code = return_code, stdout = "", stderr = "", runtime_seconds = 1.0,
                         success = success, tool_name = "iqtree", started_at = now, ended_at = now,
                         peak_memory_mb = 0.0, user_cpu_seconds = 0.0, system_cpu_seconds = 0.0)

def _argv(tmp_path: Path, alignment: Path, output_prefix: Path, output_newick: Path) -> list[str]:
    return [
        "run_iqtree.py", "--run-id", "test_run",
        "--alignment", str(alignment), "--output-prefix", str(output_prefix), "--output-newick", str(output_newick),
        "--threads", "2", "--working-dir", str(tmp_path), "--output-dir", str(tmp_path / "out"),
        "--log-file", str(tmp_path / "logs" / "iqtree.log"), "--global-log-file", str(tmp_path / "logs" / "global.log"),
    ]

def _write_alignment(path: Path, ids: list[str]) -> None:
    path.write_text("".join(f">{i}\nACGTACGT\n" for i in ids))

def test_successful_run_with_valid_tree(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    alignment = tmp_path / "aligned.fasta"
    _write_alignment(alignment, ["a", "b", "c"])
    prefix = tmp_path / "tree" / "run"
    newick = tmp_path / "tree.nwk"

    def _stub_run(self: IqtreeRunner) -> CommandResult:
        self.treefile.parent.mkdir(parents = True, exist_ok = True)
        self.treefile.write_text("(a:0.1,(b:0.1,c:0.1):0.1);\n")
        return _fake_result()

    monkeypatch.setattr(IqtreeRunner, "run", _stub_run)
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, alignment, prefix, newick))

    assert run_iqtree.main() == 0
    result = json.loads((tmp_path / "out" / "iqtree_result.json").read_text())
    assert result["success"] is True
    assert result["tip_count"] == 3
    assert newick.is_file()

def test_malformed_newick_tree_is_rejected(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    alignment = tmp_path / "aligned.fasta"
    _write_alignment(alignment, ["a", "b", "c"])
    prefix = tmp_path / "tree" / "run"

    def _stub_run(self: IqtreeRunner) -> CommandResult:
        self.treefile.parent.mkdir(parents = True, exist_ok = True)
        self.treefile.write_text("this is not a newick tree at all (((\n")
        return _fake_result()

    monkeypatch.setattr(IqtreeRunner, "run", _stub_run)
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, alignment, prefix, tmp_path / "tree.nwk"))

    assert run_iqtree.main() == 1
    result = json.loads((tmp_path / "out" / "iqtree_result.json").read_text())
    assert result["success"] is False
    assert "malformed newick" in result["reason"]

def test_tip_count_mismatch_is_rejected(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    alignment = tmp_path / "aligned.fasta"
    _write_alignment(alignment, ["a", "b", "c"])   # 3 sequences
    prefix = tmp_path / "tree" / "run"

    def _stub_run(self: IqtreeRunner) -> CommandResult:
        self.treefile.parent.mkdir(parents = True, exist_ok = True)
        self.treefile.write_text("(a:0.1,b:0.1);\n")   # only 2 tips
        return _fake_result()

    monkeypatch.setattr(IqtreeRunner, "run", _stub_run)
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, alignment, prefix, tmp_path / "tree.nwk"))

    assert run_iqtree.main() == 1
    result = json.loads((tmp_path / "out" / "iqtree_result.json").read_text())
    assert "tip count" in result["reason"]

def test_missing_treefile_is_recorded_not_a_crash(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    alignment = tmp_path / "aligned.fasta"
    _write_alignment(alignment, ["a", "b", "c"])
    prefix = tmp_path / "tree" / "run"

    def _raise(self):
        raise SystemExit(1)
    monkeypatch.setattr(IqtreeRunner, "run", _raise)
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, alignment, prefix, tmp_path / "tree.nwk"))

    assert run_iqtree.main() == 1   # must not raise

def test_tool_failure_is_recorded_not_fabricated(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    alignment = tmp_path / "aligned.fasta"
    _write_alignment(alignment, ["a", "b", "c"])
    monkeypatch.setattr(IqtreeRunner, "run", lambda self: _fake_result(success = False, return_code = 1))
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, alignment, tmp_path / "tree" / "run", tmp_path / "tree.nwk"))

    assert run_iqtree.main() == 1
    result = json.loads((tmp_path / "out" / "iqtree_result.json").read_text())
    assert result["success"] is False

def test_missing_alignment_is_recorded_not_a_crash(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, tmp_path / "does_not_exist.fasta", tmp_path / "tree" / "run", tmp_path / "tree.nwk"))
    assert run_iqtree.main() == 1

def test_too_few_sequences_is_skipped_with_a_readable_reason_not_invoked(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Regression: IQ-TREE's own real minimum is 3 sequences; running
    it with fewer used to surface as a cryptic "iqtree exited 2" in
    the report. The tool must not even be invoked in this case."""
    alignment = tmp_path / "aligned.fasta"
    _write_alignment(alignment, ["a", "b"])

    def _fail_if_called(self):
        raise AssertionError("IqtreeRunner.run() must not be called when too few sequences are present.")
    monkeypatch.setattr(IqtreeRunner, "run", _fail_if_called)
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, alignment, tmp_path / "tree" / "run", tmp_path / "tree.nwk"))

    assert run_iqtree.main() == 0
    result = json.loads((tmp_path / "out" / "iqtree_result.json").read_text())
    assert result["success"] is False
    assert result["tip_count"] == 2
    assert "at least 3 assemblies are required" in result["reason"]
    assert "but only 2 passed validation" in result["reason"]
