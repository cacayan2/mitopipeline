"""
test_iqtree.py

Unit tests for the IqtreeRunner class (#347): alignment input, model,
bootstrap parameters, output tree, and failure handling.
"""

# Imports
import logging
from pathlib import Path
import pytest
from mitopipeline.api.iqtree import IqtreeRunner

def _make_runner(tmp_path: Path, **overrides) -> IqtreeRunner:
    alignment = tmp_path / "aligned.fasta"
    alignment.write_text(">a\nACGT\n>b\nACGT\n")
    defaults = dict(
        working_dir = tmp_path,
        alignment_fasta = alignment,
        output_prefix = tmp_path / "tree" / "test_run",
        sample_id = "test_run",
        logger = logging.getLogger("test"),
    )
    defaults.update(overrides)
    return IqtreeRunner(**defaults)

def test_build_command_includes_real_flags(tmp_path: Path) -> None:
    runner = _make_runner(tmp_path)
    command = runner.build_command()
    assert command[0] == "iqtree3"
    assert "-s" in command and command[command.index("-s") + 1] == str(runner.alignment_fasta)
    assert "-m" in command and command[command.index("-m") + 1] == "MFP"
    assert "-B" in command and command[command.index("-B") + 1] == "1000"
    assert "--alrt" in command and command[command.index("--alrt") + 1] == "1000"
    assert "--prefix" in command and command[command.index("--prefix") + 1] == str(runner.output_prefix)

def test_treefile_path_is_derived_from_prefix(tmp_path: Path) -> None:
    runner = _make_runner(tmp_path)
    assert runner.treefile == Path(f"{runner.output_prefix}.treefile")

def test_validate_inputs_rejects_missing_alignment(tmp_path: Path) -> None:
    runner = _make_runner(tmp_path)
    runner.alignment_fasta.unlink()
    with pytest.raises(SystemExit):
        runner.validate_inputs()

def test_validate_inputs_rejects_low_bootstrap_replicates(tmp_path: Path) -> None:
    runner = _make_runner(tmp_path, ultrafast_bootstrap = 100)
    with pytest.raises(SystemExit):
        runner.validate_inputs()

def test_validate_inputs_rejects_low_sh_alrt_replicates(tmp_path: Path) -> None:
    runner = _make_runner(tmp_path, sh_alrt = 100)
    with pytest.raises(SystemExit):
        runner.validate_inputs()

def test_validate_inputs_rejects_invalid_threads(tmp_path: Path) -> None:
    runner = _make_runner(tmp_path, threads = 0)
    with pytest.raises(SystemExit):
        runner.validate_inputs()

def test_validate_inputs_passes_for_valid_configuration(tmp_path: Path) -> None:
    runner = _make_runner(tmp_path)
    runner.validate_inputs()   # does not raise

def test_validate_outputs_rejects_missing_treefile(tmp_path: Path) -> None:
    runner = _make_runner(tmp_path)
    runner.output_prefix.parent.mkdir(parents = True, exist_ok = True)
    with pytest.raises(SystemExit):
        runner.validate_outputs()

def test_validate_outputs_rejects_empty_treefile(tmp_path: Path) -> None:
    runner = _make_runner(tmp_path)
    runner.output_prefix.parent.mkdir(parents = True, exist_ok = True)
    runner.treefile.touch()
    with pytest.raises(SystemExit):
        runner.validate_outputs()

def test_validate_outputs_passes_with_real_content(tmp_path: Path) -> None:
    runner = _make_runner(tmp_path)
    runner.output_prefix.parent.mkdir(parents = True, exist_ok = True)
    runner.treefile.write_text("(a:0.1,b:0.1);\n")
    runner.validate_outputs()   # does not raise
