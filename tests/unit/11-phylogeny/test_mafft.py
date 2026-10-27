"""
test_mafft.py

Unit tests for the MafftRunner class (#346): input FASTA, threads,
output alignment, and exit code.
"""

# Imports
import logging
from pathlib import Path
import pytest
from mitopipeline.api.mafft import MafftRunner

def _make_runner(tmp_path: Path, **overrides) -> MafftRunner:
    input_fasta = tmp_path / "combined.fasta"
    input_fasta.write_text(">a\nACGT\n>b\nACGT\n")
    defaults = dict(
        working_dir = tmp_path,
        input_fasta = input_fasta,
        output_fasta = tmp_path / "aligned.fasta",
        sample_id = "test_run",
        logger = logging.getLogger("test"),
    )
    defaults.update(overrides)
    return MafftRunner(**defaults)

def test_build_command_uses_auto_strategy_by_default(tmp_path: Path) -> None:
    runner = _make_runner(tmp_path)
    command = runner.build_command()
    assert command[0] == "bash"
    assert command[1] == "-c"
    assert "mafft" in command[2]
    assert "--auto" in command[2]
    assert "--nuc" in command[2]

def test_build_command_redirects_to_the_real_output_path(tmp_path: Path) -> None:
    runner = _make_runner(tmp_path)
    command = runner.build_command()
    assert f"> {runner.output_fasta}" in command[2] or str(runner.output_fasta) in command[2]

def test_build_command_honors_custom_strategy(tmp_path: Path) -> None:
    runner = _make_runner(tmp_path, strategy = "--localpair")
    command = runner.build_command()
    assert "--localpair" in command[2]
    assert "--auto" not in command[2]

def test_validate_inputs_rejects_missing_fasta(tmp_path: Path) -> None:
    runner = _make_runner(tmp_path)
    runner.input_fasta.unlink()
    with pytest.raises(SystemExit):
        runner.validate_inputs()

def test_validate_inputs_rejects_empty_fasta(tmp_path: Path) -> None:
    runner = _make_runner(tmp_path)
    runner.input_fasta.write_text("")
    with pytest.raises(SystemExit):
        runner.validate_inputs()

def test_validate_inputs_rejects_invalid_threads(tmp_path: Path) -> None:
    runner = _make_runner(tmp_path, threads = 0)
    with pytest.raises(SystemExit):
        runner.validate_inputs()

def test_validate_inputs_passes_for_valid_configuration(tmp_path: Path) -> None:
    runner = _make_runner(tmp_path)
    runner.validate_inputs()   # does not raise

def test_validate_outputs_rejects_missing_alignment(tmp_path: Path) -> None:
    runner = _make_runner(tmp_path)
    with pytest.raises(SystemExit):
        runner.validate_outputs()

def test_validate_outputs_rejects_empty_alignment(tmp_path: Path) -> None:
    runner = _make_runner(tmp_path)
    runner.output_fasta.parent.mkdir(parents = True, exist_ok = True)
    runner.output_fasta.touch()
    with pytest.raises(SystemExit):
        runner.validate_outputs()

def test_validate_outputs_passes_with_real_content(tmp_path: Path) -> None:
    runner = _make_runner(tmp_path)
    runner.output_fasta.parent.mkdir(parents = True, exist_ok = True)
    runner.output_fasta.write_text(">a\nACGT\n>b\nACGT\n")
    runner.validate_outputs()   # does not raise
