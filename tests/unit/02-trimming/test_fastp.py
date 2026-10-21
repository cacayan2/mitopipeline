"""
test_fastp.py

Unit tests for the FastpRunner class.
"""

# Imports
import logging
from pathlib import Path
import pytest
from mitopipeline.api.fastp import FastpRunner
from mitopipeline.models.sample.sample import Sample

@pytest.fixture
def paired_sample(tmp_path: Path) -> Sample:
    """
    Creates a Sample backed by two real (empty) FASTQ files, so input
    validation has something to find on disk.

    Args:
        tmp_path (Path): The pytest-provided temporary directory.

    Returns:
        Sample: The constructed sample.
    """
    r1 = tmp_path / "raw_R1.fastq.gz"
    r2 = tmp_path / "raw_R2.fastq.gz"
    r1.touch()
    r2.touch()
    return Sample(sample_id = "SAMPLE_01", r1 = r1, r2 = r2)

def _make_runner(tmp_path: Path, sample: Sample, **overrides) -> FastpRunner:
    """
    Builds a FastpRunner with sensible defaults, overridable per test.

    Args:
        tmp_path (Path): The pytest-provided temporary directory.
        sample (Sample): The sample to run fastp on.
        **overrides: Keyword arguments overriding the defaults.

    Returns:
        FastpRunner: The constructed runner.
    """
    defaults = dict(
        working_dir = tmp_path,
        output_dir = tmp_path / "trimmed",
        sample = sample,
        logger = logging.getLogger("test"),
    )
    defaults.update(overrides)
    return FastpRunner(**defaults)

def test_build_command_includes_required_flags(tmp_path: Path, paired_sample: Sample) -> None:
    """
    build_command should wire up fastp's required I/O flags and the
    sample-scoped output filenames, with no stray/duplicate flags.
    """
    runner = _make_runner(tmp_path, paired_sample, threads = 8)
    command = runner.build_command()
    assert command[0] == "fastp"
    assert command.count("--threads") == 0
    assert "--thread" in command
    assert command[command.index("--thread") + 1] == "8"
    assert "--in1" in command and command[command.index("--in1") + 1] == str(paired_sample.r1)
    assert "--in2" in command and command[command.index("--in2") + 1] == str(paired_sample.r2)
    assert str(runner.output_r1).endswith("SAMPLE_01.R1.trimmed.fastq.gz")
    assert str(runner.output_r2).endswith("SAMPLE_01.R2.trimmed.fastq.gz")

def test_build_command_creates_output_directory(tmp_path: Path, paired_sample: Sample) -> None:
    """
    build_command should create the output directory if it does not exist.
    """
    runner = _make_runner(tmp_path, paired_sample)
    assert not runner.output_dir.exists()
    runner.build_command()
    assert runner.output_dir.is_dir()

def test_build_command_applies_value_and_boolean_options(tmp_path: Path, paired_sample: Sample) -> None:
    """
    Tool options should be appended for provided values, and boolean
    flags only appended when explicitly True.
    """
    runner = _make_runner(tmp_path, paired_sample, tool_options = {
        "length_required": 75,
        "trim_poly_g": True,
        "cut_front": False,
    })
    command = runner.build_command()
    assert "--length_required" in command
    assert command[command.index("--length_required") + 1] == "75"
    assert "--trim_poly_g" in command
    assert "--cut_front" not in command

def test_validate_inputs_exits_on_missing_read_file(tmp_path: Path) -> None:
    """
    A missing R1/R2 file should exit rather than let fastp fail later.
    """
    sample = Sample(sample_id = "SAMPLE_02", r1 = tmp_path / "missing_R1.fastq.gz", r2 = tmp_path / "missing_R2.fastq.gz")
    runner = _make_runner(tmp_path, sample)
    with pytest.raises(SystemExit):
        runner.validate_inputs()

def test_validate_inputs_exits_on_zero_threads(tmp_path: Path, paired_sample: Sample) -> None:
    """
    A non-positive thread count should exit.
    """
    runner = _make_runner(tmp_path, paired_sample, threads = 0)
    with pytest.raises(SystemExit):
        runner.validate_inputs()

def test_validate_inputs_passes_for_valid_sample(tmp_path: Path, paired_sample: Sample) -> None:
    """
    A sample with real input files and a valid thread count should pass
    validation without exiting.
    """
    runner = _make_runner(tmp_path, paired_sample, threads = 1)
    runner.validate_inputs()

def test_validate_outputs_exits_on_missing_output(tmp_path: Path, paired_sample: Sample) -> None:
    """
    validate_outputs should exit if any expected output file is missing.
    """
    runner = _make_runner(tmp_path, paired_sample)
    runner.output_dir.mkdir(parents = True, exist_ok = True)
    with pytest.raises(SystemExit):
        runner.validate_outputs()

def test_validate_outputs_exits_on_empty_output(tmp_path: Path, paired_sample: Sample) -> None:
    """
    validate_outputs should exit if an expected output file exists but is
    empty, since an empty file means fastp did not actually produce it.
    """
    runner = _make_runner(tmp_path, paired_sample)
    runner.output_dir.mkdir(parents = True, exist_ok = True)
    for path in (runner.output_r1, runner.output_r2, runner.output_html, runner.output_json):
        path.touch()
    with pytest.raises(SystemExit):
        runner.validate_outputs()

def test_validate_outputs_passes_when_all_outputs_present(tmp_path: Path, paired_sample: Sample) -> None:
    """
    validate_outputs should pass once every expected output file exists
    and is non-empty.
    """
    runner = _make_runner(tmp_path, paired_sample)
    runner.output_dir.mkdir(parents = True, exist_ok = True)
    for path in (runner.output_r1, runner.output_r2, runner.output_html, runner.output_json):
        path.write_text("data")
    runner.validate_outputs()
