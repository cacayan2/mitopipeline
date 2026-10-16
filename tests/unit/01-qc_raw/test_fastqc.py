"""
test_fastqc.py

Unit tests for the FastQCRunner class.
"""

# Imports
import logging
from pathlib import Path
import pytest
from mitopipeline.api.fastqc import FastQCRunner
from mitopipeline.models.sample.sample import Sample

@pytest.fixture
def paired_sample(tmp_path: Path) -> Sample:
    """
    Creates a Sample backed by two real (empty) FASTQ files.

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

def _make_runner(tmp_path: Path, sample: Sample, **overrides) -> FastQCRunner:
    """
    Builds a FastQCRunner with sensible defaults, overridable per test.

    Args:
        tmp_path (Path): The pytest-provided temporary directory.
        sample (Sample): The sample to run FastQC on.
        **overrides: Keyword arguments overriding the defaults.

    Returns:
        FastQCRunner: The constructed runner.
    """
    defaults = dict(
        working_dir = tmp_path,
        output_dir = tmp_path / "qc",
        sample = sample,
        r1_output_stem = "R1",
        r2_output_stem = "R2",
        logger = logging.getLogger("test"),
    )
    defaults.update(overrides)
    return FastQCRunner(**defaults)

def test_init_sets_working_dir_on_base_tool(tmp_path: Path, paired_sample: Sample) -> None:
    """
    A regression check for a typo'd "workind_dir" keyword that would have
    made BaseTool.__init__ raise, and left working_dir unset entirely.
    """
    runner = _make_runner(tmp_path, paired_sample)
    assert runner.working_dir == tmp_path

def test_build_command_includes_required_flags(tmp_path: Path, paired_sample: Sample) -> None:
    """
    build_command should wire up fastqc's required flags and inputs.
    """
    runner = _make_runner(tmp_path, paired_sample, threads = 4)
    command = runner.build_command()
    assert command[0] == "fastqc"
    assert "--threads" in command
    assert command[command.index("--threads") + 1] == "4"
    assert str(paired_sample.r1) in command
    assert str(paired_sample.r2) in command
    assert "-o" in command

def test_build_command_creates_output_directory(tmp_path: Path, paired_sample: Sample) -> None:
    """
    build_command should create the output directory if it does not exist.
    """
    runner = _make_runner(tmp_path, paired_sample)
    assert not runner.output_dir.exists()
    runner.build_command()
    assert runner.output_dir.is_dir()

def test_validate_inputs_exits_on_missing_read_file(tmp_path: Path) -> None:
    """
    A missing R1/R2 file should exit.
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

def test_validate_inputs_exits_on_too_many_threads(tmp_path: Path, paired_sample: Sample) -> None:
    """
    A thread count exceeding the machine's CPU count should exit --
    regression check for a missing sys.exit after the error was logged.
    """
    import os
    runner = _make_runner(tmp_path, paired_sample, threads = (os.cpu_count() or 1) + 1000)
    with pytest.raises(SystemExit):
        runner.validate_inputs()

def test_postprocess_outputs_renames_native_outputs(tmp_path: Path, paired_sample: Sample) -> None:
    """
    postprocess_outputs should move FastQC's native per-input-filename
    outputs to the configured R1/R2 stems.
    """
    runner = _make_runner(tmp_path, paired_sample)
    runner.output_dir.mkdir(parents = True, exist_ok = True)
    (runner.output_dir / "raw_R1_fastqc.html").write_text("a")
    (runner.output_dir / "raw_R1_fastqc.zip").write_text("b")
    (runner.output_dir / "raw_R2_fastqc.html").write_text("c")
    (runner.output_dir / "raw_R2_fastqc.zip").write_text("d")
    runner.postprocess_outputs()
    assert (runner.output_dir / "R1_fastqc.html").read_text() == "a"
    assert (runner.output_dir / "R2_fastqc.zip").read_text() == "d"

def test_postprocess_outputs_exits_on_missing_native_output(tmp_path: Path, paired_sample: Sample) -> None:
    """
    postprocess_outputs should exit if a native FastQC output never
    appeared, rather than crash later trying to move a nonexistent file.
    """
    runner = _make_runner(tmp_path, paired_sample)
    runner.output_dir.mkdir(parents = True, exist_ok = True)
    with pytest.raises(SystemExit):
        runner.postprocess_outputs()

def test_validate_outputs_exits_on_missing_output(tmp_path: Path, paired_sample: Sample) -> None:
    """
    validate_outputs should exit if an expected (renamed) output is missing.
    """
    runner = _make_runner(tmp_path, paired_sample)
    runner.output_dir.mkdir(parents = True, exist_ok = True)
    with pytest.raises(SystemExit):
        runner.validate_outputs()

def test_validate_outputs_passes_when_all_outputs_present(tmp_path: Path, paired_sample: Sample) -> None:
    """
    validate_outputs should pass once every expected renamed output exists
    and is non-empty.
    """
    runner = _make_runner(tmp_path, paired_sample)
    runner.output_dir.mkdir(parents = True, exist_ok = True)
    for stem in ("R1", "R2"):
        (runner.output_dir / f"{stem}_fastqc.html").write_text("x")
        (runner.output_dir / f"{stem}_fastqc.zip").write_text("y")
    runner.validate_outputs()

@pytest.mark.parametrize("name,expected", [
    ("sample_R1.fastq.gz", "sample_R1"),
    ("sample_R1.fq.gz", "sample_R1"),
    ("sample_R1.fastq", "sample_R1"),
    ("sample_R1.fq", "sample_R1"),
    ("sample_R1.bam", "sample_R1"),
])
def test_strip_fastq_suffix(name: str, expected: str) -> None:
    """
    _strip_fastq_suffix should remove every recognized FASTQ suffix, and
    fall back to Path.stem (stripping just the final extension) for an
    unrecognized one.
    """
    assert FastQCRunner._strip_fastq_suffix(Path(name)) == expected
