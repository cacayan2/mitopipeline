"""
test_getorganelle.py

Unit tests for the GetOrganelleRunner class.
"""

# Imports
import logging
from pathlib import Path
import pytest
from mitopipeline.api.getorganelle import GetOrganelleRunner
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
    r1 = tmp_path / "trimmed_R1.fastq.gz"
    r2 = tmp_path / "trimmed_R2.fastq.gz"
    r1.touch()
    r2.touch()
    return Sample(sample_id = "SAMPLE_01", r1 = r1, r2 = r2)

def _make_runner(tmp_path: Path, sample: Sample, **overrides) -> GetOrganelleRunner:
    """
    Builds a GetOrganelleRunner with sensible defaults, overridable per test.

    Args:
        tmp_path (Path): The pytest-provided temporary directory.
        sample (Sample): The sample to assemble.
        **overrides: Keyword arguments overriding the defaults.

    Returns:
        GetOrganelleRunner: The constructed runner.
    """
    defaults = dict(
        working_dir = tmp_path,
        output_dir = tmp_path / "assembly",
        sample = sample,
        organelle_type = "animal_mt",
        logger = logging.getLogger("test"),
    )
    defaults.update(overrides)
    return GetOrganelleRunner(**defaults)

def test_build_command_without_seed_omits_seed_flag(tmp_path: Path, paired_sample: Sample) -> None:
    """
    A run with no "seed_file" tool option is a valid, intentional
    no-seed run -- it must fall back to GetOrganelle's own default
    database for the organelle type, not append a bare/empty -s flag.
    """
    runner = _make_runner(tmp_path, paired_sample)
    command = runner.build_command()
    assert "-s" not in command
    assert "-F" in command and command[command.index("-F") + 1] == "animal_mt"
    assert "-1" in command and command[command.index("-1") + 1] == str(paired_sample.r1)

def test_build_command_with_seed_includes_seed_flag(tmp_path: Path, paired_sample: Sample) -> None:
    """
    A "seed_file" tool option should be passed through as GetOrganelle's
    -s flag.
    """
    seed_path = tmp_path / "seed.fasta"
    runner = _make_runner(tmp_path, paired_sample, tool_options = {"seed_file": str(seed_path)})
    command = runner.build_command()
    assert "-s" in command
    assert command[command.index("-s") + 1] == str(seed_path)

def test_build_command_creates_output_directory(tmp_path: Path, paired_sample: Sample) -> None:
    """
    build_command should create the output directory if it does not exist.
    """
    runner = _make_runner(tmp_path, paired_sample)
    assert not runner.output_dir.exists()
    runner.build_command()
    assert runner.output_dir.is_dir()

def test_validate_inputs_exits_on_unsupported_organelle_type(tmp_path: Path, paired_sample: Sample) -> None:
    """
    An organelle type GetOrganelle itself doesn't support should exit
    before any command is built.
    """
    runner = _make_runner(tmp_path, paired_sample, organelle_type = "not_a_real_type")
    with pytest.raises(SystemExit):
        runner.validate_inputs()

def test_validate_inputs_exits_on_missing_read_file(tmp_path: Path) -> None:
    """
    A missing R1/R2 file should exit.
    """
    sample = Sample(sample_id = "SAMPLE_02", r1 = tmp_path / "missing_R1.fastq.gz", r2 = tmp_path / "missing_R2.fastq.gz")
    runner = _make_runner(tmp_path, sample)
    with pytest.raises(SystemExit):
        runner.validate_inputs()

def test_validate_inputs_passes_for_valid_sample(tmp_path: Path, paired_sample: Sample) -> None:
    """
    A valid organelle type, real input files, and a valid thread count
    should pass validation without exiting.
    """
    runner = _make_runner(tmp_path, paired_sample, threads = 1)
    runner.validate_inputs()

def _write_native_outputs(output_dir: Path, prefix: str = "run") -> None:
    """
    Simulates GetOrganelle's own native, versioned output filenames plus
    one piece of scratch state it normally leaves behind.

    Args:
        output_dir (Path): The GetOrganelle output directory.
        prefix (str): A filename prefix, standing in for GetOrganelle's own naming.
    """
    output_dir.mkdir(parents = True, exist_ok = True)
    (output_dir / f"{prefix}.path_sequence.fasta").write_text(">seq\nACGT\n")
    (output_dir / f"{prefix}.selected_graph.gfa").write_text("H\tVN:Z:1.0\n")
    (output_dir / "assembly.log").write_text("log")
    (output_dir / "scratch_round_1").mkdir()
    (output_dir / "scratch_round_1" / "temp.bin").write_text("x")

def test_postprocess_outputs_standardizes_and_cleans_up(tmp_path: Path, paired_sample: Sample) -> None:
    """
    postprocess_outputs should copy GetOrganelle's versioned native
    outputs to the fixed data.fasta/graph.gfa pair, and remove every
    other scratch file/directory GetOrganelle left behind.
    """
    runner = _make_runner(tmp_path, paired_sample)
    _write_native_outputs(runner.output_dir, prefix = "k75")
    runner.postprocess_outputs()
    assert runner.standardized_fasta.read_text() == ">seq\nACGT\n"
    assert runner.standardized_gfa.read_text() == "H\tVN:Z:1.0\n"
    remaining = {p.name for p in runner.output_dir.iterdir()}
    assert remaining == {"data.fasta", "graph.gfa", "assembly.log"}

def test_postprocess_outputs_exits_on_ambiguous_native_fasta(tmp_path: Path, paired_sample: Sample) -> None:
    """
    More than one candidate native FASTA means GetOrganelle produced
    something unexpected; this must exit rather than silently pick one.
    """
    runner = _make_runner(tmp_path, paired_sample)
    _write_native_outputs(runner.output_dir, prefix = "k75")
    (runner.output_dir / "k95.path_sequence.fasta").write_text(">seq2\nGGGG\n")
    with pytest.raises(SystemExit):
        runner.postprocess_outputs()

def test_validate_outputs_exits_on_missing_output(tmp_path: Path, paired_sample: Sample) -> None:
    """
    validate_outputs should exit if the standardized FASTA/GFA are missing.
    """
    runner = _make_runner(tmp_path, paired_sample)
    runner.output_dir.mkdir(parents = True, exist_ok = True)
    with pytest.raises(SystemExit):
        runner.validate_outputs()

def test_validate_outputs_passes_when_all_outputs_present(tmp_path: Path, paired_sample: Sample) -> None:
    """
    validate_outputs should pass once both standardized outputs exist
    and are non-empty.
    """
    runner = _make_runner(tmp_path, paired_sample)
    runner.output_dir.mkdir(parents = True, exist_ok = True)
    runner.standardized_fasta.write_text(">seq\nACGT\n")
    runner.standardized_gfa.write_text("H\tVN:Z:1.0\n")
    runner.validate_outputs()
