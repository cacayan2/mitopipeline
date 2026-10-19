"""
test_blast.py

Unit tests for the BlastRunner class.
"""

# Imports
import logging
from pathlib import Path
from unittest.mock import patch
import pytest
from mitopipeline.api.blast import BlastRunner, DEFAULT_OUTPUT_FORMAT, BLAST_OUTPUT_COLUMNS
from mitopipeline.models.sample.sample import Sample

@pytest.fixture
def sample(tmp_path: Path) -> Sample:
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

def _make_runner(tmp_path: Path, sample: Sample, **overrides) -> BlastRunner:
    """
    Builds a BlastRunner with sensible defaults, overridable per test.

    Args:
        tmp_path (Path): The pytest-provided temporary directory.
        sample (Sample): The sample to BLAST.
        **overrides: Keyword arguments overriding the defaults.

    Returns:
        BlastRunner: The constructed runner.
    """
    query_fasta = tmp_path / "query.fasta"
    query_fasta.touch()
    defaults = dict(
        working_dir = tmp_path,
        query_fasta = query_fasta,
        database = "test_db",
        output_file = tmp_path / "blast_out.tsv",
        sample = sample,
        logger = logging.getLogger("test"),
    )
    defaults.update(overrides)
    return BlastRunner(**defaults)

def test_default_output_format_requests_all_parser_columns() -> None:
    """
    The default -outfmt must explicitly list every column
    blast_hits.parse_blast_tsv expects -- bare "6" alone only yields
    the first 12 of BLAST's default columns and silently drops
    qcovs/staxids/sscinames/stitle, breaking the parser on every row.
    """
    for column in BLAST_OUTPUT_COLUMNS:
        assert column in DEFAULT_OUTPUT_FORMAT

def test_build_command_local_mode_includes_threads(tmp_path: Path, sample: Sample) -> None:
    """
    A local-mode run should pass -num_threads and omit -remote.
    """
    runner = _make_runner(tmp_path, sample, threads = 8)
    command = runner.build_command()
    assert "-remote" not in command
    assert "-num_threads" in command
    assert command[command.index("-num_threads") + 1] == "8"

def test_build_command_remote_mode_omits_threads(tmp_path: Path, sample: Sample) -> None:
    """
    A remote-mode run should pass -remote and omit -num_threads.
    """
    runner = _make_runner(tmp_path, sample, mode = "remote")
    command = runner.build_command()
    assert "-remote" in command
    assert "-num_threads" not in command

def test_build_command_omits_optional_filters_by_default(tmp_path: Path, sample: Sample) -> None:
    """
    perc_identity, query_coverage, and word_size should not appear in
    the command unless explicitly configured.
    """
    runner = _make_runner(tmp_path, sample)
    command = runner.build_command()
    assert "-perc_identity" not in command
    assert "-qcov_hsp_perc" not in command
    assert "-word_size" not in command

def test_build_command_includes_configured_optional_filters(tmp_path: Path, sample: Sample) -> None:
    """
    Explicitly configured optional filters should appear with their
    values.
    """
    runner = _make_runner(tmp_path, sample, perc_identity = 90.0, query_coverage = 80.0, word_size = 11)
    command = runner.build_command()
    assert command[command.index("-perc_identity") + 1] == "90.0"
    assert command[command.index("-qcov_hsp_perc") + 1] == "80.0"
    assert command[command.index("-word_size") + 1] == "11"

def test_build_command_creates_output_directory(tmp_path: Path, sample: Sample) -> None:
    """
    build_command should create the output file's parent directory if
    it does not exist.
    """
    runner = _make_runner(tmp_path, sample, output_file = tmp_path / "nested" / "out.tsv")
    assert not runner.output_file.parent.exists()
    runner.build_command()
    assert runner.output_file.parent.is_dir()

def test_validate_inputs_exits_on_missing_query(tmp_path: Path, sample: Sample) -> None:
    """
    A missing query FASTA should exit.
    """
    runner = _make_runner(tmp_path, sample, query_fasta = tmp_path / "missing.fasta")
    with pytest.raises(SystemExit):
        runner.validate_inputs()

def test_validate_inputs_exits_on_unsupported_mode(tmp_path: Path, sample: Sample) -> None:
    """
    An unsupported mode should exit before any command is built.
    """
    runner = _make_runner(tmp_path, sample, mode = "not_a_real_mode")
    with pytest.raises(SystemExit):
        runner.validate_inputs()

def test_validate_inputs_exits_on_unsupported_task(tmp_path: Path, sample: Sample) -> None:
    """
    An unsupported BLAST task should exit before any command is built.
    """
    runner = _make_runner(tmp_path, sample, task = "not_a_real_task")
    with pytest.raises(SystemExit):
        runner.validate_inputs()

def test_validate_inputs_exits_on_invalid_local_database(tmp_path: Path, sample: Sample) -> None:
    """
    A local-mode run whose database doesn't resolve via blastdbcmd
    should exit.
    """
    runner = _make_runner(tmp_path, sample)
    with patch("subprocess.run") as mock_run:
        mock_run.return_value.returncode = 1
        mock_run.return_value.stderr = "no such database"
        with pytest.raises(SystemExit):
            runner.validate_inputs()

def test_validate_inputs_skips_database_check_in_remote_mode(tmp_path: Path, sample: Sample) -> None:
    """
    A remote-mode run never calls blastdbcmd -- there's no local
    database to check.
    """
    runner = _make_runner(tmp_path, sample, mode = "remote")
    with patch("subprocess.run") as mock_run:
        runner.validate_inputs()
        mock_run.assert_not_called()

def test_validate_inputs_exits_on_too_many_threads(tmp_path: Path, sample: Sample) -> None:
    """
    A thread count exceeding the machine's CPU count should exit.
    """
    runner = _make_runner(tmp_path, sample, threads = 10_000)
    with patch("subprocess.run") as mock_run:
        mock_run.return_value.returncode = 0
        with pytest.raises(SystemExit):
            runner.validate_inputs()

def test_validate_outputs_exits_on_missing_output(tmp_path: Path, sample: Sample) -> None:
    """
    A missing output file after a run should exit.
    """
    runner = _make_runner(tmp_path, sample)
    with pytest.raises(SystemExit):
        runner.validate_outputs()

def test_validate_outputs_accepts_empty_output(tmp_path: Path, sample: Sample) -> None:
    """
    An empty output file is a valid "zero hits" result, not a
    failure -- only a genuinely missing file should exit.
    """
    runner = _make_runner(tmp_path, sample)
    runner.output_file.parent.mkdir(parents = True, exist_ok = True)
    runner.output_file.touch()
    runner.validate_outputs()
