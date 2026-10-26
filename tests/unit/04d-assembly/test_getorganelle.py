"""
test_getorganelle.py

Unit tests for the GetOrganelleRunner class.
"""

# Imports
import json
import logging
from pathlib import Path
import pytest
from mitopipeline.api.getorganelle import (
    GetOrganelleRunner,
    build_runner_from_seed_result,
    discover_assembly_artifacts,
)
from mitopipeline.models.sample.sample import Sample
from mitopipeline.utils.seed.models import SeedSelectionResult

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

def test_build_command_does_not_create_output_directory(tmp_path: Path, paired_sample: Sample) -> None:
    """
    build_command must NOT create the output directory: GetOrganelle
    refuses to run if its -o directory already exists ("Please use a
    new output directory, or use '--continue'/'--overwrite'"), and
    creates it itself on a successful run.
    """
    runner = _make_runner(tmp_path, paired_sample)
    assert not runner.output_dir.exists()
    runner.build_command()
    assert not runner.output_dir.exists()

def test_build_command_creates_missing_parent_directories(tmp_path: Path, paired_sample: Sample) -> None:
    """GetOrganelle's own directory creation is a bare os.mkdir(), not
    recursive -- confirmed against a real failure (FileNotFoundError)
    on a nested output_dir like a rescue attempt's "assembly/attempt0"
    where "assembly/" didn't exist yet. The leaf itself still must not
    exist (see the test above), but every directory ABOVE it must."""
    runner = _make_runner(tmp_path, paired_sample)
    runner.output_dir = tmp_path / "assembly" / "attempt0"
    assert not runner.output_dir.parent.exists()
    runner.build_command()
    assert runner.output_dir.parent.is_dir()
    assert not runner.output_dir.exists()

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
    assert remaining == {"data.fasta", "graph.gfa", "assembly.log", runner.recruited_reads_qc_record_path.name}

def test_postprocess_outputs_captures_assembly_status_before_cleanup(tmp_path: Path, paired_sample: Sample) -> None:
    """
    A caller that needs the circular/fragmented/missing classification
    after run() returns must be able to read it from self.assembly_status,
    since by then postprocess_outputs() has already renamed/deleted the
    native *.path_sequence.fasta/*.selected_graph.gfa files
    discover_assembly_artifacts() globs for -- calling it again afterwards
    would always find zero candidates, regardless of the real outcome.
    """
    runner = _make_runner(tmp_path, paired_sample)
    assert runner.assembly_status is None
    (runner.output_dir / "k75.complete.path_sequence.fasta").parent.mkdir(parents = True, exist_ok = True)
    (runner.output_dir / "k75.complete.path_sequence.fasta").write_text(">seq\nACGT\n")
    (runner.output_dir / "k75.complete.selected_graph.gfa").write_text("H\tVN:Z:1.0\n")
    runner.postprocess_outputs()
    assert runner.assembly_status is not None
    assert runner.assembly_status.status == "circular"
    assert runner.assembly_status.candidate_count == 1
    # The native filenames are gone by now -- confirming a second,
    # post-hoc discover_assembly_artifacts() call would be too late.
    assert not list(runner.output_dir.glob("*.path_sequence.fasta"))

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

def test_discover_assembly_artifacts_missing_when_no_candidates(tmp_path: Path) -> None:
    """
    No native FASTA candidate at all means the assembly failed outright.
    """
    output_dir = tmp_path / "assembly"
    output_dir.mkdir()
    result = discover_assembly_artifacts(output_dir, logging.getLogger("test"))
    assert result.status == "missing"
    assert result.fasta_path is None
    assert result.gfa_path is None
    assert result.candidate_count == 0

def test_discover_assembly_artifacts_fragmented_when_no_complete_candidate(tmp_path: Path) -> None:
    """
    A candidate exists but GetOrganelle never marked one complete --
    status is "fragmented", not "circular".
    """
    output_dir = tmp_path / "assembly"
    _write_native_outputs(output_dir, prefix = "k75")
    result = discover_assembly_artifacts(output_dir, logging.getLogger("test"))
    assert result.status == "fragmented"
    assert result.fasta_path == output_dir / "k75.path_sequence.fasta"
    assert result.candidate_count == 1

def test_discover_assembly_artifacts_circular_when_complete_candidate_present(tmp_path: Path) -> None:
    """
    A ".complete." candidate is preferred and reported as "circular",
    even when a fragmented candidate also exists.
    """
    output_dir = tmp_path / "assembly"
    _write_native_outputs(output_dir, prefix = "k75")
    (output_dir / "k95.complete.path_sequence.fasta").write_text(">seq2\nGGGG\n")
    (output_dir / "k95.complete.selected_graph.gfa").write_text("H\tVN:Z:1.0\n")
    result = discover_assembly_artifacts(output_dir, logging.getLogger("test"))
    assert result.status == "circular"
    assert result.fasta_path == output_dir / "k95.complete.path_sequence.fasta"
    assert result.gfa_path == output_dir / "k95.complete.selected_graph.gfa"
    assert result.candidate_count == 2

def test_discover_assembly_artifacts_multiple_complete_candidates_picks_one_deterministically(tmp_path: Path) -> None:
    """
    Two ".complete." candidates is an unusual but real possibility
    (e.g. two GetOrganelle rounds both circularized). Still reported as
    "circular" -- ambiguity between two already-complete assemblies is
    not the same failure mode as having none -- but the choice must be
    deterministic (lexical order, same tiebreak glob() already applies
    everywhere else) rather than arbitrary, so repeat runs over the same
    directory never silently pick a different candidate.
    """
    output_dir = tmp_path / "assembly"
    _write_native_outputs(output_dir, prefix = "k75")
    (output_dir / "k85.complete.path_sequence.fasta").write_text(">seq1\nAAAA\n")
    (output_dir / "k85.complete.selected_graph.gfa").write_text("H\tVN:Z:1.0\n")
    (output_dir / "k95.complete.path_sequence.fasta").write_text(">seq2\nGGGG\n")
    (output_dir / "k95.complete.selected_graph.gfa").write_text("H\tVN:Z:1.0\n")

    first = discover_assembly_artifacts(output_dir, logging.getLogger("test"))
    second = discover_assembly_artifacts(output_dir, logging.getLogger("test"))

    assert first.status == "circular"
    assert first.candidate_count == 3
    assert first.fasta_path == second.fasta_path  # deterministic across repeat calls
    assert first.fasta_path == output_dir / "k85.complete.path_sequence.fasta"  # lexically first

def _seed_result(fasta_path: Path | None, fallback_used: bool) -> SeedSelectionResult:
    """
    Builds a minimal SeedSelectionResult for wiring tests.
    """
    return SeedSelectionResult(
        sample_id = "SAMPLE_01", scientific_name = "Danio rerio", gene = "12S_rRNA",
        accession = "ACC123" if not fallback_used else "GETORGANELLE_DEFAULT",
        source = "mitofish" if not fallback_used else "fallback_default_database",
        fasta_path = fasta_path, fallback_used = fallback_used,
    )

def test_build_runner_from_seed_result_wires_resolved_seed(tmp_path: Path, paired_sample: Sample) -> None:
    """
    A genuine (non-fallback) seed result's FASTA is wired in as the
    -s seed file.
    """
    seed_path = tmp_path / "seed.fasta"
    seed_result = _seed_result(fasta_path = seed_path, fallback_used = False)
    runner = build_runner_from_seed_result(seed_result, tmp_path, tmp_path / "assembly", paired_sample)
    assert runner.tool_options["seed_file"] == str(seed_path)
    command = runner.build_command()
    assert "-s" in command
    assert command[command.index("-s") + 1] == str(seed_path)

def test_build_runner_from_seed_result_leaves_fallback_unseeded(tmp_path: Path, paired_sample: Sample) -> None:
    """
    A fallback seed result (no seed found) produces an unseeded
    runner -- GetOrganelle falls back to its own default database.
    """
    seed_result = _seed_result(fasta_path = None, fallback_used = True)
    runner = build_runner_from_seed_result(seed_result, tmp_path, tmp_path / "assembly", paired_sample)
    assert "seed_file" not in runner.tool_options
    command = runner.build_command()
    assert "-s" not in command

def test_build_runner_from_seed_result_defaults_to_animal_mt(tmp_path: Path, paired_sample: Sample) -> None:
    """
    organelle_type defaults to "animal_mt" when not explicitly overridden.
    """
    seed_result = _seed_result(fasta_path = None, fallback_used = True)
    runner = build_runner_from_seed_result(seed_result, tmp_path, tmp_path / "assembly", paired_sample)
    assert runner.organelle_type == "animal_mt"

def test_build_runner_from_seed_result_does_not_mutate_caller_options(tmp_path: Path, paired_sample: Sample) -> None:
    """
    The caller's own tool_options dict must never be mutated -- a
    resolved seed is merged into a copy.
    """
    caller_options = {"max_rounds": 10}
    seed_result = _seed_result(fasta_path = tmp_path / "seed.fasta", fallback_used = False)
    build_runner_from_seed_result(seed_result, tmp_path, tmp_path / "assembly", paired_sample, tool_options = caller_options)
    assert caller_options == {"max_rounds": 10}

def test_init_defaults_keep_temp_files_without_mutating_caller_dict(tmp_path: Path, paired_sample: Sample) -> None:
    """
    keep_temp_files defaults to True (needed so GetOrganelle's own
    recruited-reads files survive for QC), via a COPY of tool_options --
    the caller's own dict must come back untouched.
    """
    caller_options: dict = {}
    runner = _make_runner(tmp_path, paired_sample, tool_options = caller_options)
    assert runner.tool_options["keep_temp_files"] is True
    assert caller_options == {}

def test_init_respects_explicit_keep_temp_files_false(tmp_path: Path, paired_sample: Sample) -> None:
    """
    An explicit keep_temp_files=False from the caller is not
    overridden -- setdefault only fills in when the key is absent.
    """
    runner = _make_runner(tmp_path, paired_sample, tool_options = {"keep_temp_files": False})
    assert runner.tool_options["keep_temp_files"] is False

def test_build_command_uses_the_real_keep_temp_flag(tmp_path: Path, paired_sample: Sample) -> None:
    """
    GetOrganelle's real flag is --keep-temp (confirmed against its own
    --help / argparse dest), not the --keep-temp-files this used to
    emit -- which GetOrganelle would have rejected as unrecognized.
    """
    runner = _make_runner(tmp_path, paired_sample)
    command = runner.build_command()
    assert "--keep-temp" in command
    assert "--keep-temp-files" not in command

def test_recruited_reads_paths_default_prefix(tmp_path: Path, paired_sample: Sample) -> None:
    """
    With no "prefix" tool option, the recruited-reads files are named
    exactly as GetOrganelle itself names them with its own default
    (empty) --prefix.
    """
    runner = _make_runner(tmp_path, paired_sample)
    r1, r2 = runner._recruited_reads_paths()
    assert r1 == runner.output_dir / "extended_1_paired.fq"
    assert r2 == runner.output_dir / "extended_2_paired.fq"

def test_recruited_reads_paths_honors_custom_prefix(tmp_path: Path, paired_sample: Sample) -> None:
    """
    A configured "prefix" tool option is reflected in the expected
    recruited-reads filenames, matching GetOrganelle's own --prefix
    behavior.
    """
    runner = _make_runner(tmp_path, paired_sample, tool_options = {"prefix": "sample01."})
    r1, r2 = runner._recruited_reads_paths()
    assert r1.name == "sample01.extended_1_paired.fq"
    assert r2.name == "sample01.extended_2_paired.fq"

def test_postprocess_outputs_skips_qc_and_cleans_up_when_recruited_reads_absent(tmp_path: Path, paired_sample: Sample) -> None:
    """
    No real recruited-reads files (e.g. an older GetOrganelle, or
    keep_temp_files didn't take effect) must not fail the whole
    assembly over a missing supplementary QC artifact -- cleanup still
    runs normally, no FastQC directory gets created, but a record of
    the absence is still written explicitly (#413), not just logged.
    """
    runner = _make_runner(tmp_path, paired_sample)
    _write_native_outputs(runner.output_dir, prefix = "k75")
    runner.postprocess_outputs()
    remaining = {p.name for p in runner.output_dir.iterdir()}
    assert remaining == {"data.fasta", "graph.gfa", "assembly.log", runner.recruited_reads_qc_record_path.name}
    assert not runner.recruited_reads_fastqc_dir.exists()
    record = json.loads(runner.recruited_reads_qc_record_path.read_text())
    assert record["qc_ran"] is False
    assert record["reason"] == "recruited-reads files not found"

def test_postprocess_outputs_retains_recruited_reads_fastqc_dir(tmp_path: Path, paired_sample: Sample) -> None:
    """
    If FastQC output already exists in the recruited-reads QC
    directory (simulating a completed QC run), the general cleanup
    pass must not remove it, same as the standardized FASTA/GFA/log.
    """
    runner = _make_runner(tmp_path, paired_sample)
    _write_native_outputs(runner.output_dir, prefix = "k75")
    runner.recruited_reads_fastqc_dir.mkdir()
    (runner.recruited_reads_fastqc_dir / "recruited_R1_fastqc.html").write_text("<html></html>")
    runner.postprocess_outputs()
    remaining = {p.name for p in runner.output_dir.iterdir()}
    assert "recruited_reads_fastqc" in remaining
    assert (runner.recruited_reads_fastqc_dir / "recruited_R1_fastqc.html").is_file()

# ---------------------------------------------------------------------------
# #413: recruited-reads QC artifact contract -- nonempty/pairing validation,
# and the record written for every outcome (success, empty, malformed,
# pairing mismatch, and FastQC tool failure).
# ---------------------------------------------------------------------------

def _fastq_text(n_records: int) -> str:
    """Builds a minimal, well-formed FASTQ body with n_records records."""
    return "".join(f"@read{i}\nACGT\n+\nIIII\n" for i in range(n_records))

def _fake_fastqc_result(success: bool = True, return_code: int = 0):
    from datetime import datetime
    from mitopipeline.models.command_result.command_result import CommandResult
    now = datetime.now()
    return CommandResult(
        command = [], return_code = return_code, stdout = "", stderr = "", runtime_seconds = 0.5,
        success = success, tool_name = "fastqc", started_at = now, ended_at = now,
        peak_memory_mb = 0.0, user_cpu_seconds = 0.0, system_cpu_seconds = 0.0,
    )

def test_recruited_reads_qc_record_written_on_success(tmp_path: Path, paired_sample: Sample, monkeypatch: pytest.MonkeyPatch) -> None:
    from mitopipeline.api.fastqc import FastQCRunner
    runner = _make_runner(tmp_path, paired_sample)
    runner.output_dir.mkdir(parents = True)
    r1, r2 = runner._recruited_reads_paths()
    r1.write_text(_fastq_text(3))
    r2.write_text(_fastq_text(3))

    monkeypatch.setattr(FastQCRunner, "run", lambda self: _fake_fastqc_result(success = True))
    runner._run_recruited_reads_qc()

    record = json.loads(runner.recruited_reads_qc_record_path.read_text())
    assert record["qc_ran"] is True
    assert record["r1_read_count"] == 3
    assert record["r2_read_count"] == 3
    assert record["paired"] is True
    assert record["fastqc_dir"] == str(runner.recruited_reads_fastqc_dir)
    assert record["reason"] == "qc completed"

def test_recruited_reads_qc_record_pairing_mismatch_still_runs_qc(tmp_path: Path, paired_sample: Sample, monkeypatch: pytest.MonkeyPatch) -> None:
    from mitopipeline.api.fastqc import FastQCRunner
    runner = _make_runner(tmp_path, paired_sample)
    runner.output_dir.mkdir(parents = True)
    r1, r2 = runner._recruited_reads_paths()
    r1.write_text(_fastq_text(5))
    r2.write_text(_fastq_text(3))

    monkeypatch.setattr(FastQCRunner, "run", lambda self: _fake_fastqc_result(success = True))
    runner._run_recruited_reads_qc()

    record = json.loads(runner.recruited_reads_qc_record_path.read_text())
    assert record["r1_read_count"] == 5
    assert record["r2_read_count"] == 3
    assert record["paired"] is False
    assert record["qc_ran"] is True   # pairing mismatch does not block QC itself

def test_recruited_reads_qc_skipped_when_files_empty(tmp_path: Path, paired_sample: Sample, monkeypatch: pytest.MonkeyPatch) -> None:
    from mitopipeline.api.fastqc import FastQCRunner
    runner = _make_runner(tmp_path, paired_sample)
    runner.output_dir.mkdir(parents = True)
    r1, r2 = runner._recruited_reads_paths()
    r1.touch(); r2.touch()

    def _must_not_run(self):
        raise AssertionError("FastQC must not run on empty recruited-reads files")
    monkeypatch.setattr(FastQCRunner, "run", _must_not_run)

    runner._run_recruited_reads_qc()
    record = json.loads(runner.recruited_reads_qc_record_path.read_text())
    assert record["qc_ran"] is False
    assert record["paired"] is None
    assert record["reason"] == "recruited-reads file(s) empty"

def test_recruited_reads_qc_malformed_fastq_runs_qc_with_unknown_pairing(tmp_path: Path, paired_sample: Sample, monkeypatch: pytest.MonkeyPatch) -> None:
    from mitopipeline.api.fastqc import FastQCRunner
    runner = _make_runner(tmp_path, paired_sample)
    runner.output_dir.mkdir(parents = True)
    r1, r2 = runner._recruited_reads_paths()
    r1.write_text("@read0\nACGT\n+\n")   # 3 lines -- not a multiple of 4
    r2.write_text(_fastq_text(2))

    monkeypatch.setattr(FastQCRunner, "run", lambda self: _fake_fastqc_result(success = True))
    runner._run_recruited_reads_qc()

    record = json.loads(runner.recruited_reads_qc_record_path.read_text())
    assert record["r1_read_count"] is None
    assert record["r2_read_count"] == 2
    assert record["paired"] is None
    assert record["qc_ran"] is True   # malformed counting does not block QC itself

def test_recruited_reads_qc_tool_failure_is_recorded_not_raised(tmp_path: Path, paired_sample: Sample, monkeypatch: pytest.MonkeyPatch) -> None:
    from mitopipeline.api.fastqc import FastQCRunner
    runner = _make_runner(tmp_path, paired_sample)
    runner.output_dir.mkdir(parents = True)
    r1, r2 = runner._recruited_reads_paths()
    r1.write_text(_fastq_text(1)); r2.write_text(_fastq_text(1))

    monkeypatch.setattr(FastQCRunner, "run", lambda self: _fake_fastqc_result(success = False, return_code = 1))
    runner._run_recruited_reads_qc()   # must not raise

    record = json.loads(runner.recruited_reads_qc_record_path.read_text())
    assert record["qc_ran"] is False
    assert record["reason"] == "fastqc exited 1"
    assert record["fastqc_dir"] is None

def test_recruited_reads_qc_tool_systemexit_is_recorded_not_raised(tmp_path: Path, paired_sample: Sample, monkeypatch: pytest.MonkeyPatch) -> None:
    from mitopipeline.api.fastqc import FastQCRunner
    runner = _make_runner(tmp_path, paired_sample)
    runner.output_dir.mkdir(parents = True)
    r1, r2 = runner._recruited_reads_paths()
    r1.write_text(_fastq_text(1)); r2.write_text(_fastq_text(1))

    def _raise(self):
        raise SystemExit(1)
    monkeypatch.setattr(FastQCRunner, "run", _raise)
    runner._run_recruited_reads_qc()   # must not raise

    record = json.loads(runner.recruited_reads_qc_record_path.read_text())
    assert record["qc_ran"] is False
    assert record["reason"] == "fastqc failed validation"

def test_recruited_reads_qc_record_is_attempt_specific(tmp_path: Path, paired_sample: Sample) -> None:
    """
    Two runners for different rescue attempts (different output_dir)
    must never share a record path -- each attempt's output_dir is
    already distinct, so the record naturally never collides.
    """
    attempt1 = _make_runner(tmp_path, paired_sample, output_dir = tmp_path / "rescue_attempt_1")
    attempt2 = _make_runner(tmp_path, paired_sample, output_dir = tmp_path / "rescue_attempt_2")
    assert attempt1.recruited_reads_qc_record_path != attempt2.recruited_reads_qc_record_path
    assert attempt1.recruited_reads_qc_record_path.parent == attempt1.output_dir
    assert attempt2.recruited_reads_qc_record_path.parent == attempt2.output_dir
