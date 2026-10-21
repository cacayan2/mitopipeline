"""
test_getorganelle.py

Unit tests for the GetOrganelleRunner class.
"""

# Imports
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
    runs normally and no FastQC directory gets created.
    """
    runner = _make_runner(tmp_path, paired_sample)
    _write_native_outputs(runner.output_dir, prefix = "k75")
    runner.postprocess_outputs()
    remaining = {p.name for p in runner.output_dir.iterdir()}
    assert remaining == {"data.fasta", "graph.gfa", "assembly.log"}
    assert not runner.recruited_reads_fastqc_dir.exists()

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
