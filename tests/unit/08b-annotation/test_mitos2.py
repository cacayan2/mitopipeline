"""
test_mitos2.py

Unit tests for the Mitos2Runner class (#409): builds the real
runmitos.py command for FASTA input, genetic code, taxonomy
(refdir/refseqver), output paths, and nonzero-exit behavior.
"""

# Imports
import logging
from pathlib import Path
import pytest
from mitopipeline.api.mitos2 import Mitos2Runner

def _make_runner(tmp_path: Path, **overrides) -> Mitos2Runner:
    input_fasta = tmp_path / "SAMPLE_01.fasta"
    input_fasta.write_text(">SAMPLE_01\nACGT\n")
    refdir = tmp_path / "resources" / "mitos"
    (refdir / "refseq89m").mkdir(parents = True)
    defaults = dict(
        working_dir = tmp_path,
        output_dir = tmp_path / "annotation",
        sample_id = "SAMPLE_01",
        input_fasta = input_fasta,
        genetic_code = 2,
        refdir = refdir,
        refseqver = "refseq89m",
        logger = logging.getLogger("test"),
    )
    defaults.update(overrides)
    return Mitos2Runner(**defaults)

def test_build_command_includes_fasta_code_and_taxonomy(tmp_path: Path) -> None:
    runner = _make_runner(tmp_path)
    command = runner.build_command()
    assert command[0] == "runmitos.py"
    assert "-i" in command and command[command.index("-i") + 1] == str(runner.input_fasta)
    assert "-c" in command and command[command.index("-c") + 1] == "2"
    assert "-R" in command and command[command.index("-R") + 1] == str(runner.refdir)
    assert "-r" in command and command[command.index("-r") + 1] == "refseq89m"
    assert "-o" in command and command[command.index("-o") + 1] == str(runner.output_dir)

def test_build_command_creates_the_output_directory(tmp_path: Path) -> None:
    # Unlike GetOrganelle, runmitos.py requires -o to already exist --
    # it errors out ("no such directory") rather than creating it.
    runner = _make_runner(tmp_path)
    assert not runner.output_dir.exists()
    runner.build_command()
    assert runner.output_dir.is_dir()

def test_build_command_clears_a_preexisting_output_directory(tmp_path: Path) -> None:
    # A stale scratch file from a prior run against this same sample
    # (e.g. mitfi-global/*.cmout) genuinely crashes the next run --
    # confirmed directly against the real tool. build_command() must
    # wipe the directory before recreating it, not just mkdir over it.
    runner = _make_runner(tmp_path)
    runner.output_dir.mkdir(parents = True)
    (runner.output_dir / "mitfi-global").mkdir()
    (runner.output_dir / "mitfi-global" / "stale.cmout").write_text("stale")
    (runner.output_dir / "sequence.fas-1234").write_text("stale")
    (runner.output_dir / "SAMPLE_01.annotation.gff").write_text("stale prior result")

    runner.build_command()

    assert runner.output_dir.is_dir()
    assert list(runner.output_dir.iterdir()) == []

def test_default_treats_sequence_as_circular(tmp_path: Path) -> None:
    # runmitos.py's own default is circular=True; --linear is the
    # flag that flips it to False. Never passing --linear means
    # accepting that default.
    runner = _make_runner(tmp_path)
    command = runner.build_command()
    assert "--linear" not in command

def test_non_circular_passes_linear_flag(tmp_path: Path) -> None:
    runner = _make_runner(tmp_path, circular = False)
    command = runner.build_command()
    assert "--linear" in command

def test_boolean_tool_options_map_to_real_flags(tmp_path: Path) -> None:
    runner = _make_runner(tmp_path, tool_options = {
        "noplots": True, "best": True, "ncbicode": True, "zip_output": True,
    })
    command = runner.build_command()
    assert "--noplots" in command
    assert "--best" in command
    assert "--ncbicode" in command
    assert "--zip" in command

def test_false_boolean_tool_options_are_omitted(tmp_path: Path) -> None:
    runner = _make_runner(tmp_path, tool_options = {"noplots": False})
    command = runner.build_command()
    assert "--noplots" not in command

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

def test_validate_inputs_rejects_invalid_genetic_code(tmp_path: Path) -> None:
    runner = _make_runner(tmp_path, genetic_code = 0)
    with pytest.raises(SystemExit):
        runner.validate_inputs()

def test_validate_inputs_rejects_missing_taxonomy_reference_dir(tmp_path: Path) -> None:
    # The real tool checks os.path.join(refdir, refseqver) is a directory
    # (confirmed against the installed runmitos.py source) -- a refseqver
    # that doesn't exist under refdir must fail here, not at runtime.
    runner = _make_runner(tmp_path, refseqver = "does_not_exist")
    with pytest.raises(SystemExit):
        runner.validate_inputs()

def test_validate_inputs_rejects_empty_refseqver(tmp_path: Path) -> None:
    runner = _make_runner(tmp_path, refseqver = "")
    with pytest.raises(SystemExit):
        runner.validate_inputs()

def test_validate_inputs_passes_for_valid_configuration(tmp_path: Path) -> None:
    runner = _make_runner(tmp_path)
    runner.validate_inputs()   # does not raise

def test_validate_outputs_rejects_missing_gff(tmp_path: Path) -> None:
    runner = _make_runner(tmp_path)
    runner.output_dir.mkdir(parents = True)
    with pytest.raises(SystemExit):
        runner.validate_outputs()

def test_validate_outputs_rejects_empty_gff(tmp_path: Path) -> None:
    runner = _make_runner(tmp_path)
    runner.output_dir.mkdir(parents = True)
    runner.result_gff.touch()
    runner.result_fas.write_text(">SAMPLE_01\nACGT\n")
    with pytest.raises(SystemExit):
        runner.validate_outputs()

def test_validate_outputs_rejects_missing_fas_even_with_real_gff(tmp_path: Path) -> None:
    runner = _make_runner(tmp_path)
    runner.output_dir.mkdir(parents = True)
    runner.result_gff.write_text("##gff-version 3\n")
    with pytest.raises(SystemExit):
        runner.validate_outputs()

def test_validate_outputs_passes_with_real_content(tmp_path: Path) -> None:
    runner = _make_runner(tmp_path)
    runner.output_dir.mkdir(parents = True)
    runner.result_gff.write_text("##gff-version 3\nSAMPLE_01\tMITOS\tgene\t1\t10\t.\t+\t.\tID=nad1\n")
    runner.result_fas.write_text(">nad1\nACGT\n")
    runner.validate_outputs()   # does not raise

def test_result_paths_carry_the_sample_id_prefix(tmp_path: Path) -> None:
    runner = _make_runner(tmp_path)
    assert runner.result_gff.name == "SAMPLE_01.annotation.gff"
    assert runner.result_fas.name == "SAMPLE_01.annotation.fas"
    assert runner.result_faa.name == "SAMPLE_01.annotation.faa"
    assert runner.result_bed.name == "SAMPLE_01.annotation.bed"
    assert runner.result_table.name == "SAMPLE_01.annotation.seq"
    assert runner.result_geneorder.name == "SAMPLE_01.annotation.geneorder"

def test_postprocess_outputs_renames_every_real_native_file(tmp_path: Path) -> None:
    """runmitos.py itself offers no output-prefix option -- it always
    writes fixed names like "result.gff" -- so postprocess_outputs()
    must rename whatever it actually produced into the sample-prefixed
    paths every other consumer reads via the runner's own attributes."""
    runner = _make_runner(tmp_path)
    runner.output_dir.mkdir(parents = True)
    for native_name in ("result.gff", "result.fas", "result.faa", "result.bed",
                        "result.seq", "result.geneorder", "ignored.mitos", "stst.dat"):
        (runner.output_dir / native_name).write_text("x")

    runner.postprocess_outputs()

    assert runner.result_gff.is_file()
    assert runner.result_fas.is_file()
    assert not (runner.output_dir / "result.gff").exists()
    assert not (runner.output_dir / "ignored.mitos").exists()
    assert (runner.output_dir / "SAMPLE_01.annotation.ignored.mitos").is_file()
    assert (runner.output_dir / "SAMPLE_01.annotation.stst.dat").is_file()

def test_postprocess_outputs_skips_native_files_that_were_never_produced(tmp_path: Path) -> None:
    """Optional MITOS2 outputs (e.g. result.png, only with certain
    options) must not crash postprocess_outputs() when absent."""
    runner = _make_runner(tmp_path)
    runner.output_dir.mkdir(parents = True)
    runner.postprocess_outputs()   # does not raise with nothing present
    assert not runner.result_gff.exists()

def test_postprocess_outputs_is_idempotent_on_an_already_renamed_directory(tmp_path: Path) -> None:
    runner = _make_runner(tmp_path)
    runner.output_dir.mkdir(parents = True)
    runner.result_gff.write_text("already renamed")
    runner.postprocess_outputs()   # does not raise re-running against renamed files
    assert runner.result_gff.read_text() == "already renamed"
