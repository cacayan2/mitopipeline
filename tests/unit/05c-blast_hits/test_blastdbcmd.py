"""
test_blastdbcmd.py

Unit tests for the BlastDbCmdRunner class (#401): retrieves one
accession's FASTA from a local BLAST database.
"""

# Imports
import logging
from pathlib import Path
import pytest
from mitopipeline.api.blastdbcmd import BlastDbCmdRunner

def _make_runner(tmp_path: Path, **overrides) -> BlastDbCmdRunner:
    defaults = dict(
        working_dir = tmp_path,
        output_dir = tmp_path / "blast_reference",
        sample_id = "SAMPLE_01",
        accession = "NC_025767.1",
        database = "actinopterygii_mt_refs",
        logger = logging.getLogger("test"),
    )
    defaults.update(overrides)
    return BlastDbCmdRunner(**defaults)

def test_build_command(tmp_path: Path) -> None:
    runner = _make_runner(tmp_path)
    command = runner.build_command()
    assert command[0] == "blastdbcmd"
    assert "-db" in command and command[command.index("-db") + 1] == "actinopterygii_mt_refs"
    assert "-entry" in command and command[command.index("-entry") + 1] == "NC_025767.1"
    assert "-out" in command and command[command.index("-out") + 1] == str(runner.output_fasta)

def test_output_fasta_includes_sample_and_accession(tmp_path: Path) -> None:
    runner = _make_runner(tmp_path)
    assert runner.output_fasta.name == "SAMPLE_01.NC_025767.1.reference.fasta"

def test_output_fasta_sanitizes_unsafe_accession_characters(tmp_path: Path) -> None:
    # Older-style NCBI identifiers can contain '|' -- never embed raw
    # into a path.
    runner = _make_runner(tmp_path, accession = "gi|12345|ref|NC_000001.1|")
    assert "|" not in runner.output_fasta.name

def test_validate_inputs_rejects_empty_accession(tmp_path: Path) -> None:
    runner = _make_runner(tmp_path, accession = "")
    with pytest.raises(SystemExit):
        runner.validate_inputs()

def test_validate_inputs_rejects_empty_database(tmp_path: Path) -> None:
    runner = _make_runner(tmp_path, database = "")
    with pytest.raises(SystemExit):
        runner.validate_inputs()

def test_validate_inputs_passes_for_valid_accession_and_database(tmp_path: Path) -> None:
    runner = _make_runner(tmp_path)
    runner.validate_inputs()   # does not raise

def test_validate_outputs_rejects_missing_fasta(tmp_path: Path) -> None:
    runner = _make_runner(tmp_path)
    runner.output_dir.mkdir(parents = True)
    with pytest.raises(SystemExit):
        runner.validate_outputs()

def test_validate_outputs_rejects_empty_fasta(tmp_path: Path) -> None:
    runner = _make_runner(tmp_path)
    runner.output_dir.mkdir(parents = True)
    runner.output_fasta.touch()
    with pytest.raises(SystemExit):
        runner.validate_outputs()

def test_validate_outputs_passes_with_real_content(tmp_path: Path) -> None:
    runner = _make_runner(tmp_path)
    runner.output_dir.mkdir(parents = True)
    runner.output_fasta.write_text(">NC_025767.1\nACGT\n")
    runner.validate_outputs()   # does not raise
