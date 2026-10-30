"""
test_table2asn.py

Unit tests for the Table2AsnRunner class: builds the real table2asn
command for its own template/FASTA/feature-table/output flags, and
rejects a basename mismatch between the FASTA and feature table
before ever shelling out.
"""

# Imports
import logging
from pathlib import Path
import pytest
from mitopipeline.api.table2asn import Table2AsnRunner

def _make_runner(tmp_path: Path, **overrides) -> Table2AsnRunner:
    fasta = tmp_path / "sample1.fsa"
    fasta.write_text(">sample1\nACGT\n")
    tbl = tmp_path / "sample1.tbl"
    tbl.write_text(">Feature sample1\n1\t4\tsource\n")
    template = tmp_path / "template.sbt"
    template.write_text("Submit-block ::= {}\n")
    defaults = dict(
        working_dir = tmp_path, fasta_path = fasta, feature_table_path = tbl,
        submission_template_path = template, output_basename = tmp_path / "output" / "sample1",
        logger = logging.getLogger("test"),
    )
    defaults.update(overrides)
    return Table2AsnRunner(**defaults)

def test_build_command_includes_template_fasta_tbl_and_output(tmp_path: Path) -> None:
    runner = _make_runner(tmp_path)
    command = runner.build_command()
    assert command[0] == "table2asn"
    assert "-t" in command and command[command.index("-t") + 1] == str(runner.submission_template_path.resolve())
    assert "-i" in command and command[command.index("-i") + 1] == str(runner.fasta_path.resolve())
    assert "-f" in command and command[command.index("-f") + 1] == str(runner.feature_table_path.resolve())
    assert "-o" in command and command[command.index("-o") + 1] == str(runner.output_basename.with_suffix(".sqn").resolve())
    assert "-V" in command and command[command.index("-V") + 1] == "vb"
    assert "-Z" in command

def test_build_command_creates_the_output_directory(tmp_path: Path) -> None:
    runner = _make_runner(tmp_path)
    runner.build_command()
    assert runner.output_basename.parent.is_dir()

def test_version_command_uses_the_configured_binary_path(tmp_path: Path) -> None:
    runner = _make_runner(tmp_path, binary_path = "/opt/bin/table2asn")
    assert runner.version_command() == ["/opt/bin/table2asn", "-version"]

def test_validate_inputs_rejects_mismatched_basenames(tmp_path: Path) -> None:
    other_tbl = tmp_path / "different_name.tbl"
    other_tbl.write_text(">Feature x\n1\t4\tsource\n")
    runner = _make_runner(tmp_path, feature_table_path = other_tbl)
    with pytest.raises(SystemExit):
        runner.validate_inputs()

def test_validate_inputs_rejects_missing_template(tmp_path: Path) -> None:
    runner = _make_runner(tmp_path, submission_template_path = tmp_path / "missing.sbt")
    with pytest.raises(SystemExit):
        runner.validate_inputs()
