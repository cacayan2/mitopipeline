"""
test_minimap2.py

Unit tests for the Minimap2Runner class (#423).
"""

# Imports
import logging
from pathlib import Path
import pytest
from mitopipeline.api.minimap2 import Minimap2Runner
from mitopipeline.models.sample.sample import Sample

@pytest.fixture
def paired_sample(tmp_path: Path) -> Sample:
    r1 = tmp_path / "rescue_R1.fastq.gz"
    r2 = tmp_path / "rescue_R2.fastq.gz"
    r1.touch()
    r2.touch()
    return Sample(sample_id = "SAMPLE_01", r1 = r1, r2 = r2)

@pytest.fixture
def reference_fasta(tmp_path: Path) -> Path:
    ref = tmp_path / "reference.fasta"
    ref.write_text(">ref\nACGT\n")
    return ref

def _make_runner(tmp_path: Path, sample: Sample, reference_fasta: Path, **overrides) -> Minimap2Runner:
    defaults = dict(
        working_dir = tmp_path,
        output_dir = tmp_path / "mapping",
        sample = sample,
        reference_fasta = reference_fasta,
        logger = logging.getLogger("test"),
    )
    defaults.update(overrides)
    return Minimap2Runner(**defaults)

def test_build_command_wires_up_reference_and_reads(tmp_path: Path, paired_sample: Sample, reference_fasta: Path) -> None:
    runner = _make_runner(tmp_path, paired_sample, reference_fasta, threads = 6)
    command = runner.build_command()
    assert command[0] == "minimap2"
    assert "-ax" in command and command[command.index("-ax") + 1] == "sr"
    assert "-t" in command and command[command.index("-t") + 1] == "6"
    assert str(reference_fasta) in command
    assert str(paired_sample.r1) in command
    assert str(paired_sample.r2) in command
    assert str(runner.output_sam).endswith("SAMPLE_01.rescue_mapping.sam")

def test_build_command_creates_output_directory(tmp_path: Path, paired_sample: Sample, reference_fasta: Path) -> None:
    runner = _make_runner(tmp_path, paired_sample, reference_fasta)
    assert not runner.output_dir.exists()
    runner.build_command()
    assert runner.output_dir.is_dir()

def test_validate_inputs_rejects_missing_reference(tmp_path: Path, paired_sample: Sample) -> None:
    runner = _make_runner(tmp_path, paired_sample, tmp_path / "missing_reference.fasta")
    with pytest.raises(SystemExit):
        runner.validate_inputs()

def test_validate_inputs_rejects_missing_reads(tmp_path: Path, reference_fasta: Path) -> None:
    sample = Sample(sample_id = "SAMPLE_02", r1 = tmp_path / "missing_R1.fastq.gz", r2 = tmp_path / "missing_R2.fastq.gz")
    runner = _make_runner(tmp_path, sample, reference_fasta)
    with pytest.raises(SystemExit):
        runner.validate_inputs()

def test_validate_inputs_rejects_nonpositive_threads(tmp_path: Path, paired_sample: Sample, reference_fasta: Path) -> None:
    runner = _make_runner(tmp_path, paired_sample, reference_fasta, threads = 0)
    with pytest.raises(SystemExit):
        runner.validate_inputs()

def test_validate_inputs_passes_with_real_files(tmp_path: Path, paired_sample: Sample, reference_fasta: Path) -> None:
    runner = _make_runner(tmp_path, paired_sample, reference_fasta)
    runner.validate_inputs()  # must not raise

def test_validate_outputs_rejects_missing_sam(tmp_path: Path, paired_sample: Sample, reference_fasta: Path) -> None:
    runner = _make_runner(tmp_path, paired_sample, reference_fasta)
    runner.output_dir.mkdir(parents = True)
    with pytest.raises(SystemExit):
        runner.validate_outputs()

def test_validate_outputs_rejects_empty_sam(tmp_path: Path, paired_sample: Sample, reference_fasta: Path) -> None:
    runner = _make_runner(tmp_path, paired_sample, reference_fasta)
    runner.output_dir.mkdir(parents = True)
    runner.output_sam.touch()
    with pytest.raises(SystemExit):
        runner.validate_outputs()

def test_validate_outputs_passes_with_real_sam(tmp_path: Path, paired_sample: Sample, reference_fasta: Path) -> None:
    runner = _make_runner(tmp_path, paired_sample, reference_fasta)
    runner.output_dir.mkdir(parents = True)
    runner.output_sam.write_text("@HD\tVN:1.6\n")
    runner.validate_outputs()  # must not raise

def test_custom_preset_is_used_instead_of_default(tmp_path: Path, paired_sample: Sample, reference_fasta: Path) -> None:
    runner = _make_runner(tmp_path, paired_sample, reference_fasta, preset = "map-ont")
    command = runner.build_command()
    assert command[command.index("-ax") + 1] == "map-ont"
