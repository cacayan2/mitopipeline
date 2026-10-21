"""
test_samtools.py

Unit tests for the SamtoolsRunner class (#424): sort/index/fastq
command construction and missing-BAM failure cases.
"""

# Imports
import logging
from pathlib import Path
import pytest
from mitopipeline.api.samtools import SamtoolsRunner
from mitopipeline.models.sample.sample import Sample

@pytest.fixture
def sample() -> Sample:
    return Sample(sample_id = "SAMPLE_01", r1 = Path("unused_R1.fastq.gz"), r2 = Path("unused_R2.fastq.gz"))

def _make_runner(tmp_path: Path, sample: Sample, operation: str, input_path: Path, **overrides) -> SamtoolsRunner:
    defaults = dict(
        working_dir = tmp_path,
        output_dir = tmp_path / "rescue",
        sample = sample,
        operation = operation,
        input_path = input_path,
        logger = logging.getLogger("test"),
    )
    defaults.update(overrides)
    return SamtoolsRunner(**defaults)

# ---- sort ----

def test_sort_build_command(tmp_path: Path, sample: Sample) -> None:
    sam = tmp_path / "mapping.sam"
    sam.write_text("@HD\n")
    runner = _make_runner(tmp_path, sample, "sort", sam, threads = 4)
    command = runner.build_command()
    assert command[:2] == ["samtools", "sort"]
    assert "-@" in command and command[command.index("-@") + 1] == "4"
    assert "-o" in command and command[command.index("-o") + 1] == str(runner.sorted_bam)
    assert command[-1] == str(sam)

def test_sort_validate_outputs_rejects_missing_bam(tmp_path: Path, sample: Sample) -> None:
    sam = tmp_path / "mapping.sam"
    sam.write_text("@HD\n")
    runner = _make_runner(tmp_path, sample, "sort", sam)
    runner.output_dir.mkdir(parents = True)
    with pytest.raises(SystemExit):
        runner.validate_outputs()

# ---- index ----

def test_index_build_command(tmp_path: Path, sample: Sample) -> None:
    bam = tmp_path / "sorted.bam"
    bam.write_bytes(b"BAM")
    runner = _make_runner(tmp_path, sample, "index", bam)
    command = runner.build_command()
    assert command[:2] == ["samtools", "index"]
    assert command[-1] == str(bam)

def test_index_validate_outputs_rejects_missing_index(tmp_path: Path, sample: Sample) -> None:
    bam = tmp_path / "sorted.bam"
    bam.write_bytes(b"BAM")
    runner = _make_runner(tmp_path, sample, "index", bam)
    runner.output_dir.mkdir(parents = True)
    with pytest.raises(SystemExit):
        runner.validate_outputs()

def test_index_validate_outputs_passes_with_real_index(tmp_path: Path, sample: Sample) -> None:
    bam = tmp_path / "sorted.bam"
    bam.write_bytes(b"BAM")
    runner = _make_runner(tmp_path, sample, "index", bam)
    runner.output_dir.mkdir(parents = True)
    runner.bam_index.write_bytes(b"BAI")
    runner.validate_outputs()  # must not raise

# ---- fastq (rescued-read extraction) ----

def test_fastq_build_command_discards_unpaired_reads(tmp_path: Path, sample: Sample) -> None:
    bam = tmp_path / "sorted.bam"
    bam.write_bytes(b"BAM")
    runner = _make_runner(tmp_path, sample, "fastq", bam)
    command = runner.build_command()
    assert command[:2] == ["samtools", "fastq"]
    assert "-n" in command  # keep paired read names intact
    assert "-1" in command and command[command.index("-1") + 1] == str(runner.rescued_r1)
    assert "-2" in command and command[command.index("-2") + 1] == str(runner.rescued_r2)
    assert "-0" in command  # unpaired reads discarded
    assert "-s" in command  # singletons discarded
    assert command[-1] == str(bam)

def test_fastq_validate_outputs_requires_both_mates(tmp_path: Path, sample: Sample) -> None:
    bam = tmp_path / "sorted.bam"
    bam.write_bytes(b"BAM")
    runner = _make_runner(tmp_path, sample, "fastq", bam)
    runner.output_dir.mkdir(parents = True)
    runner.rescued_r1.write_text("@read\nACGT\n+\nIIII\n")
    # r2 intentionally left missing
    with pytest.raises(SystemExit):
        runner.validate_outputs()

def test_fastq_validate_outputs_passes_with_both_mates(tmp_path: Path, sample: Sample) -> None:
    bam = tmp_path / "sorted.bam"
    bam.write_bytes(b"BAM")
    runner = _make_runner(tmp_path, sample, "fastq", bam)
    runner.output_dir.mkdir(parents = True)
    runner.rescued_r1.write_text("@read\nACGT\n+\nIIII\n")
    runner.rescued_r2.write_text("@read\nACGT\n+\nIIII\n")
    runner.validate_outputs()  # must not raise

# ---- shared validation / failure cases ----

def test_validate_inputs_rejects_unknown_operation(tmp_path: Path, sample: Sample) -> None:
    bam = tmp_path / "sorted.bam"
    bam.write_bytes(b"BAM")
    runner = _make_runner(tmp_path, sample, "merge", bam)  # not a supported operation
    with pytest.raises(SystemExit):
        runner.validate_inputs()

def test_validate_inputs_rejects_missing_bam(tmp_path: Path, sample: Sample) -> None:
    """The missing-BAM failure case #424 specifically asks for."""
    runner = _make_runner(tmp_path, sample, "index", tmp_path / "does_not_exist.bam")
    with pytest.raises(SystemExit):
        runner.validate_inputs()

def test_validate_inputs_rejects_nonpositive_threads(tmp_path: Path, sample: Sample) -> None:
    bam = tmp_path / "sorted.bam"
    bam.write_bytes(b"BAM")
    runner = _make_runner(tmp_path, sample, "index", bam, threads = -1)
    with pytest.raises(SystemExit):
        runner.validate_inputs()
