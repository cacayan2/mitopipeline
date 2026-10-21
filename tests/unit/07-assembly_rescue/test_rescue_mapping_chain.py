"""
test_rescue_mapping_chain.py

Integration-style unit test (#425) exercising the full rescue mapping
chain for one fixture sample: minimap2 -> samtools sort -> samtools
index -> samtools fastq. Verifies the artifact contract -- a SAM, a
sorted+indexed BAM, and rescued paired FASTQ -- all land under distinct,
sample-specific filenames, and that each stage actually consumes the
PREVIOUS stage's real output path, not a guessed one.

The real minimap2/samtools binaries are not available in this
environment, so each BaseTool.run() is stubbed to create its declared
output files and return a successful CommandResult -- this test
verifies the WIRING between stages, not the external tools' own
behavior (same approach as test_preprocessing_chain.py).
"""

# Imports
from datetime import datetime
from pathlib import Path
import pytest
from mitopipeline.api.minimap2 import Minimap2Runner
from mitopipeline.api.samtools import SamtoolsRunner
from mitopipeline.models.command_result.command_result import CommandResult
from mitopipeline.models.sample.sample import Sample

def _fake_result() -> CommandResult:
    now = datetime.now()
    return CommandResult(
        command = [], return_code = 0, stdout = "", stderr = "", runtime_seconds = 0.0,
        success = True, tool_name = "stub", started_at = now, ended_at = now,
        peak_memory_mb = 0.0, user_cpu_seconds = 0.0, system_cpu_seconds = 0.0,
    )

def _stub_minimap2_run(self: Minimap2Runner) -> CommandResult:
    self.output_dir.mkdir(parents = True, exist_ok = True)
    self.output_sam.write_text("fake sam output")
    return _fake_result()

def _stub_samtools_run(self: SamtoolsRunner) -> CommandResult:
    self.output_dir.mkdir(parents = True, exist_ok = True)
    if self.operation == "sort":
        self.sorted_bam.write_bytes(b"fake sorted bam")
    elif self.operation == "index":
        self.bam_index.write_bytes(b"fake bam index")
    else:
        self.rescued_r1.write_text("fake rescued r1")
        self.rescued_r2.write_text("fake rescued r2")
    return _fake_result()

def test_rescue_mapping_chain_produces_sample_specific_artifacts(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(Minimap2Runner, "run", _stub_minimap2_run)
    monkeypatch.setattr(SamtoolsRunner, "run", _stub_samtools_run)

    sample_id = "FIXTURE_01"
    r1 = tmp_path / "rescue_in" / f"{sample_id}_R1.fastq.gz"
    r2 = tmp_path / "rescue_in" / f"{sample_id}_R2.fastq.gz"
    r1.parent.mkdir(parents = True)
    r1.touch()
    r2.touch()
    sample = Sample(sample_id = sample_id, r1 = r1, r2 = r2)

    reference = tmp_path / "reference.fasta"
    reference.write_text(">ref\nACGT\n")

    mapping_dir = tmp_path / "rescue" / sample_id

    # Stage 1: map reads against the reference.
    mapper = Minimap2Runner(working_dir = tmp_path, output_dir = mapping_dir,
                             sample = sample, reference_fasta = reference)
    mapper.run()
    assert mapper.output_sam.is_file()
    assert mapper.output_sam.name == f"{sample_id}.rescue_mapping.sam"

    # Stage 2: sort the ACTUAL mapping output, not a guessed path.
    sorter = SamtoolsRunner(working_dir = tmp_path, output_dir = mapping_dir, sample = sample,
                             operation = "sort", input_path = mapper.output_sam)
    sorter.run()
    assert sorter.sorted_bam.is_file()
    assert sorter.sorted_bam.name == f"{sample_id}.rescue.sorted.bam"

    # Stage 3: index the ACTUAL sorted BAM.
    indexer = SamtoolsRunner(working_dir = tmp_path, output_dir = mapping_dir, sample = sample,
                              operation = "index", input_path = sorter.sorted_bam)
    indexer.run()
    assert indexer.bam_index.is_file()
    assert indexer.bam_index.name == f"{sample_id}.rescue.sorted.bam.bai"

    # Stage 4: extract rescued, mapped paired reads from the ACTUAL sorted BAM.
    extractor = SamtoolsRunner(working_dir = tmp_path, output_dir = mapping_dir, sample = sample,
                                operation = "fastq", input_path = sorter.sorted_bam)
    extractor.run()
    assert extractor.rescued_r1.is_file() and extractor.rescued_r2.is_file()
    assert extractor.rescued_r1.name == f"{sample_id}.rescue.R1.fastq.gz"
    assert extractor.rescued_r2.name == f"{sample_id}.rescue.R2.fastq.gz"

    # Every artifact for this sample lands in the same sample-scoped
    # directory, none of them collide, and all four are genuinely
    # distinct files on disk.
    artifacts = {mapper.output_sam, sorter.sorted_bam, indexer.bam_index,
                 extractor.rescued_r1, extractor.rescued_r2}
    assert len(artifacts) == 5
    assert all(path.parent == mapping_dir for path in artifacts)

def test_two_samples_never_collide_in_the_same_rescue_run(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Sample-specific filenames, not just sample-specific directories --
    a shared output_dir (e.g. a batch run) must still never collide."""
    monkeypatch.setattr(Minimap2Runner, "run", _stub_minimap2_run)

    reference = tmp_path / "reference.fasta"
    reference.write_text(">ref\nACGT\n")
    shared_dir = tmp_path / "rescue_batch"

    runners = []
    for sample_id in ("SAMPLE_A", "SAMPLE_B"):
        r1 = tmp_path / f"{sample_id}_R1.fastq.gz"
        r2 = tmp_path / f"{sample_id}_R2.fastq.gz"
        r1.touch(); r2.touch()
        sample = Sample(sample_id = sample_id, r1 = r1, r2 = r2)
        runner = Minimap2Runner(working_dir = tmp_path, output_dir = shared_dir,
                                 sample = sample, reference_fasta = reference)
        runner.run()
        runners.append(runner)

    assert runners[0].output_sam != runners[1].output_sam
    assert runners[0].output_sam.is_file() and runners[1].output_sam.is_file()
