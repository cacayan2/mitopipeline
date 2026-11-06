"""
test_run_rescue_mapping.py

Fixture tests for the rescue-mapping exec (#426): connects the rescue
trigger decision to the real minimap2/samtools chain, and records the
outcome for #406 to read next.
"""

# Imports
import importlib
import json
import sys
from pathlib import Path
import pytest
from mitopipeline.api.minimap2 import Minimap2Runner
from mitopipeline.api.samtools import SamtoolsRunner
from mitopipeline.models.command_result.command_result import CommandResult

run_rescue_mapping = importlib.import_module("mitopipeline.exec.07-assembly_rescue.run_rescue_mapping")

def _fake_result(success: bool = True, return_code: int = 0) -> CommandResult:
    from datetime import datetime
    now = datetime.now()
    return CommandResult(
        command = [], return_code = return_code, stdout = "", stderr = "", runtime_seconds = 0.0,
        success = success, tool_name = "stub", started_at = now, ended_at = now,
        peak_memory_mb = 0.0, user_cpu_seconds = 0.0, system_cpu_seconds = 0.0,
    )

def _stub_minimap2_run(self: Minimap2Runner) -> CommandResult:
    self.output_dir.mkdir(parents = True, exist_ok = True)
    self.output_sam.write_text("fake sam")
    return _fake_result()

def _stub_samtools_run(self: SamtoolsRunner) -> CommandResult:
    self.output_dir.mkdir(parents = True, exist_ok = True)
    if self.operation == "sort":
        self.sorted_bam.write_bytes(b"fake bam")
    elif self.operation == "index":
        self.bam_index.write_bytes(b"fake index")
    else:
        self.rescued_r1.write_text("fake r1")
        self.rescued_r2.write_text("fake r2")
    return _fake_result()

def _argv(tmp_path: Path, sample_id: str, status: str, attempt: int = 0, max_attempts: int = 2,
          recruitment: str = "candidate_mapped") -> list[str]:
    r1 = tmp_path / f"{sample_id}_R1.fastq.gz"
    r2 = tmp_path / f"{sample_id}_R2.fastq.gz"
    r1.touch(); r2.touch()
    reference = tmp_path / "reference.fasta"
    reference.write_text(">ref\nACGT\n")
    return [
        "run_rescue_mapping.py",
        "--sample-id", sample_id,
        "--r1", str(r1), "--r2", str(r2),
        "--reference-fasta", str(reference),
        "--recruitment", recruitment,
        "--status", status,
        "--attempt", str(attempt),
        "--max-attempts", str(max_attempts),
        "--output-dir", str(tmp_path / "rescue"),
        "--working-dir", str(tmp_path),
        "--log-file", str(tmp_path / "logs" / "rescue.log"),
        "--global-log-file", str(tmp_path / "logs" / "global.log"),
    ]

def test_circular_status_is_never_rescued(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    sample_id = "SAMPLE_01"
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, "circular"))
    assert run_rescue_mapping.main() == 0
    result = json.loads((tmp_path / "rescue" / f"{sample_id}.rescue_mapping.json").read_text())
    assert result["attempted"] is False
    assert result["rescued_r1"] is None

def test_declined_decision_records_the_real_status_it_was_given(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    # A declined attempt (e.g. already circular) never runs reassembly,
    # so its own decision record is the only place the real status
    # survives for the next attempt to read -- see
    # 07-assembly_rescue.smk's _previous_status_and_fasta().
    sample_id = "SAMPLE_01B"
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, "circular"))
    assert run_rescue_mapping.main() == 0
    result = json.loads((tmp_path / "rescue" / f"{sample_id}.rescue_mapping.json").read_text())
    assert result["status"] == "circular"

def test_fragmented_status_triggers_the_real_chain(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(Minimap2Runner, "run", _stub_minimap2_run)
    monkeypatch.setattr(SamtoolsRunner, "run", _stub_samtools_run)
    sample_id = "SAMPLE_02"
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, "fragmented"))

    assert run_rescue_mapping.main() == 0

    result = json.loads((tmp_path / "rescue" / f"{sample_id}.rescue_mapping.json").read_text())
    assert result["attempted"] is True
    assert result["rescued_r1"].endswith(f"{sample_id}.rescue.R1.fastq")
    assert result["rescued_r2"].endswith(f"{sample_id}.rescue.R2.fastq")
    assert Path(result["rescued_r1"]).is_file()

def test_attempts_exhausted_is_not_rescued(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    sample_id = "SAMPLE_03"
    # attempt 2 of max_attempts 2 -- every profile already tried.
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, "fragmented", attempt = 2, max_attempts = 2))
    assert run_rescue_mapping.main() == 0
    result = json.loads((tmp_path / "rescue" / f"{sample_id}.rescue_mapping.json").read_text())
    assert result["attempted"] is False
    assert "exhausted" in result["reason"]

def test_chain_failure_is_recorded_not_an_uncaught_crash(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    def _failing_run(self):
        raise SystemExit(1)
    monkeypatch.setattr(Minimap2Runner, "run", _failing_run)
    sample_id = "SAMPLE_04"
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, "fragmented"))

    assert run_rescue_mapping.main() == 1

    result = json.loads((tmp_path / "rescue" / f"{sample_id}.rescue_mapping.json").read_text())
    assert result["attempted"] is True
    assert result["rescued_r1"] is None

def test_a_failed_but_non_raising_stage_is_also_caught(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """BaseTool.run() does not raise when the underlying command
    exits nonzero -- it just returns a CommandResult with
    success=False, so a stage failure must be caught by checking
    that flag, the same as a SystemExit is caught."""
    monkeypatch.setattr(Minimap2Runner, "run", _stub_minimap2_run)
    monkeypatch.setattr(SamtoolsRunner, "run", lambda self: _fake_result(success = False, return_code = 1))
    sample_id = "SAMPLE_06"
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, "fragmented"))

    assert run_rescue_mapping.main() == 1

    result = json.loads((tmp_path / "rescue" / f"{sample_id}.rescue_mapping.json").read_text())
    assert result["attempted"] is True
    assert result["rescued_r1"] is None
    # These files must not exist just because the record claims
    # attempted=True -- nothing here is fabricated.
    assert not (tmp_path / "rescue" / f"{sample_id}.rescue.R1.fastq").exists()

def test_successful_chain_cleans_up_its_own_sam_and_bam_scratch(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """The raw SAM, sorted BAM, and its .bai index are intermediate
    mapping artifacts -- nothing downstream of rescue reads them back
    (confirmed against the whole codebase and this rule's own declared
    Snakemake outputs), so a successful chain must not leave them
    behind. Previously kept forever "for provenance"; found to total
    ~470GB across the live run's already-resolved samples."""
    monkeypatch.setattr(Minimap2Runner, "run", _stub_minimap2_run)
    monkeypatch.setattr(SamtoolsRunner, "run", _stub_samtools_run)
    sample_id = "SAMPLE_10"
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, "fragmented"))

    assert run_rescue_mapping.main() == 0

    rescue_dir = tmp_path / "rescue"
    assert not (rescue_dir / f"{sample_id}.rescue_mapping.sam").exists()
    assert not (rescue_dir / f"{sample_id}.rescue.sorted.bam").exists()
    assert not (rescue_dir / f"{sample_id}.rescue.sorted.bam.bai").exists()
    # The real output -- the extracted FASTQ -- must survive the cleanup.
    assert (rescue_dir / f"{sample_id}.rescue.R1.fastq").is_file()

def test_failed_chain_also_cleans_up_whatever_scratch_it_managed_to_write(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """A failure partway through (minimap2 and sort succeed, index
    fails) must not leave the SAM/BAM it did manage to produce behind
    either -- a failed attempt is just as useless to keep as a
    successful one once nothing downstream will ever read it."""
    monkeypatch.setattr(Minimap2Runner, "run", _stub_minimap2_run)

    def _sort_then_fail_index(self: SamtoolsRunner):
        if self.operation == "sort":
            self.output_dir.mkdir(parents = True, exist_ok = True)
            self.sorted_bam.write_bytes(b"fake bam")
            return _fake_result()
        return _fake_result(success = False, return_code = 1)
    monkeypatch.setattr(SamtoolsRunner, "run", _sort_then_fail_index)

    sample_id = "SAMPLE_11"
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, "fragmented"))

    assert run_rescue_mapping.main() == 1

    rescue_dir = tmp_path / "rescue"
    assert not (rescue_dir / f"{sample_id}.rescue_mapping.sam").exists()
    assert not (rescue_dir / f"{sample_id}.rescue.sorted.bam").exists()

def test_status_not_in_trigger_list_is_not_rescued(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    sample_id = "SAMPLE_05"
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, "missing"))
    assert run_rescue_mapping.main() == 0
    result = json.loads((tmp_path / "rescue" / f"{sample_id}.rescue_mapping.json").read_text())
    assert result["attempted"] is False

def test_full_pool_recruitment_skips_the_mapping_chain_entirely(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """full_pool must never touch Minimap2Runner/SamtoolsRunner at
    all -- stubbing them to raise proves the real code path really
    is skipped, not just coincidentally unexercised."""
    def _must_not_run(self):
        raise AssertionError("full_pool recruitment must not invoke the mapping chain")
    monkeypatch.setattr(Minimap2Runner, "run", _must_not_run)
    monkeypatch.setattr(SamtoolsRunner, "run", _must_not_run)

    sample_id = "SAMPLE_07"
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, "fragmented", recruitment = "full_pool"))

    assert run_rescue_mapping.main() == 0

    result = json.loads((tmp_path / "rescue" / f"{sample_id}.rescue_mapping.json").read_text())
    assert result["attempted"] is True
    assert result["recruitment"] == "full_pool"
    rescued_r1 = Path(result["rescued_r1"])
    assert rescued_r1.name == f"{sample_id}.rescue.R1.fastq"
    assert rescued_r1.is_file()

def test_full_pool_recruitment_decompresses_gzipped_input(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """The original trimmed reads are always .gz; the rescued-read
    naming convention (shared with candidate_mapped's own samtools
    fastq output) is always plain text -- full_pool must actually
    decompress on the way in, not just copy the bytes verbatim under
    a misleading .fastq extension."""
    sample_id = "SAMPLE_08"
    argv = _argv(tmp_path, sample_id, "fragmented", recruitment = "full_pool")
    r1_path = Path(argv[argv.index("--r1") + 1])
    r2_path = Path(argv[argv.index("--r2") + 1])
    import gzip
    with gzip.open(r1_path, "wb") as handle:
        handle.write(b"@read1\nACGTACGT\n+\nIIIIIIII\n")
    with gzip.open(r2_path, "wb") as handle:
        handle.write(b"@read1\nTTTTTTTT\n+\nIIIIIIII\n")
    monkeypatch.setattr(sys, "argv", argv)

    assert run_rescue_mapping.main() == 0

    result = json.loads((tmp_path / "rescue" / f"{sample_id}.rescue_mapping.json").read_text())
    rescued_r1 = Path(result["rescued_r1"])
    assert rescued_r1.read_text() == "@read1\nACGTACGT\n+\nIIIIIIII\n"

def test_full_pool_recruitment_respects_the_rescue_decision(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """A status that doesn't warrant rescue must still short-circuit
    before full_pool's own passthrough, same as candidate_mapped."""
    sample_id = "SAMPLE_09"
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, "circular", recruitment = "full_pool"))
    assert run_rescue_mapping.main() == 0
    result = json.loads((tmp_path / "rescue" / f"{sample_id}.rescue_mapping.json").read_text())
    assert result["attempted"] is False
    assert result["rescued_r1"] is None
