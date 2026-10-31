"""
test_run_rescue_reassembly.py

Fixture tests for the rescue-reassembly exec (#406): re-runs
GetOrganelle on rescued reads with one getorganelle_profiles entry's
tuning, through the same GetOrganelleRunner the primary attempt uses.
"""

# Imports
import importlib
import json
import sys
from datetime import datetime
from pathlib import Path
import pytest
from mitopipeline.api.getorganelle import GetOrganelleRunner
from mitopipeline.models.command_result.command_result import CommandResult

run_rescue_reassembly = importlib.import_module("mitopipeline.exec.07-assembly_rescue.run_rescue_reassembly")

def _fake_result(success: bool = True, return_code: int = 0) -> CommandResult:
    now = datetime.now()
    return CommandResult(
        command = [], return_code = return_code, stdout = "", stderr = "", runtime_seconds = 1.5,
        success = success, tool_name = "getorganelle", started_at = now, ended_at = now,
        peak_memory_mb = 0.0, user_cpu_seconds = 0.0, system_cpu_seconds = 0.0,
    )

def _argv(tmp_path: Path, sample_id: str, output_dir: Path, **extra) -> list[str]:
    r1 = tmp_path / f"{sample_id}.rescue.R1.fastq"
    r2 = tmp_path / f"{sample_id}.rescue.R2.fastq"
    r1.touch(); r2.touch()
    argv = [
        "run_rescue_reassembly.py",
        "--sample-id", sample_id,
        "--r1", str(r1), "--r2", str(r2),
        "--profile-name", extra.get("profile_name", "reduced_coverage"),
        "--output-dir", str(output_dir),
        "--working-dir", str(tmp_path),
        "--threads", "2",
        "--log-file", str(tmp_path / "logs" / "rescue_reassembly.log"),
        "--global-log-file", str(tmp_path / "logs" / "global.log"),
    ]
    if "kmer" in extra:
        argv += ["--kmer", extra["kmer"]]
    if "target_coverage" in extra:
        argv += ["--target-coverage", str(extra["target_coverage"])]
    if "seed_file" in extra:
        argv += ["--seed-file", str(extra["seed_file"])]
    if extra.get("reseed"):
        argv += ["--reseed"]
    if "exclude_gene" in extra:
        argv += ["--exclude-gene", extra["exclude_gene"]]
    if "scientific_name" in extra:
        argv += ["--scientific-name", extra["scientific_name"]]
    if "database_dir" in extra:
        argv += ["--database-dir", str(extra["database_dir"])]
    return argv

def _stub_run_success(self) -> CommandResult:
    self.output_dir.mkdir(parents = True, exist_ok = True)
    (self.output_dir / "k85.complete.path_sequence.fasta").write_text(">seq\nACGT\n")
    (self.output_dir / "k85.complete.selected_graph.gfa").write_text("H\n")
    from mitopipeline.api.getorganelle import discover_assembly_artifacts
    self.assembly_status = discover_assembly_artifacts(self.output_dir, self.logger)
    return _fake_result()

def test_successful_rescue_writes_status_with_profile_name(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    sample_id = "SAMPLE_01"
    output_dir = tmp_path / "rescue_attempt_1"

    def _stub_run(self: GetOrganelleRunner) -> CommandResult:
        self.output_dir.mkdir(parents = True, exist_ok = True)
        (self.output_dir / "k85.complete.path_sequence.fasta").write_text(">seq\nACGT\n")
        (self.output_dir / "k85.complete.selected_graph.gfa").write_text("H\n")
        from mitopipeline.api.getorganelle import discover_assembly_artifacts
        self.assembly_status = discover_assembly_artifacts(self.output_dir, self.logger)
        return _fake_result()

    monkeypatch.setattr(GetOrganelleRunner, "run", _stub_run)
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, output_dir,
                                            kmer = "21,45,65,85,105", target_coverage = 100))

    assert run_rescue_reassembly.main() == 0
    status = json.loads((output_dir / f"{sample_id}.assembly_status.json").read_text())
    assert status == {"sample_id": sample_id, "status": "circular", "success": True,
                       "profile_name": "reduced_coverage"}

def test_tool_options_carry_the_profile_tuning_through(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    sample_id = "SAMPLE_02"
    output_dir = tmp_path / "rescue_attempt_1"
    seen = {}

    def _stub_run(self: GetOrganelleRunner) -> CommandResult:
        seen.update(self.tool_options)
        self.output_dir.mkdir(parents = True, exist_ok = True)
        (self.output_dir / "k85.path_sequence.fasta").write_text(">seq\nACGT\n")
        (self.output_dir / "k85.selected_graph.gfa").write_text("H\n")
        from mitopipeline.api.getorganelle import discover_assembly_artifacts
        self.assembly_status = discover_assembly_artifacts(self.output_dir, self.logger)
        return _fake_result()

    monkeypatch.setattr(GetOrganelleRunner, "run", _stub_run)
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, output_dir,
                                            kmer = "21,33,55,85,115", target_coverage = 100))

    assert run_rescue_reassembly.main() == 0
    assert seen.get("kmer") == "21,33,55,85,115"
    assert seen.get("target_coverage") == 100.0

def test_no_native_candidate_is_classified_missing_and_fails(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    sample_id = "SAMPLE_03"
    output_dir = tmp_path / "rescue_attempt_1"

    def _stub_run(self: GetOrganelleRunner) -> CommandResult:
        self.output_dir.mkdir(parents = True, exist_ok = True)
        raise SystemExit(1)

    monkeypatch.setattr(GetOrganelleRunner, "run", _stub_run)
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, output_dir))

    assert run_rescue_reassembly.main() == 1
    status = json.loads((output_dir / f"{sample_id}.assembly_status.json").read_text())
    assert status["status"] == "missing"
    assert status["success"] is False

def test_two_rescue_attempts_never_collide(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Different attempts must write to different output dirs -- the
    caller's job (a Snakemake rule/orchestrator), but this confirms
    the script itself never assumes or hardcodes a shared path."""
    sample_id = "SAMPLE_04"
    attempt1_dir = tmp_path / "rescue_attempt_1"
    attempt2_dir = tmp_path / "rescue_attempt_2"

    def _stub_run(self: GetOrganelleRunner) -> CommandResult:
        self.output_dir.mkdir(parents = True, exist_ok = True)
        (self.output_dir / "k85.path_sequence.fasta").write_text(">seq\nACGT\n")
        (self.output_dir / "k85.selected_graph.gfa").write_text("H\n")
        from mitopipeline.api.getorganelle import discover_assembly_artifacts
        self.assembly_status = discover_assembly_artifacts(self.output_dir, self.logger)
        return _fake_result()

    monkeypatch.setattr(GetOrganelleRunner, "run", _stub_run)

    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, attempt1_dir, profile_name = "reduced_coverage"))
    assert run_rescue_reassembly.main() == 0
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, attempt2_dir, profile_name = "alternate_kmers"))
    assert run_rescue_reassembly.main() == 0

    status1 = json.loads((attempt1_dir / f"{sample_id}.assembly_status.json").read_text())
    status2 = json.loads((attempt2_dir / f"{sample_id}.assembly_status.json").read_text())
    assert status1["profile_name"] == "reduced_coverage"
    assert status2["profile_name"] == "alternate_kmers"

def test_seed_file_is_reused_and_recorded(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """candidate_guided/full_pool_relaxed_graph both reuse the primary
    attempt's own seed as-is -- no reseed logic involved."""
    sample_id = "SAMPLE_10"
    output_dir = tmp_path / "rescue_attempt_1"
    seed = tmp_path / "primary_seed.fasta"
    seed.write_text(">seed\nACGT\n")
    seen = {}

    def _stub_run(self):
        seen.update(self.tool_options)
        return _stub_run_success(self)

    monkeypatch.setattr(GetOrganelleRunner, "run", _stub_run)
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, output_dir, seed_file = seed))

    assert run_rescue_reassembly.main() == 0
    assert seen.get("seed_file") == str(seed)
    status = json.loads((output_dir / f"{sample_id}.assembly_status.json").read_text())
    assert status["seed_file"] == str(seed)

def test_reseed_derives_a_different_locus_and_uses_it(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """full_pool_alternate_seed's own real behavior: a different real
    locus is selected (never the excluded one) and actually handed to
    GetOrganelle as the seed."""
    sample_id = "SAMPLE_11"
    output_dir = tmp_path / "rescue_attempt_2"
    reseeded_fasta = tmp_path / "reseeded.fasta"
    reseeded_fasta.write_text(">reseed\nACGT\n")

    from mitopipeline.utils.seed.models import SeedSelectionResult
    def _fake_select_seed(sample_id, scientific_name, database_dir, output_dir, logger, loci):
        assert "COXI" not in loci  # the excluded gene must never be retried
        return SeedSelectionResult(sample_id = sample_id, scientific_name = scientific_name, gene = loci[0],
                                   accession = "FAKE123", source = "mitofish", fasta_path = reseeded_fasta,
                                   fallback_used = False)

    seen = {}
    def _stub_run(self):
        seen.update(self.tool_options)
        return _stub_run_success(self)

    monkeypatch.setattr(run_rescue_reassembly, "select_seed", _fake_select_seed)
    monkeypatch.setattr(GetOrganelleRunner, "run", _stub_run)
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, output_dir, reseed = True,
                                            exclude_gene = "COXI", scientific_name = "Danio rerio",
                                            database_dir = tmp_path / "mitofish"))

    assert run_rescue_reassembly.main() == 0
    assert seen.get("seed_file") == str(reseeded_fasta)
    record = json.loads((output_dir / f"{sample_id}.rescue_seed_selection.json").read_text())
    assert record["excluded_gene"] == "COXI"
    assert record["gene"] != "COXI"

def test_reseed_without_scientific_name_falls_back_to_unseeded(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """A misconfigured reseed profile (missing --scientific-name/
    --database-dir) must degrade to a real unseeded run, never crash
    the whole attempt over what the primary attempt's own seed
    selection already proved is non-fatal."""
    sample_id = "SAMPLE_12"
    output_dir = tmp_path / "rescue_attempt_2"
    seen = {}

    def _stub_run(self):
        seen.update(self.tool_options)
        return _stub_run_success(self)

    monkeypatch.setattr(GetOrganelleRunner, "run", _stub_run)
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, output_dir, reseed = True, exclude_gene = "COXI"))

    assert run_rescue_reassembly.main() == 0
    assert "seed_file" not in seen

def test_reseed_takes_priority_over_an_explicit_seed_file(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """If a caller ever passes both (shouldn't happen from the real
    Snakemake rule, which treats them as mutually exclusive per
    profile), --reseed must win -- it's the more specific, attempt-
    local instruction."""
    sample_id = "SAMPLE_13"
    output_dir = tmp_path / "rescue_attempt_2"
    stale_seed = tmp_path / "stale.fasta"
    stale_seed.write_text(">stale\nACGT\n")
    reseeded_fasta = tmp_path / "reseeded.fasta"
    reseeded_fasta.write_text(">reseed\nACGT\n")

    from mitopipeline.utils.seed.models import SeedSelectionResult
    def _fake_select_seed(sample_id, scientific_name, database_dir, output_dir, logger, loci):
        return SeedSelectionResult(sample_id = sample_id, scientific_name = scientific_name, gene = loci[0],
                                   accession = "FAKE123", source = "mitofish", fasta_path = reseeded_fasta,
                                   fallback_used = False)

    seen = {}
    def _stub_run(self):
        seen.update(self.tool_options)
        return _stub_run_success(self)

    monkeypatch.setattr(run_rescue_reassembly, "select_seed", _fake_select_seed)
    monkeypatch.setattr(GetOrganelleRunner, "run", _stub_run)
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, output_dir, seed_file = stale_seed, reseed = True,
                                            scientific_name = "Danio rerio", database_dir = tmp_path / "mitofish"))

    assert run_rescue_reassembly.main() == 0
    assert seen.get("seed_file") == str(reseeded_fasta)
