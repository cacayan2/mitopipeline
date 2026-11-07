"""
test_pipeline_launcher.py

Unit tests for the pipeline launcher's config loading: load_config()
runs before any PipelineJob exists (the job's own output_root/job_id
come from the config itself), so it must never require one. Also
covers _job_configuration() -- every real config file nests
output_root/manifest/job_id/input_directory under its own
"job_configuration:" heading, not at the top level -- and
run_snakemake()'s own command construction.
"""

# Imports
import logging
from pathlib import Path
import subprocess
import sys
import pytest
from mitopipeline.launcher import pipeline_launcher
from mitopipeline.launcher.pipeline_launcher import (
    _configured_job_id, _job_configuration, _rerender_run_report_with_final_status,
    discover_fastq_pairs, load_config, parse_args, prepare_runtime_manifest, run_snakemake,
    write_validated_manifest,
)
from mitopipeline.models.sample.sample import Sample

def test_loads_a_real_yaml_config_without_a_job(tmp_path: Path) -> None:
    config_path = tmp_path / "config.yaml"
    config_path.write_text("job_configuration:\n  job_id: test_job\n  output_root: /tmp/out\n")

    config = load_config(config_path)
    assert config["job_configuration"]["job_id"] == "test_job"

def test_missing_config_file_exits_without_a_job(tmp_path: Path) -> None:
    with pytest.raises(SystemExit):
        load_config(tmp_path / "does_not_exist.yaml")

def test_empty_config_file_exits_without_a_job(tmp_path: Path) -> None:
    config_path = tmp_path / "config.yaml"
    config_path.write_text("")

    with pytest.raises(SystemExit):
        load_config(config_path)

def test_job_configuration_reads_the_nested_block() -> None:
    config = {"job_configuration": {"output_root": "/tmp/out", "job_id": "j1"}, "stages": {}}
    assert _job_configuration(config) == {"output_root": "/tmp/out", "job_id": "j1"}

def test_job_configuration_falls_back_to_a_flat_config() -> None:
    config = {"output_root": "/tmp/out", "job_id": "j1"}
    assert _job_configuration(config) == config

def test_configured_job_id_is_importable_and_normalizes_blank_values() -> None:
    # Regression: this function used to be defined after the module's
    # own "if __name__ == '__main__'" guard, so running this file as
    # a script raised NameError before ever reaching main().
    assert _configured_job_id(None) is None
    assert _configured_job_id("null") is None
    assert _configured_job_id("  ") is None
    assert _configured_job_id("real_job") == "real_job"

def test_conda_is_enabled_by_default_without_the_no_conda_flag(monkeypatch: pytest.MonkeyPatch) -> None:
    # Regression: --no-conda used action="store_false" on itself,
    # which inverted its own meaning -- conda was off by default and
    # only turned on if you passed the flag that's supposed to
    # disable it. Every real Snakemake invocation from this launcher
    # silently ran without --use-conda, so conda-managed tool
    # binaries (GetOrganelle, BLAST, MITOS2) were never on $PATH.
    monkeypatch.setattr(sys, "argv", ["pipeline_launcher.py", "--config", "config.yaml"])
    args = parse_args()
    assert args.no_conda is False
    assert not args.no_conda  # use_conda = not args.no_conda must be True by default

def test_no_conda_flag_disables_conda(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(sys, "argv", ["pipeline_launcher.py", "--config", "config.yaml", "--no-conda"])
    args = parse_args()
    assert args.no_conda is True

class _FakeProcess:
    def __init__(self, *args, **kwargs) -> None:
        self.stdout = iter(["dry-run output\n"])
        self.returncode = 0
    def wait(self) -> None:
        return None

def test_run_snakemake_resolves_the_binary_against_sys_executable(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    # Regression: a bare "snakemake" command string only resolves via
    # $PATH, which subprocess.Popen does not get for free just because
    # this launcher was itself invoked through a venv's python binary.
    captured_commands: list[list[str]] = []

    def _fake_popen(command: list[str], **kwargs) -> _FakeProcess:
        captured_commands.append(command)
        return _FakeProcess()

    monkeypatch.setattr(subprocess, "Popen", _fake_popen)

    run_snakemake(
        runtime_config = tmp_path / "runtime_config.yaml", cores = 2, use_conda = True,
        logger = logging.getLogger("test"), concurrent_slots = 8,
    )

    command = captured_commands[0]
    assert command[0] == str(Path(sys.executable).parent / "snakemake")
    assert "--conda-not-block-search-path-envvars" in command
    # Regression: "entrez = 1" (with spaces) is not a valid Snakemake
    # --resources identifier and raised a real WorkflowError.
    assert command[command.index("--resources") + 1] == "entrez=1"
    # Regression (2026-11-06): resuming with a large backlog let a
    # burst of GetOrganelle-class jobs (primary assembly, rescue
    # reassembly) launch within seconds of each other, and the kernel
    # logged real memory pressure while many reassembly attempts died
    # before writing any log at all. This caps their concurrency
    # regardless of idle --cores, matching the rules' own
    # `resources: getorganelle_slots = 1`. concurrent_slots is passed
    # explicitly here for a deterministic assertion -- the real default
    # (omitted) is now live-detected from whatever machine runs it, see
    # test_run_snakemake_defaults_concurrent_slots_to_a_live_detection below.
    assert "getorganelle_slots=8" in command
    # Regression: the default 5s --latency-wait was too tight once
    # --use-conda's real per-job activation overhead was ever actually
    # exercised, raising a real MissingOutputException.
    assert command[command.index("--latency-wait") + 1] == "30"
    # Regression: a real interruption leaves Snakemake unable to trust
    # an in-flight job's own output (a real IncompleteFilesException)
    # without this flag -- an unattended launcher must resume on its
    # own, not require a manual --cleanup-metadata/--rerun-incomplete step.
    assert "--rerun-incomplete" in command

def test_run_snakemake_defaults_concurrent_slots_to_a_live_detection(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Omitting concurrent_slots must never fall back to a number
    hardcoded for one specific machine -- it has to come from a real,
    live detection of whatever machine is actually running this."""
    captured_commands: list[list[str]] = []
    monkeypatch.setattr(subprocess, "Popen", lambda command, **kwargs: (captured_commands.append(command), _FakeProcess())[1])
    monkeypatch.setattr(pipeline_launcher, "recommended_concurrent_slots", lambda: 42)

    run_snakemake(runtime_config = tmp_path / "runtime_config.yaml", cores = 2, use_conda = True,
                  logger = logging.getLogger("test"))

    assert "getorganelle_slots=42" in captured_commands[0]

def test_rerun_triggers_flag_is_parsed_and_defaults_to_none(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(sys, "argv", ["pipeline_launcher.py", "--config", "config.yaml"])
    assert parse_args().rerun_triggers is None
    monkeypatch.setattr(sys, "argv", ["pipeline_launcher.py", "--config", "config.yaml", "--rerun-triggers", "mtime"])
    assert parse_args().rerun_triggers == ["mtime"]

def test_run_snakemake_passes_rerun_triggers_through_only_when_given(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Resuming after output directories were moved needs Snakemake's own
    --rerun-triggers mtime, or every finished job reruns because its
    recorded input paths changed; omitted, Snakemake's default stands."""
    captured_commands: list[list[str]] = []
    monkeypatch.setattr(subprocess, "Popen", lambda command, **kwargs: (captured_commands.append(command), _FakeProcess())[1])

    run_snakemake(runtime_config = tmp_path / "runtime_config.yaml", cores = 2, use_conda = True,
                  logger = logging.getLogger("test"))
    assert "--rerun-triggers" not in captured_commands[0]

    run_snakemake(runtime_config = tmp_path / "runtime_config.yaml", cores = 2, use_conda = True,
                  logger = logging.getLogger("test"), rerun_triggers = ["mtime"])
    command = captured_commands[1]
    assert command[command.index("--rerun-triggers") + 1] == "mtime"

class _FakeDiscoveryJob:
    """Only what discover_fastq_pairs()/write_*_manifest() touch: a logger."""
    def __init__(self) -> None:
        self.job_logger = logging.getLogger("test_discovery")

def _touch_fastqs(directory: Path, names: list[str]) -> None:
    directory.mkdir(parents = True, exist_ok = True)
    for name in names:
        (directory / name).write_bytes(b"")

def test_discover_fastq_pairs_pairs_illumina_and_sra_style_names(tmp_path: Path) -> None:
    """Regression: discover_fastq_pairs() crashed on its first line
    (`fastq_files = list[Path] = []`, a chained assignment into the list
    type, a real TypeError) so a manifest-less run could never start. Now
    it pairs Illumina "_R1_001"/"_R2_001", bare "_R1"/"_R2" and SRA
    "_1"/"_2" names, keying each pair on everything before the token."""
    _touch_fastqs(tmp_path, [
        "54_FMNH_124727_CGTGTAACTT-TTGTCACGTT_L001_R1_001.fastq.gz",
        "54_FMNH_124727_CGTGTAACTT-TTGTCACGTT_L001_R2_001.fastq.gz",
        "SRR000001_1.fastq.gz", "SRR000001_2.fastq.gz",
        "tiny_R1.fq", "tiny_R2.fq",
    ])
    records = discover_fastq_pairs(tmp_path, _FakeDiscoveryJob())
    by_id = {r["sample_id"]: r for r in records}
    assert sorted(by_id) == ["54_FMNH_124727_CGTGTAACTT-TTGTCACGTT_L001", "SRR000001", "tiny"]
    assert by_id["SRR000001"]["r1"].endswith("SRR000001_1.fastq.gz")
    assert by_id["SRR000001"]["r2"].endswith("SRR000001_2.fastq.gz")
    assert by_id["tiny"]["r1"].endswith("tiny_R1.fq")
    assert all(r["source"] == "discovered_fastq" for r in records)

def test_discover_fastq_pairs_does_not_mistake_an_inner_underscore_digit_for_the_mate_token(tmp_path: Path) -> None:
    """The old substring check split "FMNH_124727..." at the first "_1"; the
    mate token must sit right before the extension (or "_001")."""
    _touch_fastqs(tmp_path, ["FMNH_124727_R1.fastq.gz", "FMNH_124727_R2.fastq.gz"])
    records = discover_fastq_pairs(tmp_path, _FakeDiscoveryJob())
    assert [r["sample_id"] for r in records] == ["FMNH_124727"]

def test_discover_fastq_pairs_rejects_a_missing_mate(tmp_path: Path) -> None:
    _touch_fastqs(tmp_path, ["a_R1.fastq.gz", "a_R2.fastq.gz", "b_R1.fastq.gz"])
    with pytest.raises(SystemExit):
        discover_fastq_pairs(tmp_path, _FakeDiscoveryJob())

def test_discover_fastq_pairs_rejects_an_empty_or_missing_directory(tmp_path: Path) -> None:
    with pytest.raises(SystemExit):
        discover_fastq_pairs(tmp_path / "missing", _FakeDiscoveryJob())
    _touch_fastqs(tmp_path / "unpaired", ["single.fastq.gz"])
    with pytest.raises(SystemExit):
        discover_fastq_pairs(tmp_path / "unpaired", _FakeDiscoveryJob())

def test_validated_manifest_always_carries_an_organism_column(tmp_path: Path) -> None:
    """The rule layer reads SAMPLE_TABLE["organism"] directly, so the
    runtime manifest must always have the column: derived from genus/
    species when no explicit organism was given, blank when nothing is
    known (a discovered sample)."""
    fastq = tmp_path / "x_R1.fastq.gz"; fastq.write_bytes(b"")
    samples = [
        Sample(sample_id = "named", r1 = fastq, r2 = fastq, metadata = {"organism": "Danio rerio"}),
        Sample(sample_id = "derived", r1 = fastq, r2 = fastq, genus = "Danio", species = "rerio", metadata = {}),
        Sample(sample_id = "nameless", r1 = fastq, r2 = fastq, metadata = {"source": "discovered_fastq"}),
    ]
    manifest = tmp_path / "runtime_manifest.tsv"
    write_validated_manifest(manifest, samples, _FakeDiscoveryJob())

    import pandas as pd
    table = pd.read_csv(manifest, sep = "\t", dtype = str, keep_default_na = False).set_index("sample_id")
    assert table.loc["named", "organism"] == "Danio rerio"
    assert table.loc["derived", "organism"] == "Danio rerio"
    assert table.loc["nameless", "organism"] == ""

def test_prepare_runtime_manifest_runs_end_to_end_without_a_manifest(tmp_path: Path) -> None:
    """The whole manifest-less path: discover from input_directory, write,
    re-parse through the same validator a user manifest goes through, and
    land a runtime manifest the Snakefile can load."""
    input_dir = tmp_path / "raw"
    _touch_fastqs(input_dir, ["s1_R1_001.fastq.gz", "s1_R2_001.fastq.gz"])
    job = _FakeDiscoveryJob(); job.job_dir = tmp_path / "job"; job.job_dir.mkdir()
    config = {"job_configuration": {"input_directory": str(input_dir)}}

    runtime_manifest, samples = prepare_runtime_manifest(config, job)
    assert runtime_manifest.is_file()
    assert [s.sample_id for s in samples] == ["s1"]
    assert samples[0].scientific_name is None
    header = runtime_manifest.read_text().splitlines()[0].split("\t")
    assert {"sample_id", "r1", "r2", "organism"} <= set(header)

class _FakeJob:
    """A minimal stand-in exposing only what _rerender_run_report_with_final_status() actually touches -- never a real PipelineJob, which creates real directories and log files as a side effect of construction."""
    def __init__(self, job_dir: Path) -> None:
        self.job_dir = job_dir
        self.job_id = "job1"
        self.job_logger = logging.getLogger("test_rerender")

class _FakeSample:
    def __init__(self, sample_id: str) -> None:
        self.sample_id = sample_id

def test_rerender_run_report_is_a_noop_when_reporting_never_ran(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Regression guard: a run with reporting disabled (or a run that
    failed before ever reaching it) must not attempt to re-invoke a
    report-rendering script against a report that was never rendered
    in the first place."""
    calls: list = []
    monkeypatch.setattr(subprocess, "run", lambda *a, **k: calls.append((a, k)))

    job = _FakeJob(tmp_path / "job")
    _rerender_run_report_with_final_status({}, job, tmp_path / "manifest.tsv", [_FakeSample("s1")])

    assert calls == []

def test_rerender_run_report_invokes_the_real_exec_script_with_the_right_arguments(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Regression: the Snakemake rule that first renders the run-level
    report runs mid-DAG, before this launcher's own mark_completed()/
    mark_failed() call -- so its own first pass can only show the
    job's mid-run status, never the real terminal one. This second
    pass, invoked once current_metadata.json carries the real final
    status, must call the same real exec script with the same real
    arguments the Snakemake rule itself would use."""
    job_dir = tmp_path / "job"
    report_dir = job_dir / "13-report_pdf"
    report_dir.mkdir(parents = True)
    (report_dir / "run_report.done").touch()

    captured: dict = {}
    def _fake_run(command, **kwargs):
        captured["command"] = command
        return subprocess.CompletedProcess(args = command, returncode = 0, stdout = "", stderr = "")
    monkeypatch.setattr(subprocess, "run", _fake_run)

    config = {
        "stage_options": {
            "report_export": {"pandoc": {"bin": "pandoc", "engine": "/usr/local/bin/typst"}},
            "phylogeny": {"iqtree": {"model": "MFP"}},
        },
    }
    job = _FakeJob(job_dir)
    samples = [_FakeSample("sample_a"), _FakeSample("sample_b")]
    runtime_manifest = tmp_path / "runtime_manifest.tsv"

    _rerender_run_report_with_final_status(config, job, runtime_manifest, samples)

    command = captured["command"]
    assert command[0] == sys.executable
    assert "run_run_report.py" in command[1]
    assert command[command.index("--job-id") + 1] == "job1"
    assert command.count("--sample-id") == 2
    assert "sample_a" in command and "sample_b" in command
    assert command[command.index("--runtime-manifest") + 1] == str(runtime_manifest)
    assert command[command.index("--job-dir") + 1] == str(job_dir)
    assert command[command.index("--pandoc-bin") + 1] == "pandoc"
    assert command[command.index("--pdf-engine") + 1] == "/usr/local/bin/typst"
    assert command[command.index("--phylogeny-model") + 1] == "MFP"
    assert command[command.index("--output-dir") + 1] == str(report_dir)

def test_rerender_run_report_omits_phylogeny_model_when_not_configured(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    job_dir = tmp_path / "job"
    report_dir = job_dir / "13-report_pdf"
    report_dir.mkdir(parents = True)
    (report_dir / "run_report.done").touch()

    captured: dict = {}
    def _fake_run(command, **kwargs):
        captured["command"] = command
        return subprocess.CompletedProcess(args = command, returncode = 0, stdout = "", stderr = "")
    monkeypatch.setattr(subprocess, "run", _fake_run)

    job = _FakeJob(job_dir)
    _rerender_run_report_with_final_status({}, job, tmp_path / "manifest.tsv", [_FakeSample("s1")])

    assert "--phylogeny-model" not in captured["command"]
