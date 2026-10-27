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
from mitopipeline.launcher.pipeline_launcher import _configured_job_id, _job_configuration, load_config, parse_args, run_snakemake

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
        logger = logging.getLogger("test"),
    )

    command = captured_commands[0]
    assert command[0] == str(Path(sys.executable).parent / "snakemake")
    assert "--conda-not-block-search-path-envvars" in command
    # Regression: "entrez = 1" (with spaces) is not a valid Snakemake
    # --resources identifier and raised a real WorkflowError.
    assert command[command.index("--resources") + 1] == "entrez=1"
    # Regression: the default 5s --latency-wait was too tight once
    # --use-conda's real per-job activation overhead was ever actually
    # exercised, raising a real MissingOutputException.
    assert command[command.index("--latency-wait") + 1] == "30"
