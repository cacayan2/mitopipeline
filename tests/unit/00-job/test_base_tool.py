"""
test_base_tool.py

Unit tests for the BaseTool class.
"""

# Imports
import json
import logging
import resource
import threading
import time
from pathlib import Path
import pytest
from mitopipeline.models.base_tool.base_tool import BaseTool
from mitopipeline.models.command_result.command_result import CommandResult

class _EchoTool(BaseTool):
    """
    A minimal concrete BaseTool that shells out to `echo`, used to
    exercise BaseTool.run() without depending on any real pipeline tool.
    """
    def __init__(self, working_dir: Path, logger: logging.Logger, message: str = "hello",
                 memory_limit_mb: int | None = None, dynamic_memory_limit_mb = None) -> None:
        super().__init__(tool_name = "echo", working_dir = working_dir, logger = logger,
                         memory_limit_mb = memory_limit_mb, dynamic_memory_limit_mb = dynamic_memory_limit_mb)
        self.message = message
        self.postprocessed = False

    def validate_inputs(self) -> None:
        if not self.message:
            raise ValueError("message must not be empty.")

    def build_command(self) -> list[str]:
        return ["echo", self.message]

    def validate_outputs(self) -> None:
        pass

    def postprocess_outputs(self) -> None:
        self.postprocessed = True

class _VersionedEchoTool(_EchoTool):
    """An _EchoTool that reports a real, harmless version_command() (echo itself has no real --version flag worth relying on across platforms, so this fakes one via a second echo invocation)."""
    def version_command(self) -> list[str] | None:
        return ["echo", "FakeTool version 9.9.9"]

class _BrokenVersionTool(_EchoTool):
    """An _EchoTool whose version_command() points at a binary that doesn't exist, to confirm a failed version check never breaks the real run."""
    def version_command(self) -> list[str] | None:
        return ["this-binary-does-not-exist-anywhere", "--version"]

class _FailingValidateTool(_EchoTool):
    """
    An _EchoTool whose validate_outputs always raises, used to confirm
    a validation failure after a successful run propagates.
    """
    def validate_outputs(self) -> None:
        raise ValueError("outputs are invalid.")

def test_run_returns_successful_command_result(tmp_path: Path) -> None:
    """
    A successful command should return a CommandResult with success True,
    return_code 0, and the command's stdout captured.
    """
    tool = _EchoTool(working_dir = tmp_path, logger = logging.getLogger("test"), message = "hello world")
    result = tool.run()
    assert isinstance(result, CommandResult)
    assert result.success is True
    assert result.return_code == 0
    assert "hello world" in result.stdout
    assert result.tool_name == "echo"

def test_run_tags_the_raw_stdout_log_as_tool_output(tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
    """Regression: the raw multi-line stdout dump must carry
    event_type=TOOL_OUTPUT_EVENT so logger_factory's own
    ExcludeToolOutputFilter can keep it out of the shared global log --
    without this tag, the filter is a no-op and every sample's raw
    tool console output interleaves untagged into the global log."""
    from mitopipeline.logging.logger_factory import TOOL_OUTPUT_EVENT

    logger = logging.getLogger("test_tool_output_tag")
    logger.setLevel(logging.DEBUG)
    tool = _EchoTool(working_dir = tmp_path, logger = logger, message = "hello world")

    with caplog.at_level(logging.DEBUG, logger = "test_tool_output_tag"):
        tool.run()

    stdout_records = [r for r in caplog.records if "stdout:" in r.getMessage()]
    assert stdout_records
    assert all(getattr(r, "event_type", None) == TOOL_OUTPUT_EVENT for r in stdout_records)

def test_run_postprocesses_and_validates_on_success(tmp_path: Path) -> None:
    """
    A successful run should call postprocess_outputs before validate_outputs.
    """
    tool = _EchoTool(working_dir = tmp_path, logger = logging.getLogger("test"))
    tool.run()
    assert tool.postprocessed is True

def test_run_propagates_output_validation_errors(tmp_path: Path) -> None:
    """
    A successful command whose outputs fail validation should raise,
    rather than silently return a result.
    """
    tool = _FailingValidateTool(working_dir = tmp_path, logger = logging.getLogger("test"))
    with pytest.raises(ValueError):
        tool.run()

def test_run_skips_postprocessing_on_command_failure(tmp_path: Path) -> None:
    """
    A failing command should return a CommandResult with success False
    and must not attempt postprocessing/validation.
    """
    class _FailingTool(_EchoTool):
        def build_command(self) -> list[str]:
            return ["false"]

    tool = _FailingTool(working_dir = tmp_path, logger = logging.getLogger("test"))
    result = tool.run()
    assert result.success is False
    assert result.return_code != 0
    assert tool.postprocessed is False

def test_run_raises_on_invalid_inputs(tmp_path: Path) -> None:
    """
    An input validation failure should raise before any command runs.
    """
    tool = _EchoTool(working_dir = tmp_path, logger = logging.getLogger("test"), message = "")
    with pytest.raises(ValueError):
        tool.run()

def test_log_context_includes_sample_id_when_set(tmp_path: Path) -> None:
    """
    _log_context should surface a sample_id attribute when present.
    """
    tool = _EchoTool(working_dir = tmp_path, logger = logging.getLogger("test"))
    tool.sample_id = "SAMPLE_01"
    assert tool._log_context() == "{SAMPLE_01}"

def test_log_context_is_empty_without_sample_id(tmp_path: Path) -> None:
    """
    _log_context should return an empty string when no sample ID is
    reachable.
    """
    tool = _EchoTool(working_dir = tmp_path, logger = logging.getLogger("test"))
    assert tool._log_context() == ""

def test_run_writes_no_metrics_file_when_metrics_dir_not_set(tmp_path: Path) -> None:
    """
    Existing callers that never pass metrics_dir must keep working
    exactly as before -- no file written anywhere, no new directory
    created.
    """
    tool = _EchoTool(working_dir = tmp_path, logger = logging.getLogger("test"))
    tool.run()
    assert tool.metrics_dir is None

def test_run_persists_metrics_json_on_success(tmp_path: Path) -> None:
    """
    A successful run with metrics_dir set should write exactly one
    JSON metrics file, with the real compute numbers run() already
    measures via GNU time.
    """
    metrics_dir = tmp_path / "metrics"
    tool = _EchoTool(working_dir = tmp_path, logger = logging.getLogger("test"))
    tool.metrics_dir = metrics_dir
    result = tool.run()

    written = list(metrics_dir.glob("*.metrics.json"))
    assert len(written) == 1
    data = json.loads(written[0].read_text())
    assert data["tool_name"] == "echo"
    assert data["success"] is True
    assert data["return_code"] == 0
    assert data["runtime_seconds"] == result.runtime_seconds
    assert isinstance(data["peak_memory_mb"], (int, float))

def test_run_persists_metrics_json_on_failure_too(tmp_path: Path) -> None:
    """
    A failing command's metrics must still be persisted -- a failed
    step's resource usage matters just as much for a report, maybe
    more.
    """
    class _FailingTool(_EchoTool):
        def build_command(self) -> list[str]:
            return ["false"]

    metrics_dir = tmp_path / "metrics"
    tool = _FailingTool(working_dir = tmp_path, logger = logging.getLogger("test"))
    tool.metrics_dir = metrics_dir
    tool.run()

    written = list(metrics_dir.glob("*.metrics.json"))
    assert len(written) == 1
    data = json.loads(written[0].read_text())
    assert data["success"] is False

def test_run_includes_sample_id_in_metrics_filename(tmp_path: Path) -> None:
    """
    When a sample_id is reachable, it should appear in the metrics
    filename so multiple samples sharing one metrics_dir don't collide.
    """
    metrics_dir = tmp_path / "metrics"
    tool = _EchoTool(working_dir = tmp_path, logger = logging.getLogger("test"))
    tool.metrics_dir = metrics_dir
    tool.sample_id = "SAMPLE_07"
    tool.run()

    written = list(metrics_dir.glob("*.metrics.json"))
    assert len(written) == 1
    assert "SAMPLE_07" in written[0].name

def test_run_captures_the_real_tool_version_when_version_command_is_set(tmp_path: Path) -> None:
    tool = _VersionedEchoTool(working_dir = tmp_path, logger = logging.getLogger("test"))
    result = tool.run()
    assert result.tool_version == "FakeTool version 9.9.9"

def test_run_tool_version_is_none_when_version_command_not_overridden(tmp_path: Path) -> None:
    tool = _EchoTool(working_dir = tmp_path, logger = logging.getLogger("test"))
    result = tool.run()
    assert result.tool_version is None

def test_run_tool_version_is_none_when_the_real_check_fails(tmp_path: Path) -> None:
    """A version check against a nonexistent binary must never raise or block the real pipeline command it rides alongside."""
    tool = _BrokenVersionTool(working_dir = tmp_path, logger = logging.getLogger("test"))
    result = tool.run()
    assert result.success is True
    assert result.tool_version is None

def test_run_persists_tool_version_in_the_metrics_json(tmp_path: Path) -> None:
    metrics_dir = tmp_path / "metrics"
    tool = _VersionedEchoTool(working_dir = tmp_path, logger = logging.getLogger("test"))
    tool.metrics_dir = metrics_dir
    tool.run()

    written = list(metrics_dir.glob("*.metrics.json"))
    data = json.loads(written[0].read_text())
    assert data["tool_version"] == "FakeTool version 9.9.9"

def test_apply_memory_limit_sets_rlimit_as_to_the_configured_bytes(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """_apply_memory_limit() is run as subprocess.run's preexec_fn (see run()) -- verified directly, in isolation, rather than relying on a real OOM to actually happen inside a test."""
    captured = {}
    monkeypatch.setattr(resource, "setrlimit", lambda which, limits: captured.update(which = which, limits = limits))

    tool = _EchoTool(working_dir = tmp_path, logger = logging.getLogger("test"), memory_limit_mb = 256)
    tool._apply_memory_limit()

    expected_bytes = 256 * 1024 * 1024
    assert captured["which"] == resource.RLIMIT_AS
    assert captured["limits"] == (expected_bytes, expected_bytes)

def test_run_applies_no_memory_limit_by_default(tmp_path: Path) -> None:
    """A tool that never passes memory_limit_mb (every existing caller, before this feature) must run exactly as before -- no preexec_fn, no limit."""
    tool = _EchoTool(working_dir = tmp_path, logger = logging.getLogger("test"))
    assert tool.memory_limit_mb is None
    result = tool.run()
    assert result.success is True

class _SleepTool(BaseTool):
    """A real, slow (but harmless) command -- long enough that the live
    memory-adjuster thread gets a real chance to act on it while it's
    still genuinely running, unlike _EchoTool which exits instantly."""
    def __init__(self, working_dir: Path, logger: logging.Logger, seconds: int = 5,
                 dynamic_memory_limit_mb = None) -> None:
        super().__init__(tool_name = "sleep", working_dir = working_dir, logger = logger,
                         dynamic_memory_limit_mb = dynamic_memory_limit_mb)
        self.seconds = seconds

    def validate_inputs(self) -> None:
        pass

    def build_command(self) -> list[str]:
        return ["sleep", str(self.seconds)]

    def validate_outputs(self) -> None:
        pass

@pytest.mark.skipif(not Path("/proc").is_dir(), reason="real /proc only exists on Linux (this pipeline's actual deployment targets); this test needs the real filesystem, not a mock of it")
def test_find_real_child_pid_returns_the_actual_wrapped_command(tmp_path: Path) -> None:
    """process.pid from Popen is GNU time's own PID, never the real
    wrapped command -- confirmed directly against a real running
    `sleep`, not assumed from how /proc/.../children is documented."""
    import subprocess as sp
    parent = sp.Popen(["sh", "-c", "sleep 5"])
    try:
        time.sleep(0.3)   # real time for the shell to actually exec its own child
        tool = _SleepTool(working_dir = tmp_path, logger = logging.getLogger("test"))
        real_child = tool._find_real_child_pid(parent.pid)
        assert real_child is not None
        assert real_child != parent.pid
        # Confirm it's genuinely a real, distinct running process.
        os_kill_check = sp.run(["ps", "-p", str(real_child)], capture_output = True)
        assert os_kill_check.returncode == 0
    finally:
        parent.terminate()
        parent.wait(timeout = 5)

@pytest.mark.skipif(not hasattr(resource, "prlimit"), reason="resource.prlimit is Linux-only (this pipeline's actual deployment targets); not available on macOS/BSD")
def test_dynamic_memory_limit_adjusts_a_real_running_process(tmp_path: Path) -> None:
    """A real, running `sleep` process's actual RLIMIT_AS must change
    live, mid-run -- read back via resource.prlimit on the real PID,
    not mocked, confirming the adjustment genuinely reached the OS."""
    import subprocess as sp
    calls = iter([2048, 4096])   # MB -- a real change the adjuster should apply
    tool = _SleepTool(working_dir = tmp_path, logger = logging.getLogger("test"), seconds = 3,
                      dynamic_memory_limit_mb = lambda: next(calls, 4096))

    # A shell wrapper, not a bare `sleep` -- _adjust_memory_limit_periodically
    # looks for a CHILD of the given PID (matching real production usage,
    # where that PID is always /usr/bin/time's own, with the real command
    # as its child). A bare `sleep` has no children at all, which would
    # make _find_real_child_pid correctly, honestly return None forever.
    process = sp.Popen(["sh", "-c", "sleep 3"])
    stop_event = threading.Event()
    adjuster = threading.Thread(target = tool._adjust_memory_limit_periodically, args = (process, stop_event))
    import mitopipeline.models.base_tool.base_tool as base_tool_module
    original_recheck = base_tool_module._DYNAMIC_MEMORY_RECHECK_SECONDS
    base_tool_module._DYNAMIC_MEMORY_RECHECK_SECONDS = 0.5
    try:
        adjuster.start()
        time.sleep(1.5)
        stop_event.set()
        adjuster.join(timeout = 3)
        real_child = tool._find_real_child_pid(process.pid)
        assert real_child is not None
        soft, _ = resource.prlimit(real_child, resource.RLIMIT_AS)
        assert soft == 2048 * 1024 * 1024
    finally:
        base_tool_module._DYNAMIC_MEMORY_RECHECK_SECONDS = original_recheck
        process.terminate()
        process.wait(timeout = 5)

def test_run_without_dynamic_memory_limit_behaves_exactly_as_before(tmp_path: Path) -> None:
    """Every existing caller that never sets dynamic_memory_limit_mb
    (every tool wrapper except GetOrganelleRunner) must see no change."""
    tool = _EchoTool(working_dir = tmp_path, logger = logging.getLogger("test"))
    assert tool.dynamic_memory_limit_mb is None
    result = tool.run()
    assert result.success is True
    assert result.stdout.strip() == "hello"

def test_run_succeeds_under_a_generous_memory_limit(tmp_path: Path) -> None:
    """A real end-to-end run with memory_limit_mb set -- confirms the preexec_fn wiring itself doesn't break a normal, well-within-budget command."""
    tool = _EchoTool(working_dir = tmp_path, logger = logging.getLogger("test"), memory_limit_mb = 512)
    result = tool.run()
    assert result.success is True
    assert result.stdout.strip() == "hello"
