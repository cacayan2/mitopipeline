"""
test_base_tool.py

Unit tests for the BaseTool class.
"""

# Imports
import json
import logging
from pathlib import Path
import pytest
from mitopipeline.models.base_tool.base_tool import BaseTool
from mitopipeline.models.command_result.command_result import CommandResult

class _EchoTool(BaseTool):
    """
    A minimal concrete BaseTool that shells out to `echo`, used to
    exercise BaseTool.run() without depending on any real pipeline tool.
    """
    def __init__(self, working_dir: Path, logger: logging.Logger, message: str = "hello") -> None:
        super().__init__(tool_name = "echo", working_dir = working_dir, logger = logger)
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
