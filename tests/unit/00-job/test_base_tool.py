"""
test_base_tool.py

Unit tests for the BaseTool class.
"""

# Imports
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
