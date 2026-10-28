"""
base_tool.py

Defines the base class for a pipeline tool API wrapper. Contains functionality for validating inputs and outputs and 
building/executing commands.
"""

# Imports
from __future__ import annotations
from abc import ABC, abstractmethod
from datetime import datetime
from logging import Logger
from pathlib import Path
import shlex
import  subprocess
import time
import resource
import tempfile
import json
from mitopipeline.models.command_result.command_result import CommandResult
from mitopipeline.logging.logger_factory import TOOL_OUTPUT_EVENT

class BaseTool(ABC):
    """
    Base class for external command-line tool wrappers.

    Attributes:
        tool_name (str): The name of the tool.
        working_dir (Path): The working directory for the tool.
        logger (Logger): The logger to use for the step of the Pipeline.
    """

    def __init__(self,
                 tool_name: str,
                 working_dir: Path,
                 logger: Logger,
                 metrics_dir: Path | None = None) -> None:
        """
        Instantiation method for the BaseTool class.

        Args:
            tool_name (str): The name of the tool.
            working_dir (Path): The working directory for the tool.
            logger (Logger): The logger to use for the step of the Pipeline.
            metrics_dir (Path | None): Where to persist this run's CommandResult as JSON (organized, report-ready per-step compute metrics). run() already measures CPU/memory/runtime for every tool call via GNU time. None (the default) skips persistence entirely -- existing callers that don't pass it keep working exactly as before.

        Returns:
            None
        """
        self.tool_name = tool_name
        self.working_dir = working_dir
        self.logger = logger
        self.metrics_dir = Path(metrics_dir) if metrics_dir is not None else None

    def version_command(self) -> list[str] | None:
        """
        The real command to query this tool's own version (e.g.
        ["fastp", "--version"]). None (the base default) means this
        tool's own version is never checked -- a subclass overrides
        this only when it knows the real, correct flag for its own
        wrapped tool; never guessed.

        Args:
            None

        Returns:
            list[str] | None: The version-query command, or None to skip version detection entirely.
        """
        return None

    def _detect_tool_version(self) -> str | None:
        """
        Runs this tool's own real version_command(), if provided, and
        returns the first real non-empty line of its output (stdout
        preferred; some tools print their version to stderr instead).
        Never raises -- a missing binary, an unexpected flag, or any
        other real failure here must never block the actual pipeline
        command this exists alongside; returns None in every such case.

        Args:
            None

        Returns:
            str | None: The real detected version string, or None if version_command() is None or the real check failed.
        """
        command = self.version_command()
        if command is None:
            return None
        try:
            result = subprocess.run(command, capture_output = True, text = True, timeout = 10)
        except (OSError, subprocess.TimeoutExpired):
            return None
        for stream in (result.stdout, result.stderr):
            for line in stream.splitlines():
                if line.strip():
                    return line.strip()
        return None

    @abstractmethod
    def validate_inputs(self) -> None:
        """
        validate command inputs.

        Args:
            None

        Returns:
            None
        """

    @abstractmethod
    def build_command(self) -> list[str]:
        """
        Builds the external command to be executed.
        
        Args:
            None

        Returns:
            list[str]: The external command to be executed.
        """

    @abstractmethod
    def validate_outputs(self) -> None:
        """
        Validate command outputs.

        Args:
            None

        Returns:
            None
        """

    def run(self) -> CommandResult:
        """
        Runs the external command.

        Args:
            None

        Returns:
            CommandResult: The result of the command execution.
        """
        # Obtaining the context (sample) for the logger.
        context = self._log_context()

        # Logging and executing validation portion. 
        self.logger.info(f"{context} Validating inputs.")
        self.validate_inputs()
        self.logger.info(f"{context} Suecessfully validated inputs.")

        # Loggina dnd executing build command portion. 
        self.logger.info(f"{context} Building command.")
        command = self.build_command()
        self.logger.info(f"{context} Successfully built command.")
        self.logger.debug(f"{context} Command: {shlex.join(command)}")

        # Creating a tempfile for GNU time statistics.
        with tempfile.NamedTemporaryFile(mode = "w", suffix = ".json", delete = False) as handle:
            metrics_path = Path(handle.name)
            timed_command = [
                "/usr/bin/time",
                "-f",
                (
                    '{"user_cpu_seconds": %U, '
                    '"system_cpu_seconds": %S, '
                    '"peak_memory_mb": %M}'
                ),
                "-o",
                str(metrics_path),
                "--",
                *command,
            ]

        # Logging and executing run portion.
        start_time = time.perf_counter()
        started_at = datetime.now()
        self.logger.info(f"{context} Running command.")
        self.logger.debug(f"{context} Command: {shlex.join(command)}")
        self.logger.debug(f"{context} Working directory: {self.working_dir}.")

        try:
            # Running the command. 
            completed = subprocess.run(timed_command, cwd = self.working_dir, capture_output = True, text = True)
            runtime_seconds = time.perf_counter() - start_time
            ended_at = datetime.now()

            # Reading the execution resource statistics. On a non-zero exit,
            # GNU time prepends a human-readable notice line (e.g. "Command
            # exited with non-zero status 1") before the -f JSON line, so
            # the JSON line is parsed specifically rather than the whole file.
            with metrics_path.open("r", encoding = "utf-8") as handle:
                metrics_lines = [line for line in handle.readlines() if line.strip()]
            metrics = json.loads(metrics_lines[-1]) if metrics_lines else {}
            user_cpu_seconds = metrics["user_cpu_seconds"]
            system_cpu_seconds = metrics["system_cpu_seconds"]
            peak_memory_mb = metrics["peak_memory_mb"] / 1024

        finally:
            # Removing the tempfile.
            metrics_path.unlink(missing_ok = True)

        # Obtaining the CommandResult.
        command_result = CommandResult(command = command,
                                       return_code = completed.returncode,
                                       stdout = completed.stdout,
                                       stderr = completed.stderr,
                                       runtime_seconds = runtime_seconds,
                                       success = completed.returncode == 0,
                                       tool_name = self.tool_name,
                                       started_at = started_at,
                                       ended_at = ended_at,
                                       peak_memory_mb = peak_memory_mb,
                                       user_cpu_seconds = user_cpu_seconds,
                                       system_cpu_seconds = system_cpu_seconds,
                                       tool_version = self._detect_tool_version())

        # Persisting compute metrics BEFORE the success/failure branch
        # below, so a failed run's metrics are captured too, not just a
        # successful one's.
        self._write_metrics(command_result)

        # Logic for dealing with different command result success states.
        self.logger.info(f"{context} Execution finished with return code {command_result.return_code} after {command_result.runtime_seconds} seconds.")
        if command_result.stdout:
            self.logger.debug(f"{context} stdout:\n{command_result.stdout.rstrip()}",
                              extra = {"event_type": TOOL_OUTPUT_EVENT})
        if not command_result.success:
            self.logger.error(f"{context} Execution failed with return code {command_result.return_code} after {command_result.runtime_seconds} seconds.")
            return command_result

        # Postprocessing and validating outputs.
        self.logger.info(f"{context} Postprocessing and validating outputs.")
        self.postprocess_outputs()
        self.validate_outputs()
        self.logger.info(f"{context} Successfully postprocessed and validated outputs.")

        return command_result

    def postprocess_outputs(self) -> None:
        """
        This is an optional method which normalizes/postprocesses outpputs after successful execution.

        Args:
            None

        Returns:
            None
        """

    def _write_metrics(self, command_result: CommandResult) -> None:
        """
        Persists command_result as JSON under metrics_dir, one file per
        tool invocation. A write failure is logged, never raised -- a
        missing metrics file should not be able to break an otherwise
        successful (or already-failed) pipeline step.

        Filenames include a sample label (when reachable, via the same
        _log_context() logic) and a start timestamp, so repeated calls
        from the same tool/sample across different pipeline attempts
        never collide and overwrite each other.

        Args:
            command_result (CommandResult): The result to persist.

        Returns:
            None
        """
        if self.metrics_dir is None:
            return

        sample_id = getattr(self, "sample_id", None) or getattr(getattr(self, "sample", None), "sample_id", None)
        label = f"{self.tool_name}.{sample_id}" if sample_id else self.tool_name
        timestamp = command_result.started_at.strftime("%Y%m%dT%H%M%S%f")
        metrics_path = self.metrics_dir / f"{label}.{timestamp}.metrics.json"

        try:
            self.metrics_dir.mkdir(parents = True, exist_ok = True)
            with metrics_path.open("w", encoding = "utf-8") as handle:
                json.dump(command_result.to_dict(), handle, indent = 4)
        except OSError as error:
            self.logger.warning(f"{self._log_context()} Could not write compute metrics to {metrics_path}: {error}.")

    def _log_context(self) -> str:
        """
        Adds a sample label for log messages.

        Returns:
            str: The sample label.
        """

        # Attempting to access a direct sample ID.
        sample_id = getattr(self, "sample_id", None)

        # Attempting to access the sample ID from a Sample object.
        if sample_id is None:
            sample = getattr(self, "sample", None)
            sample_id = getattr(sample, "sample_id", None)

        # Returning the formatted sample label.
        if sample_id:
            return f"{{{sample_id}}}"

        return ""
