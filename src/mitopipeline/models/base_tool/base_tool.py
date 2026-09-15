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
                 logger: Logger) -> None:
        """
        Instantiation method for the BaseTool class.

        Args:
            tool_name (str): The name of the tool.
            working_dir (Path): The working directory for the tool.
            logger (Logger): The logger to use for the step of the Pipeline.
        
        Returns:
            None
        """
        self.tool_name = tool_name
        self.working_dir = working_dir
        self.logger = logger

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

    @abstractmethod
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
            completed = subprocess.run(command, cwd = self.working_dir, capture_output = True, text = True)
            runtime_seconds = time.perf_counter() - start_time
            ended_at = datetime.now()

            # Reading the execution resource statistics.
            with metrics_path.open("r", encoding = "utf-8") as handle:
                metrics = json.load(handle)
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
                                       system_cpu_seconds = system_cpu_seconds)

        # Logic for dealing with different command result success states.
        self.logger.info(f"{context} Execution finished with return code {command_result.return_code} after {command_result.runtime_seconds} seconds.")
        if command_result.stdout: self.logger.debug(f"{context} stdout:\n{command_result.stdout.rstrip()}")
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
