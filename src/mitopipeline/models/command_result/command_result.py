"""
command_result.py

This module contains the CommandResult class, specifying the behavior of command execution.
"""

# Imports
from datetime import datetime

class CommandResult:
    """
    A class specifying the behavior associated with running a command.

    Attributes:
        command (str): The command that was executed.
        return_code (int): The return code from the command.
        stdout (str): The standard output from the command.
        stderr (str): The standard error from the command.
        runtime_seconds (float): The runtime of the command in seconds.
        success (bool): Whether the command was successful.
        tool_name (str): The name of the tool that was executed.
        started_at (datetime): The time at which the command was started.
        ended_at (datetime): The time at which the command was completed.
        peak_memory_mb (float): The peak memory usage of the command in megabytes.
        cpu_time_seconds (float): The CPU time of the command in seconds.
        user_cpu_seconds (float): The user CPU time of the command in seconds.
        system_cpu_seconds (float): The system CPU time of the command in seconds.
    """
    def __init__(self, command: str,
                 return_code: int,
                 stdout: str,
                 stderr: str,
                 runtime_seconds: float,
                 success: bool,
                 tool_name: str,
                 started_at: datetime,
                 ended_at: datetime,
                 peak_memory_mb: float,
                 user_cpu_seconds: float,
                 system_cpu_seconds: float
                 ):
        """
        Initializes a CommandResult object.

        Args: 
            command (str): The command that was executed.
            return_code (int): The return code from the command.
            stdout (str): The standard output from the command.
            stderr (str): The standard error from the command.
            runtime_seconds (float): The runtime of the command in seconds.
            success (bool): Whether the command was successful.
            tool_name (str): The name of the tool that was executed.
            started_at (datetime): The time at which the command was started.
            ended_at (datetime): The time at which the command was completed.
            user_cpu_seconds (float): The user CPU time of the command in seconds.
            system_cpu_seconds (float): The system CPU time of the command in seconds.
        Returns:
            None
        """
        self.tool_name = tool_name
        self.started_at = started_at
        self.ended_at = ended_at
        self.command = command
        self.return_code = return_code
        self.stdout = stdout
        self.stderr = stderr
        self.runtime_seconds = runtime_seconds
        self.success = success
        self.tool_name = tool_name
        self.started_at = started_at
        self.ended_at = ended_at
        self.peak_memory_mb = peak_memory_mb
        self.user_cpu_seconds = user_cpu_seconds
        self.system_cpu_seconds = system_cpu_seconds

    @property
    def cpu_time_seconds(self) -> float | None:
        """
        Property for the total CPU time of the command in seconds.

        Returns:
            float: The total CPU time of the command in seconds.
        """
        return self.user_cpu_seconds + self.system_cpu_seconds

    def to_dict(self) -> dict:
        """
        Serializes this CommandResult to a JSON-safe dictionary (Emil's
        request, 2026-10-19: per-step compute metrics, organized and
        accessible for report generation later). stdout/stderr are
        included so a persisted metrics file also carries enough
        context to be useful on its own, without needing the tool's
        separate log file too.

        Args:
            None

        Returns:
            dict: A JSON-serializable representation of this result.
        """
        return {
            "tool_name": self.tool_name,
            "command": self.command,
            "return_code": self.return_code,
            "success": self.success,
            "runtime_seconds": self.runtime_seconds,
            "started_at": self.started_at.isoformat(),
            "ended_at": self.ended_at.isoformat(),
            "peak_memory_mb": self.peak_memory_mb,
            "user_cpu_seconds": self.user_cpu_seconds,
            "system_cpu_seconds": self.system_cpu_seconds,
            "cpu_time_seconds": self.cpu_time_seconds,
            "stdout": self.stdout,
            "stderr": self.stderr,
        }