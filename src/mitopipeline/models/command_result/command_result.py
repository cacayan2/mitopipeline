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
        tool_version (str | None): The real external tool's own version string (from its own --version-equivalent output), or None when BaseTool.version_command() wasn't overridden for this tool or the real check failed.
        retry_attempts (list[dict]): Real per-attempt summaries for every attempt BEFORE this final one (empty for every command that never retried) -- see BaseTool.run()'s own oom_retry_decision. Each entry: {"attempt", "return_code", "runtime_seconds", "started_at", "ended_at", "peak_memory_mb", "user_cpu_seconds", "system_cpu_seconds"}. peak_memory_mb/cpu fields are None for a genuinely killed attempt -- GNU time never gets to write its own metrics file when the wrapped command dies by signal (see run()'s own comment on this), so there is real data to report for runtime/timing but honestly none for memory/CPU on that specific attempt. The real peak RSS at the moment of a kill is recorded separately, by launcher.oom_monitor itself (see utils.system.resources.record_oom_kill) -- not duplicated here, to avoid two sources of truth for the same real number. Deliberately kept SEPARATE from the top-level runtime_seconds/peak_memory_mb/cpu fields, never summed into them -- Emil's own explicit call, 2026-11-07: the real report-facing metadata must reflect the actual GetOrganelle instance that produced the assembly (the final attempt alone), with failed attempts available for reference, not blended into a combined total that would muddy what the final attempt itself actually cost.
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
                 system_cpu_seconds: float,
                 tool_version: str | None = None,
                 retry_attempts: list[dict] | None = None
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
            retry_attempts (list[dict] | None): See the class docstring. None (the default) becomes an empty list -- every existing caller that never retries keeps working exactly as before.
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
        self.tool_version = tool_version
        self.retry_attempts = retry_attempts if retry_attempts is not None else []

    @property
    def cpu_time_seconds(self) -> float | None:
        """
        Property for the total CPU time of the command in seconds.
        None when either component is unavailable (e.g. the command was
        killed by a signal before GNU time could record anything) --
        the type hint already promised this; the implementation just
        didn't honor it, and would otherwise crash adding None to a
        float the moment that case was ever actually reached.

        Returns:
            float | None: The total CPU time of the command in seconds, or None if unavailable.
        """
        if self.user_cpu_seconds is None or self.system_cpu_seconds is None:
            return None
        return self.user_cpu_seconds + self.system_cpu_seconds

    def to_dict(self) -> dict:
        """
        Serializes this CommandResult to a JSON-safe dictionary,
        producing organized, report-ready per-step compute metrics.
        stdout/stderr are
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
            "tool_version": self.tool_version,
            "retry_attempts": self.retry_attempts,
            "stdout": self.stdout,
            "stderr": self.stderr,
        }