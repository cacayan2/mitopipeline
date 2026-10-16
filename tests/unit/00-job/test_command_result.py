"""
test_command_result.py

Unit tests for the CommandResult class.
"""

# Imports
from datetime import datetime
from mitopipeline.models.command_result.command_result import CommandResult

def _make_command_result(**overrides) -> CommandResult:
    """
    Builds a CommandResult with sensible defaults, overridable per test.

    Args:
        **overrides: Keyword arguments overriding the defaults.

    Returns:
        CommandResult: The constructed CommandResult.
    """
    defaults = dict(
        command = ["fastp", "--in1", "r1.fastq"],
        return_code = 0,
        stdout = "ok",
        stderr = "",
        runtime_seconds = 1.5,
        success = True,
        tool_name = "fastp",
        started_at = datetime(2026, 10, 16, 0, 0, 0),
        ended_at = datetime(2026, 10, 16, 0, 0, 1),
        peak_memory_mb = 128.0,
        user_cpu_seconds = 0.4,
        system_cpu_seconds = 0.1,
    )
    defaults.update(overrides)
    return CommandResult(**defaults)

def test_init_stores_all_fields() -> None:
    """
    All constructor arguments should be stored as attributes unchanged.
    """
    command_result = _make_command_result()
    assert command_result.command == ["fastp", "--in1", "r1.fastq"]
    assert command_result.return_code == 0
    assert command_result.stdout == "ok"
    assert command_result.stderr == ""
    assert command_result.runtime_seconds == 1.5
    assert command_result.success is True
    assert command_result.tool_name == "fastp"
    assert command_result.peak_memory_mb == 128.0
    assert command_result.user_cpu_seconds == 0.4
    assert command_result.system_cpu_seconds == 0.1

def test_cpu_time_seconds_sums_user_and_system() -> None:
    """
    cpu_time_seconds should be the sum of user and system CPU time.
    """
    command_result = _make_command_result(user_cpu_seconds = 0.4, system_cpu_seconds = 0.1)
    assert command_result.cpu_time_seconds == 0.5

def test_cpu_time_seconds_is_a_class_property() -> None:
    """
    cpu_time_seconds must be reachable on every instance, not just one built
    inside __init__ -- a regression check for a prior bug where the
    property was defined (and discarded) inside the constructor body.
    """
    first = _make_command_result(user_cpu_seconds = 1.0, system_cpu_seconds = 2.0)
    second = _make_command_result(user_cpu_seconds = 5.0, system_cpu_seconds = 0.0)
    assert first.cpu_time_seconds == 3.0
    assert second.cpu_time_seconds == 5.0
