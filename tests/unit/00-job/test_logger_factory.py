"""
test_logger_factory.py

Unit tests for make_logger()'s handler wiring and the
ExcludeToolOutputFilter it attaches: tool-output-tagged records must
be kept out of the console and the shared global log, but still reach
the per-sample component log in full.
"""

# Imports
import logging
from pathlib import Path
import pytest
from mitopipeline.logging.logger_factory import PIPELINE_EVENT, TOOL_OUTPUT_EVENT, ExcludeToolOutputFilter, make_logger

def test_exclude_tool_output_filter_rejects_tool_output_records() -> None:
    record = logging.LogRecord("test", logging.DEBUG, __file__, 1, "msg", None, None)
    record.event_type = TOOL_OUTPUT_EVENT
    assert ExcludeToolOutputFilter().filter(record) is False

def test_exclude_tool_output_filter_accepts_untagged_records() -> None:
    record = logging.LogRecord("test", logging.DEBUG, __file__, 1, "msg", None, None)
    assert ExcludeToolOutputFilter().filter(record) is True

def test_exclude_tool_output_filter_accepts_explicit_pipeline_events() -> None:
    record = logging.LogRecord("test", logging.DEBUG, __file__, 1, "msg", None, None)
    record.event_type = PIPELINE_EVENT
    assert ExcludeToolOutputFilter().filter(record) is True

def test_tool_output_tagged_messages_reach_the_component_log_but_not_the_global_log(tmp_path: Path) -> None:
    """Regression: a raw tool stdout dump belongs in the per-sample
    component log (full detail for debugging that one job), but not
    in the global log shared across every parallel sample's jobs."""
    component_log = tmp_path / "component.log"
    global_log = tmp_path / "global.log"

    logger = make_logger(name = "test_tool_output_wiring", log_file_path = component_log, global_log_file_path = global_log)
    logger.debug("tagged tool console output", extra = {"event_type": TOOL_OUTPUT_EVENT})
    logger.info("ordinary pipeline message")
    for handler in logger.handlers:
        handler.flush()

    component_text = component_log.read_text(encoding = "utf-8")
    global_text = global_log.read_text(encoding = "utf-8")

    assert "tagged tool console output" in component_text
    assert "ordinary pipeline message" in component_text
    assert "tagged tool console output" not in global_text
    assert "ordinary pipeline message" in global_text
