"""
test_pipeline_job.py

Unit tests for the PipelineJob class.
"""

# Imports
import json
from pathlib import Path
import pytest
from mitopipeline.models.pipeline_job.pipeline_job import PipelineJob


def _make_job(tmp_path: Path, job_id: str | None = None) -> PipelineJob:
    """
    Builds a PipelineJob under a real temporary directory.

    Args:
        tmp_path (Path): The pytest-provided temporary directory.
        job_id (str | None): An explicit job_id, to simulate resuming an existing job.

    Returns:
        PipelineJob: The constructed job.
    """
    return PipelineJob(parent_dir = tmp_path / "jobs", job_id = job_id)


def test_mark_created_writes_valid_json_metadata(tmp_path: Path) -> None:
    """
    Instantiating a job (which calls mark_created()) must actually
    succeed and leave behind valid, re-loadable JSON -- previously
    json.dump's ident= typo crashed this on every single call.
    """
    job = _make_job(tmp_path)
    live_file = job.job_dir / "current_metadata.json"
    history_file = job.job_dir / "historical_metadata.json"
    assert live_file.is_file()
    assert history_file.is_file()
    live_data = json.loads(live_file.read_text())
    assert live_data["status"] == "created"
    assert live_data["job_id"] == job.job_id
    # start_time must be a real ISO string, not a raw datetime -- the
    # second bug that independently crashed json.dump() here.
    assert isinstance(live_data["start_time"], str)
    history_data = json.loads(history_file.read_text())
    assert len(history_data) == 1


def test_global_log_file_path_is_set(tmp_path: Path) -> None:
    """
    mark_failed()/mark_completed() both reference
    self.global_log_file_path in their log message -- it must actually
    be set as an attribute, not just handed to make_logger() and
    forgotten (an AttributeError on every real failure/completion
    before this).
    """
    job = _make_job(tmp_path)
    assert job.global_log_file_path == job.job_dir / "logs" / "00-job" / f"{job.job_id}.global.log"


def test_mark_failed_does_not_raise_and_records_reason(tmp_path: Path) -> None:
    """
    mark_failed() must run to completion (no AttributeError from the
    missing global_log_file_path) and persist the failure reason and
    log tail it was given.
    """
    job = _make_job(tmp_path)
    job.mark_failed(failure_reason = "Error in rule assembly: exit code 1", log_tail = "...tail of output...")
    live_data = json.loads((job.job_dir / "current_metadata.json").read_text())
    assert live_data["status"] == "failed"
    assert live_data["failure_reason"] == "Error in rule assembly: exit code 1"
    assert live_data["log_tail"] == "...tail of output..."


def test_mark_completed_does_not_raise(tmp_path: Path) -> None:
    """
    mark_completed() must also run to completion without the missing
    global_log_file_path attribute error.
    """
    job = _make_job(tmp_path)
    job.mark_completed()
    live_data = json.loads((job.job_dir / "current_metadata.json").read_text())
    assert live_data["status"] == "completed"
    assert live_data["failure_reason"] is None


def test_aggregate_history_with_no_history_reports_zero_attempts(tmp_path: Path) -> None:
    """
    A brand-new job's history is just its own "created" entry so far
    (no terminal outcome yet) -- reported as zero completed attempts,
    not an error.
    """
    job = _make_job(tmp_path)
    summary = job.aggregate_history()
    assert summary["attempt_count"] == 0
    assert summary["status_counts"] == {}
    assert summary["last_status"] is None
    assert summary["incomplete_attempts"] == 1


def test_aggregate_history_across_multiple_attempts(tmp_path: Path) -> None:
    """
    Resuming the same job_id several times (some failing, one
    eventually completing) must aggregate correctly: attempts counted
    by terminal outcome only (not double-counted with "created"),
    cumulative runtime summed, and the last status reflecting the most
    recent terminal outcome.
    """
    parent_dir = tmp_path / "jobs"
    job_id = "20261019_test_job"

    first = PipelineJob(parent_dir = parent_dir, job_id = job_id)
    first.mark_failed(failure_reason = "Error in rule trimming: exit code 1")

    second = PipelineJob(parent_dir = parent_dir, job_id = job_id)
    second.mark_failed(failure_reason = "Error in rule assembly: exit code 1")

    third = PipelineJob(parent_dir = parent_dir, job_id = job_id)
    third.mark_completed()

    summary = third.aggregate_history()
    assert summary["attempt_count"] == 3
    assert summary["status_counts"] == {"failed": 2, "completed": 1}
    assert summary["last_status"] == "completed"
    assert summary["incomplete_attempts"] == 0
    assert summary["cumulative_runtime_seconds"] >= 0


def test_write_metadata_appends_to_history_without_losing_prior_entries(tmp_path: Path) -> None:
    """
    Each write_metadata() call must append to historical_metadata.json,
    not overwrite it -- the whole point of the history file.
    """
    job = _make_job(tmp_path)
    job.mark_failed(failure_reason = "first failure")
    job.mark_failed(failure_reason = "second failure")
    history_data = json.loads((job.job_dir / "historical_metadata.json").read_text())
    assert len(history_data) == 3  # created + two mark_failed calls
    assert history_data[-1]["failure_reason"] == "second failure"
    assert history_data[-2]["failure_reason"] == "first failure"
