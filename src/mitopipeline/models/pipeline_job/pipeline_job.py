"""
pipeline_job.py

Pipeline job lifecycle and metadata tracking.
"""

# Imports
from __future__ import annotations
from datetime import datetime
from pathlib import Path
import json
import logging
import uuid
from mitopipeline.logging.logger_factory import make_logger

class PipelineJob:
    """
    Class which manages the lifecycle and job metadata for a single pipeline run.
    
    Attributes:
        parent_dir (Path): The parent directory of the pipeline run.
        job_id (str | None = None): The unique identifier for the pipeline job.
        samples: list[str] | None = None: A list of sample IDs associated with the pipeline job.
    """
    def __init__(self,
                 parent_dir: Path,
                 job_id: str | None = None,
                 samples: list[str] | None = None):
        """
        Initializes a pipeline job. Samples may be unknown when the launcher first creates the job.
        In that case, the job starts with an empty list and lets the launcher update samples after preparing 
        the runtime manifest.

        Args:
            parent_dir (Path): The parent directory of the pipeline run.
            job_id (str | None = None): The unique identifier for the pipeline job.
            samples (list[str] | None = None): A list of sample IDs associated with the pipeline job.
        
        Returns:
            None
        """
        # We will detect if the job is new by detecting if the output directory exists and matches
        # the name listed in the config. Set to True by default. If the job already exists, set to False
        # and generate a new job_id using UUID generation and current date and time. 
        self.is_new_job = True
        if job_id is not None:
            self.is_new_job = False
            self.job_id = job_id
        else:
            now = datetime.now()
            self.job_id = (f"{now.strftime('%Y%m%d')}_{now.strftime('H%M')}_{uuid.uuid4()}")
    
        # Now we initialize samples based on the description provided in the
        # docstring for this function.
        self.samples = (list(samples) if samples is not None else [])

        # Instantiation of runtime metadata.
        self.start_time = datetime.now()
        self.end_time = None
        self.created_time = None
        self.failed_time = None
        self.completed_time = None
        self.runtime = None
        self.parent_dir = Path(parent_dir)
        self.job_dir = (self.parent_dir / self.job_id)
        self.job_logger = self.create_output_directory()
        self.job_logger.info(f"Starting pipeline job {self.job_id} at {self.start_time}.")

        # Marking the job as created.
        self.mark_created()

    def create_output_directory(self) -> logging.Logger:
        """
        This creates the output directories required for the job and returns a job logger.

        Args:
            None

        Returns:
            logging.Logger: A logger for the job.
        """
        # Making the parent directory and the job directory if they do not exist.
        self.parent_dir.mkdir(parents = True, exist_ok = True)
        self.job_dir.mkdir(parents = True, exist_ok = True)

        # Setting up directories for all the outputs of the pipeline. 
        output_directories = [
            self.job_dir / "01-qc_raw",
            self.job_dir / "02-trimming",
            self.job_dir / "03-qc_trimmed",
            self.job_dir / "04a-getorganelle_config",
            self.job_dir / "04b-seed_sequences",
            self.job_dir / "04c-assembly",
            self.job_dir / "05a-blast_config",
            self.job_dir / "05b-blast",
            self.job_dir / "05c-blast_hits",
            self.job_dir / "06-contig_selection",
            self.job_dir / "07-assembly_rescue",
            self.job_dir / "08a-mitos2_config",
            self.job_dir / "08b-annotation",
            self.job_dir / "09-visualization",
            self.job_dir / "10-alignment_dataset",
            self.job_dir / "11-phylogeny",
            self.job_dir / "12-reporting",
            self.job_dir / "13-report_pdf",
            self.job_dir / "14a-sra_archival",
            self.job_dir / "14b-biosample_archival",
            self.job_dir / "15-archival_validation",
            self.job_dir / "16-table2asn_v"
        ]

        # Setting up all the directories required for logging of the pipeline.
        log_directories = [
            self.job_dir / "logs",
            self.job_dir / "logs" / "00-job",
            self.job_dir / "logs" / "01-qc_raw",
            self.job_dir / "logs" / "02-trimming",
            self.job_dir / "logs" / "03-qc_trimmed",
            self.job_dir / "logs" / "04-assembly",
            self.job_dir / "logs" / "05-blast",
            self.job_dir / "logs" / "06-contig_selection",
            self.job_dir / "logs" / "07-assembly_rescue",
            self.job_dir / "logs" / "08-annotation",
            self.job_dir / "logs" / "09-visualization",
            self.job_dir / "logs" / "10-alignment_dataset",
            self.job_dir / "logs" / "11-phylogeny",
            self.job_dir / "logs" / "12-reporting",
            self.job_dir / "logs" / "13-report_pdf",
            self.job_dir / "logs" / "14-archival",
            self.job_dir / "logs" / "15-archival_validation",
            self.job_dir / "logs" / "16-table2asn_v"
        ]

        # Creating the output directories.
        for directory in output_directories: directory.mkdir(parents = True, exist_ok = True)
        for directory in log_directories: directory.mkdir(parents = True, exist_ok = True)

        # Stored (not just passed to make_logger) so mark_failed()/
        # mark_completed() can reference it -- they already did, but the
        # attribute was never actually set, an AttributeError on every
        # real failure/completion until now.
        self.global_log_file_path = self.job_dir / "logs" / "00-job" / f"{self.job_id}.global.log"

        # Returning the logger for the job.
        return make_logger(
            name = f"{self.job_id}.job",
            log_file_path = self.job_dir / "logs" / f"{self.job_id}.job.log",
            global_log_file_path = self.global_log_file_path,
            console_level = logging.INFO,
            file_level = logging.DEBUG
        )

    def set_samples(self, samples: list[str]) -> None:
        """
        Defines the list of samples to be utilized for the job. 

        Args:
            samples (list[str]): A list of sample IDs to be used for the job.

        Returns:
            None
        """
        self.samples = list(samples)

    def _metadata_payload(self, event_time: datetime, status: str,
                          failure_reason: str | None = None,
                          log_tail: str | None = None) -> dict:
        """
        Returns a dictionary of metadata to be logged for the job.

        Args:
            event_time (datetime): The time at which the event occurred.
            status (str): The status of the event.
            failure_reason (str | None): A short, specific description of why the job failed (which rule, what error) -- only meaningful when status is "failed".
            log_tail (str | None): The tail of captured Snakemake stdout/stderr around the failure, for a reader who needs more than the short reason -- only meaningful when status is "failed".

        Returns:
            dict: A dictionary of metadata to be logged for the job.
        """
        # Calculates the runtime of the job.
        runtime = (event_time - self.start_time)

        return {
            "job_id": self.job_id,
            # .isoformat(), not the raw datetime -- json.dump cannot
            # serialize a datetime object at all; this crashed
            # write_metadata() on every single call, every status, not
            # just failures.
            "start_time": self.start_time.isoformat(),
            "event_time": event_time.isoformat(),
            "runtime_seconds": runtime.total_seconds(),
            "output_dir": str(self.job_dir),
            "sample_count": len(self.samples),
            "samples": self.samples,
            "status": status,
            "failure_reason": failure_reason,
            "log_tail": log_tail,
        }

    def write_metadata(self, data: dict) -> None:
        """
        For a given metadata payload (data: dict), we write the metadata to a metadata file.
        This will be done in json format.

        Args:
            data (dict): A dictionary of metadata to be logged for the job.

        Returns:
            None
        """
        # Live file represents the metadata for the CURRENT run of the pipeline. There will also be a historical metadata file
        # that is constantly appended to upon new runs of the pipeline.
        live_file = self.job_dir / "current_metadata.json"
        history_file = self.job_dir / "historical_metadata.json"

        # Here we overwrite the live file.
        with live_file.open("w", encoding = "utf-8") as handle:
            json.dump(data, handle, indent = 4)

        # For the historical file, some logic must be done:
        # 1. First we check if the historical file exists - if not, we create an empty file.
        # 2. Then we load the historical file into memory. 
        # 3. Then we append the new data to the historical file. 
        # 4. Then we move the file pointer of the handle to the beginning of the file to read from the beginning. 
        # 5. Since we have both the historical data and the new data in memory, we dump everything into the file. 
        # 6. We remove any data that comes after our dump with the handle.truncate() method, such that all the contents of the file
        # are from the historical data loaded from memory. 
        ## Such that each time we run, even if the job is new, we are essentially appending to the historical file.
        ## There is a try-except block to account for errors in reading the data - if so, then the historical file
        ## cannot be loaded and will effectively be overwritten (it's a little more complicated than that,
        ## but I hope my implementation makes it fairly clear what I mean by "essentially overwritten").
        if not history_file.exists():
            history_file.write_text("[]", encoding = "utf-8")
        with history_file.open("r+", encoding = "utf-8") as handle:
            try:
                historical_data = json.load(handle)
            except json.JSONDecodeError:
                historical_data = []
            historical_data.append(data)
            handle.seek(0)
            json.dump(historical_data, handle, indent = 4)
            handle.truncate()
        self.job_logger.info(f"Metadata for job {self.job_id} has been written to {live_file} and {history_file}.")

    def mark_created(self) -> None:
        """
        This function marks the job as created - this serves purely for the sake of generating metadata.
        
        Args:
            None

        Returns: 
            None
        """
        # First we get the current date and time, alongisde the runtime up to this point.
        self.created_time = datetime.now()

        # Then we obtain a metadata payload for this event.
        data = self._metadata_payload(
            event_time = self.created_time,
            status = "created"
        )

        # Afterwards we detect logic if the job is a new job or not and log accordingly.
        if self.is_new_job:
            self.job_logger.info(f"Job {self.job_id} is a new job - new directoy was successfully created at {self.job_dir}.")
        else:
            self.job_logger.info(f"Pipeline job {self.job_id} already exists, using existing directory at {self.job_dir} and automatically progressing pipeline from previous run.")

        self.write_metadata(data)

    def mark_failed(self, failure_reason: str | None = None,
                    log_tail: str | None = None) -> None:
        """
        Mark the pipeline job as failed.

        Args:
            failure_reason (str | None): A short, specific description of why the job failed (e.g. which rule failed and its exit code), when known.
            log_tail (str | None): The tail of captured Snakemake stdout/stderr around the failure, when captured.

        Returns:
            None
        """
        # First we get the current date and time.
        self.failed_time = datetime.now()

        # Then we obtain a metadata payload for this event.
        data = self._metadata_payload(
            event_time = self.failed_time,
            status = "failed",
            failure_reason = failure_reason,
            log_tail = log_tail
        )

        # Logging the failure of the job and writing metadata.
        self.job_logger.error(
            f"Pipeline job {self.job_id} has failed at {self.failed_time}"
            + (f": {failure_reason}" if failure_reason else ".") +
            f" Job details can be found in {self.job_dir}. "
            f"Global log file can be found at {self.global_log_file_path}.")
        self.write_metadata(data)

    def mark_completed(self) -> None:
        """
        Mark the pipeline job as completed.
        
        Args:
            None

        Returns:
            None
        """
        # First we get the current date and time.
        self.completed_time = datetime.now()

        # Then we obtain a metadata payload for this event.
        data = self._metadata_payload(
            event_time = self.completed_time,
            status = "completed"
        )

        # Logging the completion of the job and writing metadata.
        self.job_logger.info(
            f"Pipeline job {self.job_id} has completed and completed at {self.completed_time}. Job details can be found in {self.job_dir}. "
            f"Global log file can be found at {self.global_log_file_path}.")
        self.write_metadata(data)

    def aggregate_history(self) -> dict:
        """
        Summarizes this job's historical_metadata.json across every
        time it has been launched with this same job_id, covering a
        job resumed/retried several times before finally completing:
        total attempts, cumulative runtime across all of them, and a
        per-status breakdown, ready for 12a-stats_aggregation/
        13-report_pdf to include in a report.

        An empty or unreadable history is reported as zero attempts
        rather than raising -- a report generator should get a quiet,
        honest "nothing here yet", not a crash, over a job that simply
        hasn't failed or been retried before.

        Args:
            None

        Returns:
            dict: {"job_id", "attempt_count", "cumulative_runtime_seconds", "status_counts" (a dict of status -> count), "last_status" (the most recent attempt's status, or None)}.
        """
        history_file = self.job_dir / "historical_metadata.json"
        try:
            historical_data = json.loads(history_file.read_text(encoding = "utf-8"))
        except (OSError, json.JSONDecodeError):
            historical_data = []

        # Every real launch writes a "created" entry AND, once
        # run_snakemake finishes, exactly one terminal "completed" or
        # "failed" entry -- "created" is a lifecycle marker, not an
        # outcome, so it's excluded here to avoid double-counting one
        # attempt as two. (A "created" with no matching terminal entry
        # after it means the process itself died mid-run, e.g. a power
        # loss -- counted separately below rather than silently ignored.)
        terminal = [entry for entry in historical_data if entry.get("status") in ("completed", "failed")]
        created_count = sum(1 for entry in historical_data if entry.get("status") == "created")

        status_counts: dict[str, int] = {}
        cumulative_runtime_seconds = 0.0
        for entry in terminal:
            status = entry["status"]
            status_counts[status] = status_counts.get(status, 0) + 1
            cumulative_runtime_seconds += float(entry.get("runtime_seconds") or 0.0)

        return {
            "job_id": self.job_id,
            "attempt_count": len(terminal),
            "cumulative_runtime_seconds": cumulative_runtime_seconds,
            "status_counts": status_counts,
            "last_status": terminal[-1]["status"] if terminal else None,
            # Started but never reached a terminal status (process died
            # mid-run) -- a real gap a report should be able to flag.
            "incomplete_attempts": max(0, created_count - len(terminal)),
        }