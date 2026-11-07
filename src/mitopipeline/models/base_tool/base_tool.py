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
from typing import Callable
import shlex
import  subprocess
import threading
import time
import resource
import tempfile
import json
from mitopipeline.models.command_result.command_result import CommandResult
from mitopipeline.logging.logger_factory import TOOL_OUTPUT_EVENT

# How often the live memory-ceiling adjuster re-checks and re-applies a
# fresh fair-share budget to an already-running command (see
# dynamic_memory_limit_mb below). Frequent enough to let a job's real
# ceiling grow promptly once a sibling finishes and frees up RAM;
# infrequent enough that the check itself (a /proc scan) is never
# meaningful overhead against a job that can run for hours.
_DYNAMIC_MEMORY_RECHECK_SECONDS = 20

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
                 metrics_dir: Path | None = None,
                 memory_limit_mb: int | None = None,
                 dynamic_memory_limit_mb: Callable[[], int | None] | None = None,
                 dispatch_admission_wait: Callable[[], None] | None = None) -> None:
        """
        Instantiation method for the BaseTool class.

        Args:
            tool_name (str): The name of the tool.
            working_dir (Path): The working directory for the tool.
            logger (Logger): The logger to use for the step of the Pipeline.
            metrics_dir (Path | None): Where to persist this run's CommandResult as JSON (organized, report-ready per-step compute metrics). run() already measures CPU/memory/runtime for every tool call via GNU time. None (the default) skips persistence entirely -- existing callers that don't pass it keep working exactly as before.
            memory_limit_mb (int | None): A hard, INITIAL ceiling on this command's own virtual address space, in MB, enforced via RLIMIT_AS at process start (see run()). None (the default) applies no limit -- every existing tool wrapper that doesn't pass this keeps working exactly as before. Exists so one runaway process (e.g. GetOrganelle recruiting off-target reads) hits its OWN clean allocation failure instead of triggering the kernel's global OOM killer, which can claim an unrelated process on a shared machine. Never stops a command based on a PREDICTED pattern -- only once it's genuinely out of its real budget (Emil's explicit call, 2026-11-06, after an earlier version also tried to predict and pre-empt runaway growth from its own round-over-round pattern alone; removed).
            dynamic_memory_limit_mb (Callable[[], int | None] | None): A callback run() calls periodically (every _DYNAMIC_MEMORY_RECHECK_SECONDS) WHILE the command is still running, to re-apply a FRESH ceiling to the already-running process via resource.prlimit -- not a one-time value fixed at launch. Lets a command's own real headroom grow as sibling jobs finish and free up RAM (rather than being stuck forever at whatever a static initial division happened to allow), while still genuinely capping it against real, current system RAM -- Emil's own explicit point, 2026-11-06: a statically-divided ceiling can cut off a job that would have reached a real, natural conclusion if given the room that was genuinely free at the time. None (the default) skips this entirely -- only GetOrganelleRunner sets it.
            dispatch_admission_wait (Callable[[], None] | None): A blocking call made at the very start of run(), before this command is even built, let alone started -- lets a caller defer STARTING a new job entirely until real system state says there's genuinely enough headroom, rather than only adjusting an already-running job's own ceiling. No count-based cap on how many jobs may run at once; purely a live check of real current conditions (Emil's own explicit call, 2026-11-07, after a flat concurrency ceiling sized for worst-case safety measurably starved real throughput). None (the default) skips this entirely -- only GetOrganelleRunner sets it.

        Returns:
            None
        """
        self.tool_name = tool_name
        self.working_dir = working_dir
        self.logger = logger
        self.dynamic_memory_limit_mb = dynamic_memory_limit_mb
        self.dispatch_admission_wait = dispatch_admission_wait
        self.metrics_dir = Path(metrics_dir) if metrics_dir is not None else None
        self.memory_limit_mb = memory_limit_mb

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

        # Deferring the START of this command entirely until real system
        # state says there's genuinely enough headroom -- before even
        # validating inputs, so a job waiting on this never holds any of
        # its own resources (seed files, temp dirs) open while it waits.
        if self.dispatch_admission_wait is not None:
            self.dispatch_admission_wait()

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
        if self.memory_limit_mb is not None:
            self.logger.debug(f"{context} Enforcing a {self.memory_limit_mb}MB memory ceiling (RLIMIT_AS) on this command.")

        try:
            # Running the command. preexec_fn runs in the forked child,
            # after fork() but before exec() -- setting RLIMIT_AS there
            # caps only this one subprocess's own address space, so a
            # runaway job (e.g. GetOrganelle recruiting off-target
            # reads) hits its own clean allocation failure instead of
            # triggering the kernel's global OOM killer, which could
            # otherwise claim an unrelated process on a shared machine
            # (confirmed as a real risk, 2026-11-06 -- the kernel killed
            # a 32GB GetOrganelle process with no per-job ceiling in
            # place at all).
            #
            # Deliberately never terminates a command early based on its
            # own growth pattern alone (an earlier version of this did,
            # on a live-tailed round-growth heuristic) -- Emil's explicit
            # call, 2026-11-06: a predicted-but-unconfirmed "this might
            # get worse" is not reason enough to cut off a run that could
            # still reach a real, usable assembly. memory_limit_mb above
            # is the only thing allowed to stop a command early, and only
            # once it's genuinely, actually out of its real budget -- see
            # getorganelle.py's own post-hoc growth-pattern logging for
            # the (observe-only, never-acts) replacement.
            if self.dynamic_memory_limit_mb is not None:
                # A separate path from the plain subprocess.run() below,
                # used only when a caller (GetOrganelleRunner) wants its
                # ceiling re-adjusted LIVE while the command keeps
                # running, rather than fixed once at launch. Every other
                # tool wrapper never sets dynamic_memory_limit_mb, so
                # this branch never runs for them.
                process = subprocess.Popen(timed_command, cwd = self.working_dir, stdout = subprocess.PIPE,
                                           stderr = subprocess.PIPE, text = True, start_new_session = True,
                                           preexec_fn = self._apply_memory_limit if self.memory_limit_mb is not None else None)
                stop_event = threading.Event()
                adjuster = threading.Thread(target = self._adjust_memory_limit_periodically,
                                            args = (process, stop_event), daemon = True)
                adjuster.start()
                stdout, stderr = process.communicate()
                stop_event.set()
                adjuster.join(timeout = 2)
                completed = subprocess.CompletedProcess(args = timed_command, returncode = process.returncode,
                                                        stdout = stdout, stderr = stderr)
            else:
                # start_new_session=True puts the wrapped command (GNU
                # time + the real tool, plus every real descendant it
                # forks -- SPAdes, bowtie2, etc.) into its own process
                # group, separate from this exec script's own. This is
                # what lets an external killer (e.g. launcher.oom_monitor)
                # safely take out the ENTIRE real command subtree via a
                # single process-group kill without ever touching this
                # calling process itself -- the same proven pattern
                # already used for every real pipeline-stop tonight,
                # just scoped to one job instead of the whole pipeline.
                # Without this, a killed child's own orphaned
                # grandchildren (e.g. a bowtie2 or spades-core process
                # that doesn't die with its immediate parent) could keep
                # writing into this sample's own output directory even
                # after cleanup_failed_output() has already run --
                # confirmed as a real, recurring risk earlier this same
                # night (an orphaned bowtie2 process survived an
                # earlier pipeline-stop's own group kill).
                completed = subprocess.run(timed_command, cwd = self.working_dir, capture_output = True, text = True,
                                           start_new_session = True,
                                           preexec_fn = self._apply_memory_limit if self.memory_limit_mb is not None else None)
            runtime_seconds = time.perf_counter() - start_time
            ended_at = datetime.now()

            # Reading the execution resource statistics. On a non-zero exit,
            # GNU time prepends a human-readable notice line (e.g. "Command
            # exited with non-zero status 1") before the -f JSON line, so
            # the JSON line is parsed specifically rather than the whole file.
            # When the wrapped command is killed by a signal (e.g. an
            # external SIGKILL during a restart) rather than exiting
            # normally, GNU time never writes this file at all -- it stays
            # empty, and metrics are genuinely unavailable, not just hard
            # to parse. Confirmed live: this previously raised a bare
            # KeyError here, uncaught by anything downstream expecting a
            # documented failure (run()'s own success/failure branch never
            # even ran), so the caller's generic exception handler logged
            # "failed unexpectedly: 'user_cpu_seconds'" and returned
            # without writing any assembly status -- a real rescue attempt
            # silently burned with no record at all, rather than a
            # documented failure. Defaulting to None here instead lets the
            # real completed.returncode (already non-zero for a killed
            # process) drive the normal, already-correct failure path.
            with metrics_path.open("r", encoding = "utf-8") as handle:
                metrics_lines = [line for line in handle.readlines() if line.strip()]
            metrics = json.loads(metrics_lines[-1]) if metrics_lines else {}
            user_cpu_seconds = metrics.get("user_cpu_seconds")
            system_cpu_seconds = metrics.get("system_cpu_seconds")
            peak_memory_mb = metrics["peak_memory_mb"] / 1024 if "peak_memory_mb" in metrics else None

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

    def _find_real_child_pid(self, parent_pid: int) -> int | None:
        """
        The real, direct child PID of parent_pid -- used because the
        actual executed command is `timed_command` (/usr/bin/time --
        the-real-command), so process.pid from Popen is GNU time's own
        PID, never the real command's. Reads /proc's own direct-children
        listing (Linux-specific, matching every real deployment target
        this pipeline actually runs on). Never raises -- a parse/read
        failure here just means "not found yet," not a real error.

        Args:
            parent_pid (int): The immediate Popen PID (GNU time's own).

        Returns:
            int | None: The real wrapped command's PID, or None if it can't be found (e.g. too early, or already exited).
        """
        try:
            children_text = Path(f"/proc/{parent_pid}/task/{parent_pid}/children").read_text().strip()
            if children_text:
                return int(children_text.split()[0])
        except (OSError, ValueError):
            pass
        return None

    def _adjust_memory_limit_periodically(self, process: subprocess.Popen, stop_event: threading.Event) -> None:
        """
        While the command is still running, calls self.dynamic_memory_
        limit_mb() every _DYNAMIC_MEMORY_RECHECK_SECONDS and re-applies
        the fresh result to the REAL running child via resource.prlimit
        -- not a signal, not a kill, just a live ceiling adjustment on
        an already-running process. Lets the command's own real memory
        headroom grow as sibling jobs finish and free up RAM, instead
        of being stuck forever at whatever a one-time initial division
        happened to allow. Intended to run only as a background thread
        spawned by run(); never raises out to the caller (any read/
        adjust failure here must never take down the real tool
        invocation it's just watching).

        Args:
            process (subprocess.Popen): The running command's own process handle (GNU time's PID, not the real command -- see _find_real_child_pid).
            stop_event (threading.Event): Set by run() once the command has exited, to stop this thread promptly.

        Returns:
            None
        """
        context = self._log_context()
        if not hasattr(resource, "prlimit"):
            # prlimit is a real Linux-only syscall (prlimit64) -- not
            # available on macOS/BSD. Every deployment target this
            # pipeline actually runs on IS real Linux (fishlabworkstn
            # and WSL2 both), so this is a graceful skip for local
            # development only, never a silent failure on a target
            # that was actually supposed to work.
            self.logger.debug(f"{context} Live memory-limit adjustment skipped: resource.prlimit unavailable on this platform.")
            return
        real_pid: int | None = None
        while not stop_event.is_set() and process.poll() is None:
            if real_pid is None:
                real_pid = self._find_real_child_pid(process.pid)
            if real_pid is not None:
                new_limit_mb = self.dynamic_memory_limit_mb()
                if new_limit_mb is not None:
                    try:
                        new_bytes = new_limit_mb * 1024 * 1024
                        resource.prlimit(real_pid, resource.RLIMIT_AS, (new_bytes, new_bytes))
                        self.logger.debug(f"{context} Live-adjusted memory ceiling to {new_limit_mb}MB.")
                    except (ProcessLookupError, OSError):
                        # Already exited, or a transient permission/race
                        # issue -- never let an adjustment failure take
                        # down the real command it's just trying to help.
                        pass
            stop_event.wait(_DYNAMIC_MEMORY_RECHECK_SECONDS)

    def _apply_memory_limit(self) -> None:
        """
        Sets RLIMIT_AS on the calling process to self.memory_limit_mb --
        intended to run ONLY as subprocess.run's preexec_fn, i.e. inside
        the freshly-forked child, right before it execs the real
        command. POSIX-standard (resource.RLIMIT_AS), so this works
        unprivileged on any real Linux kernel this pipeline actually
        runs on, including WSL2.

        Args:
            None

        Returns:
            None
        """
        limit_bytes = self.memory_limit_mb * 1024 * 1024
        resource.setrlimit(resource.RLIMIT_AS, (limit_bytes, limit_bytes))

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
