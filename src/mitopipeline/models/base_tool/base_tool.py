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
import os
import re
import shlex
import signal
import  subprocess
import threading
import time
import resource
import tempfile
import json
from mitopipeline.models.command_result.command_result import CommandResult
from mitopipeline.logging.logger_factory import TOOL_OUTPUT_EVENT

# GetOrganelle's own real extension-round log line: "Round N: X/X AI <int>
# AW <int>" (confirmed directly against a real completed run's
# get_org.log.txt, 2026-11-06). AW is the real signal tracked -- the
# documented explosive-off-target-recruitment signature this session has
# repeatedly observed (percina_sciera, cyprinus_carpio) is an AW jump of
# multiple ORDERS OF MAGNITUDE between consecutive rounds, not the mild
# 2-3x acceleration a real, eventually-successful run can still show
# (confirmed live: a real completed sample's own round 6 jumped ~2.8x
# over round 5 and still finished normally). 20x sits well above that
# real survivable case and well below the catastrophic one, so it
# shouldn't false-trigger on ordinary accelerating growth.
_ROUND_GROWTH_PATTERN = re.compile(r"Round (\d+): \d+/\d+ AI (\d+) AW (\d+)")
_GROWTH_RATIO_THRESHOLD = 20.0

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
                 growth_monitor_log: Path | None = None) -> None:
        """
        Instantiation method for the BaseTool class.

        Args:
            tool_name (str): The name of the tool.
            working_dir (Path): The working directory for the tool.
            logger (Logger): The logger to use for the step of the Pipeline.
            metrics_dir (Path | None): Where to persist this run's CommandResult as JSON (organized, report-ready per-step compute metrics). run() already measures CPU/memory/runtime for every tool call via GNU time. None (the default) skips persistence entirely -- existing callers that don't pass it keep working exactly as before.
            memory_limit_mb (int | None): A hard ceiling on this command's own virtual address space, in MB, enforced via RLIMIT_AS (see run()). None (the default) applies no limit -- every existing tool wrapper that doesn't pass this keeps working exactly as before. Exists so one runaway process (e.g. GetOrganelle recruiting off-target reads) hits its OWN clean allocation failure instead of triggering the kernel's global OOM killer, which can claim an unrelated process on a shared machine.
            growth_monitor_log (Path | None): GetOrganelle's own native log file (get_org.log.txt) to watch live for its real per-round "AI ... AW ..." extension output. When set, run() tails this file in a background thread while the command is running and sends SIGTERM to the whole process group the moment a round's AW grows past _GROWTH_RATIO_THRESHOLD over the previous round -- a real, confirmed signature of off-target recruitment running away (see the module-level comment). Stops it EARLY, with memory use still modest, rather than letting it balloon toward an uncontrolled OOM the hard memory_limit_mb ceiling would otherwise have to absorb later. None (the default, every non-GetOrganelle caller) skips this entirely.

        Returns:
            None
        """
        self.tool_name = tool_name
        self.working_dir = working_dir
        self.logger = logger
        self.metrics_dir = Path(metrics_dir) if metrics_dir is not None else None
        self.memory_limit_mb = memory_limit_mb
        self.growth_monitor_log = Path(growth_monitor_log) if growth_monitor_log is not None else None

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
        if self.memory_limit_mb is not None:
            self.logger.debug(f"{context} Enforcing a {self.memory_limit_mb}MB memory ceiling (RLIMIT_AS) on this command.")

        try:
            if self.growth_monitor_log is not None:
                # A separate path from the plain subprocess.run() below,
                # used only when a caller (GetOrganelleRunner) actually
                # wants live round-growth monitoring -- every other tool
                # wrapper never sets growth_monitor_log, so this branch
                # never runs for them and their behavior is unchanged.
                # start_new_session=True puts this command in its OWN
                # process group (GNU time plus whatever it execs), so a
                # detected runaway round can be stopped by signalling the
                # whole group -- signalling just the top PID would only
                # ever reach /usr/bin/time itself, not the real command
                # it wraps, since GNU time doesn't forward signals to
                # its child on its own.
                process = subprocess.Popen(timed_command, cwd = self.working_dir, stdout = subprocess.PIPE,
                                           stderr = subprocess.PIPE, text = True, start_new_session = True,
                                           preexec_fn = self._apply_memory_limit if self.memory_limit_mb is not None else None)
                stop_event = threading.Event()
                monitor = threading.Thread(target = self._monitor_round_growth, args = (process, stop_event), daemon = True)
                monitor.start()
                stdout, stderr = process.communicate()
                stop_event.set()
                monitor.join(timeout = 2)
                completed = subprocess.CompletedProcess(args = timed_command, returncode = process.returncode,
                                                        stdout = stdout, stderr = stderr)
            else:
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
                completed = subprocess.run(timed_command, cwd = self.working_dir, capture_output = True, text = True,
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

    def _monitor_round_growth(self, process: subprocess.Popen, stop_event: threading.Event) -> None:
        """
        Tails self.growth_monitor_log (GetOrganelle's own get_org.log.txt)
        every 5s while the command runs, watching for a real, confirmed
        off-target-recruitment signature: a round's AW jumping past
        _GROWTH_RATIO_THRESHOLD over the previous round. On a real hit,
        signals the whole process group (not just the Popen PID -- see
        run()'s own comment on why) so the command exits NOW, with
        memory use still modest, rather than being left to balloon
        further toward an uncontrolled OOM. Intended to run only as a
        background thread spawned by run(); never raises out to the
        caller (any read/parse failure here must never take down the
        real tool invocation it's just watching).

        Args:
            process (subprocess.Popen): The running command's own process handle.
            stop_event (threading.Event): Set by run() once the command has exited, to stop this thread promptly.

        Returns:
            None
        """
        context = self._log_context()
        seen_rounds: set[int] = set()
        last_aw: int | None = None
        while not stop_event.is_set() and process.poll() is None:
            try:
                text = self.growth_monitor_log.read_text(encoding = "utf-8", errors = "replace") if self.growth_monitor_log.is_file() else ""
            except OSError:
                text = ""
            for match in _ROUND_GROWTH_PATTERN.finditer(text):
                round_num = int(match.group(1))
                if round_num in seen_rounds:
                    continue
                seen_rounds.add(round_num)
                aw = int(match.group(3))
                if last_aw is not None and last_aw > 0:
                    ratio = aw / last_aw
                    if ratio >= _GROWTH_RATIO_THRESHOLD:
                        self.logger.warning(
                            f"{context} Round {round_num}'s AW grew {ratio:.1f}x over the previous round "
                            f"({last_aw} -> {aw}) -- past the real explosive-recruitment threshold. "
                            f"Terminating now, before it can balloon further."
                        )
                        try:
                            os.killpg(os.getpgid(process.pid), signal.SIGTERM)
                        except (ProcessLookupError, PermissionError):
                            pass
                        return
                last_aw = aw
            time.sleep(5)

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
