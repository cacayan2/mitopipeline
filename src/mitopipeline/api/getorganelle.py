"""
getorganelle.py

GetOrganelle wrapper for command building and execution.
"""

# Imports
from __future__ import annotations
from datetime import datetime
from logging import Logger
from pathlib import Path
import json
import os
import shutil
import sys
import time
from dataclasses import asdict, dataclass
from mitopipeline.models.base_tool.base_tool import BaseTool
from mitopipeline.models.sample.sample import Sample
from mitopipeline.utils.seed.models import SeedSelectionResult
from mitopipeline.utils.system.resources import (
    GETORGANELLE_ADMISSION_FLOOR_MB,
    acquire_exclusive_trial,
    cgroup_memory_max_mb,
    current_cgroup_path,
    exclusive_trial_lock_path,
    job_root_from_global_log,
    oom_state_path,
    read_cgroup_oom_kill_count,
    read_oom_escalation_record,
    read_reserved_mb,
    record_oom_kill,
    recommended_concurrent_slots,
    thread_count_for_reserved_slots,
    wait_for_dispatch_admission,
)
from mitopipeline.api.fastqc import FastQCRunner

# How long run_oom_monitor's own real kill-and-record sequence (one
# synchronous step inside its own process, triggered by its own 5s
# poll cycle) gets to actually land in oom_escalation_state.json
# before _oom_retry_decision gives up waiting for it. Real margin
# above that 5s cycle, not a guess -- generous enough that a real
# write under genuine system load still has room to land.
_OOM_ESCALATION_RECORD_WAIT_SECONDS = 15.0

# The organelle types GetOrganelle itself supports via -F.
VALID_ORGANELLE_TYPES = {
    "animal_mt",
    "embplant_pt",
    "embplant_mt",
    "fungus_mt",
    "fungus_nr",
    "embplant_nr",
    "anonym",
}

class GetOrganelleRunner(BaseTool):
    """
    This class provides the functionality to run GetOrganelle and
    normalize I/O within a single sample directory.

    Attributes:
        working_dir (Path): The path to the directory of the input files.
        output_dir (Path): The directory to save the output files to.
        sample (Sample): The sample object that is being worked with.
        organelle_type (str): The GetOrganelle organelle type to assemble (e.g. "animal_mt").
        tool_options (dict): A dictionary of options to pass to the tool. A "seed_file" entry (-s) is what makes this a seeded run; its absence is a valid no-seed run, falling back to GetOrganelle's own default database for organelle_type.
        logger (Logger): The logger to use for this of the pipeline.
        tool_name (str): The name of the tool (by default will be getorganelle).
        threads (int): The number of the threads for this process to use.
        standardized_fasta (Path): The path to the standardized assembly FASTA output.
        standardized_gfa (Path): The path to the standardized assembly graph output.
    """
    def __init__(self,
                 working_dir: Path,
                 output_dir: Path,
                 sample: Sample,
                 organelle_type: str,
                 tool_options: dict | None = None,
                 threads: int = 4,
                 logger: Logger | None = None,
                 tool_name: str = "getorganelle",
                 metrics_dir: Path | None = None,
                 global_log_file: Path | None = None) -> None:
        """
        Instantiation method for the creation of the
        GetOrganelleRunner class. Calls the BaseTool instantiation methods
        and sets a few other object variables.

        Args:
            working_dir (Path): The path to the directory of the input files.
            output_dir (Path): The directory to save the output files to.
            sample (Sample): The sample object that is being worked with.
            organelle_type (str): The GetOrganelle organelle type to assemble.
            tool_options (dict): A dictionary of options to pass to the tool.
            threads (int): The number of the threads for this process to use.
            logger (Logger): The logger to use for this of the pipeline.
            tool_name (str): The name of the tool (by default will be getorganelle).
            metrics_dir (Path | None): Where to persist this run's compute metrics (CPU, memory, runtime) as JSON; None skips persistence. See BaseTool.
            global_log_file (Path | None): This job's global log file, as already passed to every exec script -- used only to locate this job's own real OOM-kill reservation state (see utils.system.resources.read_reserved_mb). None (e.g. a caller outside the real pipeline, like a test) means every dispatch just uses the normal, un-escalated floor.

        Returns:
            None
        """
        # No per-job memory ceiling of any kind (neither a one-time
        # memory_limit_mb nor a live-adjusting dynamic_memory_limit_mb) --
        # removed entirely, 2026-11-07, on Emil's explicit call. Real data
        # from the fresh run showed the live fair-share ceiling this
        # replaced was cutting off a real plurality of samples (roughly
        # half of completed primaries) specifically during pre-grouping,
        # well below what GetOrganelle actually needed, with zero
        # reliable way found tonight to predict in advance which ones
        # would need more (four separate real, validated attempts --
        # early round-growth, fastp signals, NUMT screening, and read
        # count -- each failed against real historical/live data). The
        # real backstop now lives one level up: the whole pipeline is
        # launched inside a single memory-capped systemd --user --scope
        # cgroup (see launcher), so the only thing left to cap any one
        # job's growth is genuine, real exhaustion of that shared pool --
        # and the kernel's own OOM killer is scoped to that cgroup alone,
        # never able to touch anything outside the pipeline's own
        # processes. A real OOM kill is then caught by launcher.oom_monitor,
        # which escalates this sample's own real reservation (see
        # record_oom_kill, utils.system.resources) and intelligently
        # requeued in place -- retried for real, with that same
        # escalated room, over and over, until it's proven to genuinely
        # exceed this machine's real capacity (killed again even while
        # running with the whole machine reserved to itself) -- rather
        # than guessed at or prevented in advance, and never given up
        # on after an arbitrary fixed number of tries.
        #
        # dispatch_admission_wait stays on -- it still avoids wasting the
        # shared pool on simultaneous cold starts, which is a real,
        # orthogonal efficiency concern, not a prediction about any one
        # job's eventual need. logger is captured by reference here, not
        # called yet -- by the time run() actually invokes this closure,
        # __init__ has already finished and self.logger is set.
        #
        # The admission floor itself is this sample's own real,
        # escalated reservation when one is on record (see
        # read_reserved_mb/record_oom_kill) -- GETORGANELLE_ADMISSION_
        # FLOOR_MB (the normal, un-escalated share) for every sample
        # that's never been killed, which is almost all of them.
        # global_log_file is None for a caller outside the real
        # pipeline (e.g. a test) -- that case just always uses the
        # normal floor, no escalation possible without a real job root.
        def _admission_floor_mb() -> int:
            if global_log_file is None:
                return GETORGANELLE_ADMISSION_FLOOR_MB
            job_root = job_root_from_global_log(global_log_file)
            return read_reserved_mb(sample.sample_id, oom_state_path(job_root), GETORGANELLE_ADMISSION_FLOOR_MB)

        # Real cgroup oom_kill count as of BEFORE this dispatch's very
        # first attempt even starts -- captured once, eagerly, here
        # (not inside the closure below, which would re-read it fresh
        # on every call and never be able to tell "already elevated
        # before we started" apart from "just increased because of
        # us"). None when undetectable (e.g. non-Linux, or no real
        # cgroup) -- the fallback below simply can't fire without it,
        # never guesses a baseline.
        _cgroup_path_at_start = current_cgroup_path()
        _oom_kill_baseline = read_cgroup_oom_kill_count(_cgroup_path_at_start) if _cgroup_path_at_start else None

        # Called by BaseTool.run() right after a real SIGKILL, before it
        # decides whether to requeue in place -- see oom_retry_decision's
        # own docstring on BaseTool for the full real design. The real
        # escalated reservation that triggers a retry (or not) is
        # written by launcher.oom_monitor's own record_oom_kill call,
        # which happens in a separate process (the monitor), not this
        # one -- so this polls for that real write to actually land
        # (the monitor's own kill-and-record is one synchronous step,
        # triggered by its own 5s poll cycle; this is a real wait for a
        # real write, not a guessed delay) before reading the result.
        def _oom_retry_decision() -> bool:
            if global_log_file is None:
                return False
            job_root = job_root_from_global_log(global_log_file)
            state_path = oom_state_path(job_root)
            deadline = time.monotonic() + _OOM_ESCALATION_RECORD_WAIT_SECONDS
            record = None
            while time.monotonic() < deadline:
                record = read_oom_escalation_record(sample.sample_id, state_path)
                if record is not None:
                    break
                time.sleep(1.0)
            if record is not None:
                if record.get("exceeds_machine_capacity"):
                    self.logger.warning(
                        f"{{{sample.sample_id}}} Killed while already running with the whole machine "
                        f"reserved (reserved_slots={record.get('reserved_slots')}) -- a proven, genuine "
                        f"resource limit, not retrying further."
                    )
                    return False
                return True

            # No record from our own monitor within the wait window --
            # real gap this closes, 2026-11-07: our own monitor polls
            # every few seconds, so a real spike fast enough to outrun
            # that poll can let the kernel's own cgroup OOM killer fire
            # FIRST, with nothing of ours ever recording a matching
            # escalation. Before failing closed, check the one real,
            # independent signal that still proves a genuine OOM kill
            # happened even when our own monitor never saw it coming:
            # this cgroup's own real oom_kill counter (see
            # read_cgroup_oom_kill_count's own docstring). That counter
            # is cgroup-wide, not scoped to this one sample, so this
            # only trusts a real INCREASE since this exact dispatch's
            # own start (captured once, above, before any attempt ever
            # ran) -- not an absolute count that could reflect some
            # earlier, unrelated sample's own kernel kill.
            current_count = read_cgroup_oom_kill_count(_cgroup_path_at_start) if _cgroup_path_at_start else None
            if (_oom_kill_baseline is not None and current_count is not None
                    and current_count > _oom_kill_baseline):
                ceiling_mb = cgroup_memory_max_mb()
                if ceiling_mb is not None:
                    # No real peak reading of our own exists for this
                    # kill (GNU time never got to write one, and our
                    # monitor never saw it either) -- the real cgroup
                    # ceiling itself is the honest, defensible proxy:
                    # the kernel only ever OOM-kills once a cgroup is
                    # genuinely AT that ceiling, so this sample's real
                    # peak was genuinely at least that large.
                    record_oom_kill(sample.sample_id, float(ceiling_mb), state_path, recommended_concurrent_slots())
                self.logger.warning(
                    f"{{{sample.sample_id}}} Killed by signal 9 with no record from our own monitor, but "
                    f"this cgroup's real oom_kill counter just increased ({_oom_kill_baseline} -> "
                    f"{current_count}) -- a genuine kernel-level OOM kill that beat our own proactive "
                    f"monitor's poll cycle. Escalating conservatively and retrying."
                )
                return True

            self.logger.warning(
                f"{{{sample.sample_id}}} Killed by signal 9, but no real OOM-monitor escalation record "
                f"appeared within {_OOM_ESCALATION_RECORD_WAIT_SECONDS:.0f}s, and the real cgroup "
                f"oom_kill counter shows no new kernel-level kill either -- no real basis to escalate "
                f"against, not retrying."
            )
            return False

        # Real gap closed, 2026-11-09: _admission_floor_mb above already
        # raises the real RAM floor this dispatch itself must see, once
        # escalated -- but that floor is just a point-in-time reading
        # for THIS caller; it was never actually enforced against OTHER
        # callers, so a small, unrelated sample could still see "enough
        # room for ME" at a slightly different moment and start
        # concurrently anyway (confirmed live the same night: a sample
        # escalated to the real machine ceiling was still killed with
        # an 8.4GB peak, nowhere near a genuine whole-machine need).
        # Once this sample's own real reservation has reached the
        # machine's real total capacity, it has earned a genuine solo
        # trial, not just a generous floor -- wait for (and respect)
        # any OTHER sample's own live trial, then claim this one before
        # actually invoking GetOrganelle.
        default_threads = int(threads)

        def _dispatch_admission_wait() -> None:
            lock_path = (exclusive_trial_lock_path(job_root_from_global_log(global_log_file))
                        if global_log_file is not None else None)
            wait_for_dispatch_admission(logger = self.logger, floor_mb = _admission_floor_mb(),
                                        sample_id = sample.sample_id, exclusive_lock_path = lock_path)
            if lock_path is None:
                return
            record = read_oom_escalation_record(sample.sample_id, oom_state_path(job_root_from_global_log(global_log_file)))
            if record is None:
                return
            reserved_slots = record.get("reserved_slots", 0)
            total_slots = recommended_concurrent_slots()
            # Real gap found live 2026-11-09/10: memory reservation
            # escalates gradually, one slot per kill, but threads used
            # to only ever jump once, at the very end, when reserved_slots
            # finally reached total_slots -- confirmed live, a sample
            # stuck at -t 4 through 8 real kills and several real hours
            # before finally earning that last-step jump. Threads now
            # track the same reserved_slots fraction, so the scale-up
            # happens across the climb, not just at its last step -- this
            # runs before build_command() on every attempt (including
            # in-place retries), so the real -t flag it emits reflects
            # whatever this sample has ACTUALLY earned so far.
            scaled_threads = thread_count_for_reserved_slots(reserved_slots, total_slots, default_threads)
            if scaled_threads > self.threads:
                self.threads = scaled_threads
                self.logger.info(f"{{{sample.sample_id}}} Escalated reservation ({reserved_slots}/{total_slots} "
                                 f"real slots): scaling real CPU threads to {self.threads}.")
            if reserved_slots >= total_slots:
                acquire_exclusive_trial(sample.sample_id, lock_path)

        super().__init__(tool_name = tool_name, working_dir = working_dir, logger = logger, metrics_dir = metrics_dir,
                         dispatch_admission_wait = _dispatch_admission_wait,
                         oom_retry_decision = _oom_retry_decision)
        self.output_dir = Path(output_dir)
        self.sample = sample
        self.organelle_type = organelle_type

        # Copying the caller's tool_options (never mutate what was passed
        # in, matching build_runner_from_seed_result's own stated care)
        # and forcing keep_temp_files on: without GetOrganelle's own
        # --keep-temp, the extended_*_paired.fq reads it recruits (the
        # ones actually handed to SPAdes for assembly) are deleted by
        # GetOrganelle itself before this wrapper ever sees them, so the
        # recruited-reads QC step below would have nothing to QC.
        self.tool_options = dict(tool_options or {})
        self.tool_options.setdefault("keep_temp_files", True)
        # Real GetOrganelle refuses outright against an "-o" directory
        # that already exists, no matter what's in it (confirmed live,
        # 2026-11-07: "...existed! Please use a new output directory, or
        # use '--continue'/'--overwrite'"). BaseTool.run()'s own retry
        # loop builds this command ONCE and reuses it verbatim for every
        # in-place retry (see run()'s own comment) -- so the SAME
        # output_dir, now genuinely non-empty from the killed attempt's
        # own real partial output, would hit this exact refusal on its
        # very first retry, for ANY caller (primary assembly or rescue,
        # reseeded or not), defeating the whole escalate-and-requeue
        # design before it ever gets a chance to use its extra room.
        # Harmless on a genuine first attempt (nothing to overwrite yet).
        self.tool_options.setdefault("overwrite", True)

        self.threads = int(threads)
        self.standardized_fasta = self.output_dir / "data.fasta"
        self.standardized_gfa = self.output_dir / "graph.gfa"
        # Resolved to an absolute path: FastQCRunner below is invoked
        # with this attempt's own output_dir as ITS working_dir (the
        # subprocess's cwd), so a relative recruited_reads_fastqc_dir
        # string passed as --outdir gets resolved a second time
        # against that already-nested cwd -- confirmed directly as the
        # real cause of every recruited-reads FastQC run failing in
        # 18ms with exit code 2 ("outdir does not exist"), silently,
        # on every real pipeline run so far.
        self.recruited_reads_fastqc_dir = (self.output_dir / "recruited_reads_fastqc").resolve()
        self.recruited_reads_qc_record_path = self.output_dir / f"{sample.sample_id}.recruited_reads_qc.json"

        # The classification discover_assembly_artifacts() produced,
        # captured in postprocess_outputs() before it cleans up
        # GetOrganelle's native output filenames. A caller that needs
        # the circular/fragmented/missing status after run() returns
        # must read this rather than re-deriving it, since by then the
        # native *.path_sequence.fasta/*.selected_graph.gfa files
        # discover_assembly_artifacts() looks for no longer exist.
        self.assembly_status: AssemblyArtifacts | None = None

    def validate_inputs(self) -> None:
        """
        Validates the organelle type, input files, and thread count for
        GetOrganelle. This function appropriately logs the result, exits
        otherwise.

        Args:
            None

        Returns:
            None
        """
        # Adding logger context.
        context = self._log_context()

        # Checking that the requested organelle type is supported.
        if self.organelle_type not in VALID_ORGANELLE_TYPES:
            self.logger.error(f"{context} Unsupported organelle type: {self.organelle_type}.")
            sys.exit(1)

        # Checking if input files exist.
        for path in (self.sample.r1, self.sample.r2):
            if not path.is_file():
                self.logger.error(f"{context} Input file does not exist: {path}.")
                sys.exit(1)

        # Checking that the requested thread count is valid for this machine.
        if self.threads <= 0:
            self.logger.error(f"{context} Please provide a valid number of threads (passed = {self.threads}).")
            sys.exit(1)
        if os.cpu_count() and self.threads > os.cpu_count():
            self.logger.error(f"{context} Number of threads is greater than number of cores. Please reduce the number of threads to {os.cpu_count()}.")
            sys.exit(1)

    def version_command(self) -> list[str] | None:
        """Real, confirmed flag: `get_organelle_from_reads.py --version` prints e.g. "GetOrganelle v1.7.7.1"."""
        return ["get_organelle_from_reads.py", "--version"]

    def build_command(self) -> list[str]:
        """
        Constructs the command to run GetOrganelle, seeded or unseeded
        depending on tool_options.

        Args:
            None

        Returns:
            list[str]: The command to run GetOrganelle.
        """
        # NOT pre-creating self.output_dir here: get_organelle_from_reads.py
        # refuses to run if its -o directory already exists ("Please use a
        # new output directory, or use '--continue'/'--overwrite'"), and it
        # creates that exact leaf directory itself on a successful run.
        # Pre-creating it here made every real invocation of this runner
        # fail outright.
        #
        # Its own directory creation is a bare os.mkdir(), NOT recursive --
        # confirmed directly against a real failure (FileNotFoundError on
        # a two-level-deep output_dir, e.g. a rescue attempt's
        # "assembly/attempt0" where "assembly/" didn't exist yet; the
        # single-level paths every prior real run happened to use never
        # exposed this). So the immediate PARENT must already exist, even
        # though the leaf itself must not -- hence parent.mkdir(), never
        # output_dir.mkdir().
        self.output_dir.parent.mkdir(parents = True, exist_ok = True)

        # Constructing the command.
        command = [
            "get_organelle_from_reads.py",
            "-1",
            str(self.sample.r1),
            "-2",
            str(self.sample.r2),
            "-o",
            str(self.output_dir),
            "-F",
            self.organelle_type,
            "-t",
            str(self.threads),
        ]

        # Returning the command with additional options (including -s, the
        # seed file, if a seed was selected for this sample).
        return self._add_options(command)

    def _recruited_reads_paths(self) -> tuple[Path, Path]:
        """
        The paired reads GetOrganelle recruited as mitochondrial-like
        during its own extending process -- the exact reads handed to
        SPAdes for assembly, before assembly. Confirmed against the
        installed get_organelle_from_reads.py: written as
        "{prefix}extended_1_paired.fq" / "{prefix}extended_2_paired.fq"
        in the output directory, and never removed by GetOrganelle
        itself (unlike the pre-split "{prefix}extended_1/2.fq", which
        it deletes once keep_temp_files is not set).

        Args:
            None

        Returns:
            tuple[Path, Path]: The expected (r1, r2) recruited-reads paths, whether or not they actually exist yet.
        """
        prefix = str(self.tool_options.get("prefix") or "")
        return (
            self.output_dir / f"{prefix}extended_1_paired.fq",
            self.output_dir / f"{prefix}extended_2_paired.fq",
        )

    def _run_recruited_reads_qc(self) -> None:
        """
        Runs FastQC on the reads GetOrganelle recruited as
        mitochondrial-like before assembly -- a real, early look at
        whether what's about to be assembled actually looks like
        clean, usable sequence, not just trusting the final
        assembly's own quality.

        Always writes a RecruitedReadsQCRecord -- including
        the absent/empty/malformed cases -- so a caller never has to
        re-derive what happened from log lines alone; this is a
        supplementary QC artifact, not a correctness gate on the
        assembly itself, so none of these cases fail the whole run.
        Compression is never checked here: GetOrganelle's own
        extended_*_paired.fq files are always plain text (confirmed
        against the installed get_organelle_from_reads.py), never
        gzipped, so there is nothing to validate on that axis.

        Args:
            None

        Returns:
            None
        """
        context = self._log_context()
        r1, r2 = self._recruited_reads_paths()

        if not (r1.is_file() and r2.is_file()):
            self.logger.warning(f"{context} Recruited-reads QC skipped: expected {r1.name}/{r2.name} not found in {self.output_dir} (keep_temp_files may not have taken effect).")
            self._write_recruited_reads_qc_record(RecruitedReadsQCRecord(
                sample_id = self.sample.sample_id, source_output_dir = str(self.output_dir),
                r1_path = None, r2_path = None, r1_read_count = None, r2_read_count = None,
                paired = None, qc_ran = False, fastqc_dir = None,
                reason = "recruited-reads files not found",
            ))
            return

        if r1.stat().st_size == 0 or r2.stat().st_size == 0:
            self.logger.warning(f"{context} Recruited-reads QC skipped: {r1.name}/{r2.name} present but empty.")
            self._write_recruited_reads_qc_record(RecruitedReadsQCRecord(
                sample_id = self.sample.sample_id, source_output_dir = str(self.output_dir),
                r1_path = str(r1), r2_path = str(r2), r1_read_count = None, r2_read_count = None,
                paired = None, qc_ran = False, fastqc_dir = None,
                reason = "recruited-reads file(s) empty",
            ))
            return

        r1_count = _count_fastq_records(r1)
        r2_count = _count_fastq_records(r2)
        paired = None if (r1_count is None or r2_count is None) else (r1_count == r2_count)

        if r1_count is None or r2_count is None:
            self.logger.warning(f"{context} Recruited-reads {r1.name}/{r2.name} line count is not a multiple of 4 (malformed FASTQ); running QC anyway, read counts unavailable.")
        elif not paired:
            self.logger.warning(f"{context} Recruited-reads pairing mismatch: {r1.name} has {r1_count} reads, {r2.name} has {r2_count}.")

        recruited_sample = Sample(sample_id = self.sample.sample_id, r1 = r1, r2 = r2)
        qc_runner = FastQCRunner(
            working_dir = self.output_dir,
            output_dir = self.recruited_reads_fastqc_dir,
            sample = recruited_sample,
            r1_output_stem = "recruited_R1",
            r2_output_stem = "recruited_R2",
            logger = self.logger,
            threads = self.threads,
        )
        qc_ran = False
        reason = "qc completed"
        try:
            result = qc_runner.run()
            if result.success:
                qc_ran = True
            else:
                reason = f"fastqc exited {result.return_code}"
                self.logger.warning(f"{context} FastQC on recruited reads exited with return code {result.return_code}; continuing, this is supplementary QC, not an assembly correctness check.")
        except SystemExit:
            reason = "fastqc failed validation"
            self.logger.warning(f"{context} FastQC on recruited reads failed validation; continuing, this is supplementary QC, not an assembly correctness check.")

        self._write_recruited_reads_qc_record(RecruitedReadsQCRecord(
            sample_id = self.sample.sample_id, source_output_dir = str(self.output_dir),
            r1_path = str(r1), r2_path = str(r2), r1_read_count = r1_count, r2_read_count = r2_count,
            paired = paired, qc_ran = qc_ran,
            fastqc_dir = str(self.recruited_reads_fastqc_dir) if qc_ran else None,
            reason = reason,
        ))

    def _write_recruited_reads_qc_record(self, record: RecruitedReadsQCRecord) -> None:
        """
        Persists a RecruitedReadsQCRecord as JSON in this attempt's
        own output_dir -- already attempt-specific (each rescue
        attempt gets its own output_dir), so this never collides or
        gets overwritten across rescue retries.

        Args:
            record (RecruitedReadsQCRecord): The record to persist.

        Returns:
            None
        """
        with self.recruited_reads_qc_record_path.open("w", encoding = "utf-8") as handle:
            json.dump(asdict(record), handle, indent = 4)

    def postprocess_outputs(self) -> None:
        """
        Standardizes GetOrganelle's own versioned output filenames to a
        fixed data.fasta/graph.gfa pair, runs FastQC on the reads it
        recruited before assembly, then removes everything else
        GetOrganelle left behind in the output directory.

        Args:
            None

        Returns:
            None
        """
        # Obtaining the logging context.
        context = self._log_context()

        # Classifying the result BEFORE anything below renames/deletes
        # GetOrganelle's native output -- discover_assembly_artifacts()
        # globs for the same native filenames, so it must run while
        # they still exist. Callers read this back via self.assembly_status.
        self.assembly_status = discover_assembly_artifacts(self.output_dir, self.logger)

        # discover_assembly_artifacts() above already chose the single
        # best candidate -- preferring a ".complete." one, since
        # GetOrganelle can legitimately leave behind several equally
        # valid circular path variants through a repeat region (e.g.
        # neogobius_melanostomus: 5 complete repeat_pattern candidates
        # from one real successful assembly, confirmed live 2026-11-06).
        # Re-globbing here and demanding exactly one, as this used to,
        # discarded a genuine circular result the moment GetOrganelle
        # produced more than one valid candidate.
        if self.assembly_status.fasta_path is None:
            self.logger.error(f"{context} No GetOrganelle path_sequence FASTA found in {self.output_dir}.")
            sys.exit(1)
        if self.assembly_status.gfa_path is None:
            self.logger.error(f"{context} No GetOrganelle selected_graph GFA found in {self.output_dir}.")
            sys.exit(1)

        # Standardizing the native outputs to fixed filenames.
        shutil.copy2(self.assembly_status.fasta_path, self.standardized_fasta)
        shutil.copy2(self.assembly_status.gfa_path, self.standardized_gfa)

        # QC'ing the recruited reads BEFORE the cleanup below deletes
        # them -- they're GetOrganelle's own intermediate files, not
        # part of the standardized output, so they don't survive past
        # this point otherwise.
        self._run_recruited_reads_qc()

        # Removing everything else GetOrganelle left behind.
        retained_names = {self.standardized_fasta.name, self.standardized_gfa.name, "assembly.log",
                           self.recruited_reads_fastqc_dir.name, self.recruited_reads_qc_record_path.name}
        for child in self.output_dir.iterdir():
            if child.name in retained_names:
                continue
            if child.is_dir():
                shutil.rmtree(child)
            else:
                child.unlink()

    def cleanup_failed_output(self) -> None:
        """
        Removes everything GetOrganelle left behind after a FAILED
        attempt -- no path_sequence.fasta/selected_graph.gfa exist to
        standardize in this case, so nothing in output_dir is a real
        output worth keeping. postprocess_outputs() (above) only ever
        runs on a successful BaseTool.run() -- a failed or crashed
        attempt skips it entirely, so without this, GetOrganelle's own
        multi-gigabyte working directory (SPAdes scratch, extended
        read files, bowtie2 indices) is retained forever instead of
        being pruned the same way a successful attempt's is. Confirmed
        against real failed attempts on disk: 17-34GB each, vs ~2MB
        for an equivalent successful one.

        Callers (run_getorganelle_assembly.py / run_rescue_reassembly.py)
        call this from their own failure branches, after writing this
        attempt's status record -- get_org.log.txt and every small
        *.json record already on disk (assembly_status, recruited-
        reads QC, a rescue attempt's own reseed provenance) are kept
        for diagnostics; everything else -- all of it genuine GetOrganelle
        scratch, never a declared pipeline output -- goes.

        Args:
            None

        Returns:
            None
        """
        if not self.output_dir.is_dir():
            return
        for child in self.output_dir.iterdir():
            if child.name == "get_org.log.txt" or child.suffix == ".json":
                continue
            if child.is_dir():
                shutil.rmtree(child)
            else:
                child.unlink()

    def validate_outputs(self) -> None:
        """
        Validates the standardized outputs of GetOrganelle. This function
        appropriately logs the result, exits otherwise.

        Args:
            None

        Returns:
            None
        """
        # Obtaining the logging context.
        context = self._log_context()

        # Iterating through each of the expected outputs - if they do not
        # exist or are empty, then log and exit.
        for path in (self.standardized_fasta, self.standardized_gfa):
            if not path.is_file() or path.stat().st_size == 0:
                self.logger.error(f"{context} Missing or empty assembly output: {path}")
                sys.exit(1)

    def _add_options(self, command: list[str]) -> list[str]:
        """
        Helper function which appends optional flags to an original command.

        Args:
            command (list[str]): The command to modify.

        Returns:
            list[str]: The modified command.
        """
        # A dictionary of options that require a value - the key is the
        # config name and the value is the true flag passed to GetOrganelle.
        # "seed_file" (-s) is the one that distinguishes a seeded run from
        # an unseeded one; its absence is a valid, intentional configuration.
        value_options = {
            "seed_file": "-s",
            "max_rounds": "-R",
            "kmer": "-k",
            "word_size": "-w",
            "pre_grouping": "-P",
            "max_reads": "--max-reads",
            "target_coverage": "--reduce-reads-for-coverage",
            "min_read_length": "--min-read-len",
            "max_read_length": "--max-read-len",
            "expected_max_size": "--expected-max-size",
            "expected_min_size": "--expected-min-size",
            "genes_file": "--genes",
            "exclude_fasta": "--exclude",
            "prefix": "--prefix",
            "round_output_prefix": "--out-per-round",
            "blast_path": "--which-blast",
            "bandage_path": "--which-bandage",
            "spades_path": "--which-spades",
            "disentangle_df": "--disentangle-df",
            "disentangle_time_limit": "--disentangle-time-limit",
        }

        # Iterating through each of the value options and appending values
        # that are not None to the command.
        for key, flag in value_options.items():
            value = self.tool_options.get(key)
            if value is not None:
                command.extend([flag, str(value)])

        # A dictionary of options that require a boolean value - the key is
        # the config name and the value is the true flag passed to GetOrganelle.
        boolean_options = {
            "overwrite": "--overwrite",
            "continue_run": "--continue",
            "fast_mode": "--fast",
            "reverse_lsc": "--reverse-lsc",
            "no_slim": "--no-slim",
            # The real GetOrganelle flag is --keep-temp (confirmed against
            # the installed get_organelle_from_reads.py's own --help and
            # argparse dest="keep_temp_files"); --keep-temp-files isn't a
            # recognized flag at all and GetOrganelle would reject it.
            "keep_temp_files": "--keep-temp",
            "verbose": "--verbose",
        }

        # Iterating through each of the boolean options and appending flags
        # that are True to the command.
        for key, flag in boolean_options.items():
            if self.tool_options.get(key, False):
                command.append(flag)

        # Returning the modified command.
        return command

# The status an assembly-artifact discovery can report.
AssemblyDiscoveryStatus = str   # one of "circular", "fragmented", "missing"

@dataclass(frozen = True, slots = True)
class AssemblyArtifacts:
    """
    The outcome of discovering GetOrganelle's native assembly output
    in one sample's output directory.

    Attributes:
        status (str): One of "circular" (a complete, circularized assembly was found), "fragmented" (at least one candidate was found, but none is marked complete/circular), or "missing" (no candidate at all).
        fasta_path (Path | None): The path to the chosen native FASTA candidate, or None when status is "missing".
        gfa_path (Path | None): The path to the chosen native GFA candidate, or None when status is "missing".
        candidate_count (int): How many native FASTA candidates GetOrganelle actually produced.
    """
    status: str
    fasta_path: Path | None
    gfa_path: Path | None
    candidate_count: int

def quarantine_existing_output_dir(output_dir: Path, logger: Logger) -> None:
    """
    Moves a stale output_dir aside before a fresh dispatch, if (and
    only if) one already exists -- real GetOrganelle refuses to run
    at all against an "-o" directory that already exists, regardless
    of what's inside it (its own real message: "Please use a new
    output directory, or use '--continue'/'--overwrite'"). A prior
    FAILED attempt's own cleanup_failed_output() deliberately keeps
    get_org.log.txt and every small *.json record for diagnostics, but
    never removes output_dir itself -- fine for a one-shot terminal
    failure, but a real, confirmed blocker for a sample this pipeline
    now deliberately retries in place with an escalated reservation
    (see launcher.oom_monitor) -- without this, the retry would fail
    immediately regardless of how much memory it was given.

    Provenance is protected, not discarded: the stale directory is
    renamed to a timestamped sibling (<name>.killed_<timestamp>/),
    never deleted -- the prior attempt's own log and status record
    stay fully readable for real forensic review, just out of the way
    of the fresh attempt about to use the real path.

    Called at the very start of a fresh dispatch, before this attempt
    even validates its own inputs -- a caller that's about to run
    GetOrganelleRunner against output_dir should call this first,
    always (a no-op when output_dir doesn't exist yet, the normal case
    for a sample's first-ever attempt).

    Args:
        output_dir (Path): The path this fresh attempt is about to run GetOrganelle against.
        logger (Logger): The logger to use.

    Returns:
        None
    """
    if not output_dir.exists():
        return
    timestamp = datetime.now().strftime("%Y%m%dT%H%M%S%f")
    quarantine_path = output_dir.parent / f"{output_dir.name}.killed_{timestamp}"
    logger.warning(f"Stale output directory {output_dir} already exists from a prior non-terminal attempt -- "
                   f"moving it aside to {quarantine_path} (provenance preserved) before this fresh dispatch.")
    output_dir.rename(quarantine_path)

def discover_assembly_artifacts(output_dir: Path, logger: Logger) -> AssemblyArtifacts:
    """
    Discovers GetOrganelle's own native assembly output in output_dir,
    before postprocess_outputs standardizes it to a fixed filename
    pair. GetOrganelle names a completed, circularized assembly with
    ".complete." in its path_sequence filename; anything else it
    produces is a fragmented/draft candidate. Zero candidates means
    the assembly failed outright.

    Args:
        output_dir (Path): The directory GetOrganelle wrote its native output to.
        logger (Logger): The logger to use for this step of the pipeline.

    Returns:
        AssemblyArtifacts: The discovered status and, when available, the
        single best candidate pair (a "complete" one if any exists,
        otherwise the first candidate found).
    """
    # Discovering every native FASTA/GFA candidate GetOrganelle left behind.
    fasta_candidates = sorted(output_dir.glob("*.path_sequence.fasta"))
    gfa_candidates = sorted(output_dir.glob("*.selected_graph.gfa"))

    # Handling the case where GetOrganelle produced nothing at all.
    if not fasta_candidates:
        logger.warning(f"No native assembly FASTA found in {output_dir}; assembly likely failed.")
        return AssemblyArtifacts(status = "missing", fasta_path = None, gfa_path = None, candidate_count = 0)

    # Preferring a candidate GetOrganelle itself marked complete/circular.
    complete_candidates = [path for path in fasta_candidates if ".complete." in path.name]
    chosen_fasta = complete_candidates[0] if complete_candidates else fasta_candidates[0]
    status = "circular" if complete_candidates else "fragmented"

    # Pairing the chosen FASTA with its matching GFA, by shared prefix
    # (GetOrganelle names a candidate's FASTA and GFA with the same
    # stem up to ".path_sequence."/.selected_graph.").
    stem = chosen_fasta.name.split(".path_sequence.fasta")[0]
    matching_gfa = next((path for path in gfa_candidates if path.name.startswith(stem)), None)
    chosen_gfa = matching_gfa or (gfa_candidates[0] if gfa_candidates else None)

    if status == "fragmented":
        logger.warning(f"No complete/circular GetOrganelle candidate in {output_dir}; using the best fragmented candidate ({chosen_fasta.name}).")
    else:
        logger.info(f"Found a complete GetOrganelle assembly in {output_dir} ({chosen_fasta.name}).")

    return AssemblyArtifacts(status = status, fasta_path = chosen_fasta, gfa_path = chosen_gfa,
                              candidate_count = len(fasta_candidates))

def downgrade_unverified_circular(status: AssemblyArtifacts) -> str:
    """
    Returns status.status, except "circular" is downgraded to
    "fragmented" -- for a caller reporting a FAILED or interrupted
    run's leftover native output (discover_assembly_artifacts() run
    against an output_dir after an exception or a non-zero exit, not
    a clean success). "circular" carries a real downstream
    consequence (rescue permanently stops once a sample is classified
    circular -- see utils.rescue.trigger's _NEVER_RESCUED_STATUS) that
    an unverified leftover checkpoint file from a run that never
    actually finished should never be able to grant.

    Confirmed live: a primary assembly killed mid-run (an external
    restart's kill sweep) right after GetOrganelle had already written
    a ".complete." checkpoint file internally, but before the process
    could exit normally, left behind a status.json recording
    status="circular" with success=false -- permanently blocking
    every configured rescue attempt for a sample that never actually
    produced a trustworthy result. The leftover file isn't discarded
    (still reported as "fragmented", real candidate data worth
    keeping), just never allowed to grant "circular"'s special,
    structural guarantee from an unverified run.

    Args:
        status (AssemblyArtifacts): The classification discover_assembly_artifacts() returned.

    Returns:
        str: status.status, or "fragmented" if it was "circular".
    """
    return "fragmented" if status.status == "circular" else status.status

@dataclass(frozen = True, slots = True)
class RecruitedReadsQCRecord:
    """
    The outcome of recruited-reads QC for one GetOrganelle attempt
    -- the structured artifact contract a caller reads instead
    of re-deriving read counts/pairing/provenance from log lines.

    Attributes:
        sample_id (str): The sample this record is for.
        source_output_dir (str): The GetOrganelle attempt's own output_dir this record belongs to (already attempt-specific, so repeated rescue attempts never collide).
        r1_path (str | None): The recruited-reads R1 path, or None if it was never found.
        r2_path (str | None): The recruited-reads R2 path, or None if it was never found.
        r1_read_count (int | None): R1's record count, or None if the file was empty/malformed/absent.
        r2_read_count (int | None): R2's record count, or None if the file was empty/malformed/absent.
        paired (bool | None): Whether R1 and R2 have matching read counts, or None when either count is unavailable.
        qc_ran (bool): Whether FastQC actually completed successfully on these reads.
        fastqc_dir (str | None): The FastQC output directory, or None when qc_ran is False.
        reason (str): A short, human-readable explanation of this outcome.
    """
    sample_id: str
    source_output_dir: str
    r1_path: str | None
    r2_path: str | None
    r1_read_count: int | None
    r2_read_count: int | None
    paired: bool | None
    qc_ran: bool
    fastqc_dir: str | None
    reason: str

def _count_fastq_records(path: Path) -> int | None:
    """
    Counts FASTQ records in a plain-text FASTQ file (4 lines per
    record, the universal FASTQ convention). Returns None rather than
    a count when the file's line count isn't a multiple of 4 -- that
    signals a malformed/truncated FASTQ, not a valid record count.

    Args:
        path (Path): The FASTQ file to count.

    Returns:
        int | None: The number of records, or None if the file is malformed.
    """
    with path.open("r", encoding = "utf-8", errors = "replace") as handle:
        line_count = sum(1 for _ in handle)
    if line_count % 4 != 0:
        return None
    return line_count // 4

def build_runner_from_seed_result(seed_result: SeedSelectionResult,
                                   working_dir: Path,
                                   output_dir: Path,
                                   sample: Sample,
                                   organelle_type: str = "animal_mt",
                                   tool_options: dict | None = None,
                                   threads: int = 4,
                                   logger: Logger | None = None,
                                   metrics_dir: Path | None = None) -> GetOrganelleRunner:
    """
    Configures a GetOrganelleRunner from a seed-selection result: a
    resolved seed FASTA is wired in as the -s seed file; a fallback
    result (no seed found) is left unseeded, so GetOrganelle falls
    back to its own default database for organelle_type.

    Args:
        seed_result (SeedSelectionResult): The outcome of seed_selector.select_seed() for this sample.
        working_dir (Path): The working directory for the command.
        output_dir (Path): The directory to save the output files to.
        sample (Sample): The sample object that is being worked with.
        organelle_type (str): The GetOrganelle organelle type to assemble. Defaults to "animal_mt".
        tool_options (dict | None): A base dictionary of options to pass to the tool; a resolved seed is merged into a copy of this, never mutating the caller's dict.
        threads (int): The number of threads for this process to use.
        logger (Logger | None): The logger to use for this step of the pipeline.
        metrics_dir (Path | None): Where to persist this run's compute metrics as JSON; None skips persistence. See BaseTool.

    Returns:
        GetOrganelleRunner: A runner configured with the seed result wired in.
    """
    # Copying the base tool_options so the caller's own dict is never
    # mutated by this function.
    resolved_options = dict(tool_options or {})

    # Wiring in the resolved seed FASTA only for a genuine (non-fallback) result.
    if not seed_result.fallback_used and seed_result.fasta_path is not None:
        resolved_options["seed_file"] = str(seed_result.fasta_path)

    return GetOrganelleRunner(
        working_dir = working_dir,
        output_dir = output_dir,
        sample = sample,
        organelle_type = organelle_type,
        tool_options = resolved_options,
        threads = threads,
        logger = logger,
        metrics_dir = metrics_dir,
    )
