"""
oom_monitor.py

A standalone, long-running process that watches this job's own
memory-capped cgroup for real usage approaching exhaustion, and
proactively kills the single worst-offending GetOrganelle-class
process BEFORE the kernel's own OOM killer would -- giving us the
chance to record that process's own real peak RSS first (the kernel's
OOM killer gives no such warning), and to choose exactly which
process to kill, scoped to this job's own cgroup only.

Run once per pipeline launch, independently of Snakemake itself (see
launcher/pipeline_launcher.py's own real launch command) -- read-only
with respect to Snakemake's own state; the only write this makes is
to this job's own oom_escalation_state.json (see utils.system.resources),
read back by GetOrganelleRunner on every future dispatch of a sample
this has ever had to kill.
"""

# Imports
from __future__ import annotations
import argparse
import logging
import os
import re
import signal
import time
from pathlib import Path
from mitopipeline.utils.system.resources import (
    oom_state_path,
    recommended_concurrent_slots,
    record_oom_kill,
)

# Checked this often -- cheap (a handful of /proc reads), and real
# memory growth plays out over minutes (confirmed live, 2026-11-06/07),
# so this is far more frequent than it needs to be to react in time.
_POLL_SECONDS = 5.0

# A sample's own real output directory appears verbatim in its
# get_organelle_from_reads.py invocation's own "-o" argument -- e.g.
# ".../04d-assembly/<sample>/primary" or
# ".../07-assembly_rescue/<sample>/attempt_N/reassembly" -- confirmed
# directly against real dispatched commands, not guessed. The sample
# id is always the path segment immediately after either stage
# folder, in both real shapes -- no need to anchor on what follows it.
# Capturing the sample id this way (rather than, say, --sample-id)
# works because it's read from the REAL memory-hogging process itself
# (get_organelle_from_reads.py), not an ancestor wrapper script.
_SAMPLE_ID_PATTERN = re.compile(r"/(?:04d-assembly|07-assembly_rescue)/([^/]+)/")


def _read_cgroup_memory_current(cgroup_path: Path) -> int | None:
    """Real current memory usage of this cgroup, in bytes, or None if unreadable."""
    try:
        return int((cgroup_path / "memory.current").read_text().strip())
    except (OSError, ValueError):
        return None


def _real_getorganelle_processes() -> list[tuple[int, int, str]]:
    """
    Every currently-running get_organelle_from_reads.py process on
    this machine, with its real RSS and the sample id parsed from its
    own real "-o" argument.

    Returns:
        list[tuple[int, int, str]]: (pid, rss_mb, sample_id) for each real process found; sample_id is None if it couldn't be parsed.
    """
    results = []
    proc_dir = Path("/proc")
    if not proc_dir.is_dir():
        return results
    for entry in proc_dir.iterdir():
        if not entry.name.isdigit():
            continue
        pid = int(entry.name)
        try:
            cmdline = (entry / "cmdline").read_bytes().decode("utf-8", errors = "replace")
        except OSError:
            continue
        if "get_organelle_from_reads" not in cmdline:
            continue
        match = _SAMPLE_ID_PATTERN.search(cmdline.replace("\x00", " "))
        sample_id = match.group(1) if match else None
        rss_mb = _read_rss_mb(pid)
        if rss_mb is not None:
            results.append((pid, rss_mb, sample_id))
    return results


def _read_rss_mb(pid: int) -> int | None:
    """This pid's real current resident memory, in MB, or None if unreadable (already exited, permission)."""
    try:
        with open(f"/proc/{pid}/status", "r", encoding = "utf-8") as handle:
            for line in handle:
                if line.startswith("VmRSS:"):
                    kb = int(line.split()[1])
                    return kb // 1024
    except (OSError, ValueError, IndexError):
        pass
    return None


def run_monitor(cgroup_path: Path, threshold_mb: int, job_root: Path, logger: logging.Logger) -> None:
    """
    Polls forever (see _POLL_SECONDS): once this cgroup's real current
    usage crosses threshold_mb, kills the single real
    get_organelle_from_reads.py process using the most memory right
    now, records its real peak RSS, and continues. Never raises --
    any single bad read is logged and skipped, not fatal to the whole
    monitor (this process is expected to run for the entire, possibly
    many-hour, life of the pipeline).

    Args:
        cgroup_path (Path): This job's real cgroup directory (see launcher's own systemd-run --user --scope invocation).
        threshold_mb (int): The real aggregate usage, in MB, above which this monitor proactively intervenes.
        job_root (Path): This job's own output root (see utils.system.resources.job_root_from_global_log).
        logger (logging.Logger): Where to log every real check-in and every real kill.

    Returns:
        None
    """
    state_path = oom_state_path(job_root)
    total_slots = recommended_concurrent_slots()
    logger.info(f"OOM monitor started: cgroup={cgroup_path}, threshold={threshold_mb}MB, total_slots={total_slots}.")

    while True:
        try:
            current_bytes = _read_cgroup_memory_current(cgroup_path)
            if current_bytes is None:
                logger.warning("OOM monitor: could not read real cgroup memory.current this cycle; skipping.")
                time.sleep(_POLL_SECONDS)
                continue
            current_mb = current_bytes // (1024 * 1024)

            if current_mb >= threshold_mb:
                processes = _real_getorganelle_processes()
                if not processes:
                    logger.warning(f"OOM monitor: real usage {current_mb}MB over threshold, but no real get_organelle_from_reads process found -- likely transient (SPAdes/bowtie2 child spike); skipping this cycle.")
                else:
                    pid, rss_mb, sample_id = max(processes, key = lambda item: item[1])
                    # Killing the whole process GROUP, not just this one
                    # pid: BaseTool.run() isolates the wrapped command
                    # (and every real descendant it forks -- SPAdes,
                    # bowtie2) into its own group via start_new_session=
                    # True, separate from the exec script's own group.
                    # A single-pid kill would leave any such descendant
                    # as a live orphan, free to keep writing into this
                    # sample's own output directory even after the exec
                    # script's own cleanup_failed_output() has already
                    # run -- confirmed as a real, recurring risk earlier
                    # tonight. The exec script's own process is in a
                    # DIFFERENT group and is never touched by this, so
                    # it survives to write a real killed_by_signal status
                    # and run its own cleanup, exactly as intended.
                    try:
                        pgid = os.getpgid(pid)
                    except ProcessLookupError:
                        logger.info(f"OOM monitor: pid {pid} already exited before the kill landed.")
                        pgid = None
                    if pgid is not None:
                        logger.warning(f"OOM monitor: real usage {current_mb}MB >= threshold {threshold_mb}MB. Killing process group {pgid} (sample={sample_id}, real RSS={rss_mb}MB) -- the real worst offender right now.")
                        try:
                            os.killpg(pgid, signal.SIGKILL)
                        except ProcessLookupError:
                            logger.info(f"OOM monitor: process group {pgid} already exited before the kill landed.")
                    if sample_id is not None:
                        new_record = record_oom_kill(sample_id, float(rss_mb), state_path, total_slots)
                        logger.warning(f"OOM monitor: recorded real kill for {sample_id}: {new_record}.")
                    else:
                        logger.error(f"OOM monitor: killed pid {pid} but could not parse its sample id from cmdline -- no reservation recorded; this sample will retry at the normal, un-escalated share next time.")
        except Exception as e:  # noqa: BLE001 -- this process must never die from one bad cycle
            logger.error(f"OOM monitor: unexpected error this cycle, continuing: {e}.")

        time.sleep(_POLL_SECONDS)


def main() -> int:
    parser = argparse.ArgumentParser(description = "Proactively watch this job's cgroup and intelligently kill-and-record real OOM risk before the kernel does.")
    parser.add_argument("--cgroup-path", type = Path, required = True, help = "This job's real cgroup directory.")
    parser.add_argument("--threshold-mb", type = int, required = True, help = "Real aggregate cgroup usage, in MB, above which to intervene.")
    parser.add_argument("--job-root", type = Path, required = True, help = "This job's own output root directory.")
    parser.add_argument("--log-file", type = Path, required = True, help = "Where to log.")
    args = parser.parse_args()

    logger = logging.getLogger("oom_monitor")
    logger.setLevel(logging.INFO)
    handler = logging.FileHandler(args.log_file)
    handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s"))
    logger.addHandler(handler)

    run_monitor(args.cgroup_path, args.threshold_mb, args.job_root, logger)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
