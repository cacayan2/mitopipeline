"""
resources.py

Real, live machine-resource detection for sizing GetOrganelle-class
jobs' memory use -- never a number hardcoded for one specific box.
Every reading comes from this machine, right now; a caller gets None
on anything undetectable (non-Linux, an unreadable /proc) and must
treat that as "unknown," never substitute a guess.

Why this exists: a GetOrganelle/SPAdes run that recruits off-target
reads can balloon past 30GB RSS on a single sample (confirmed live,
2026-11-06 -- the kernel OOM killer took out a 32GB `cyprinus_carpio`
process). Concurrency caps (`getorganelle_slots`) bound how many such
jobs run at once, but not how much any ONE of them is allowed to grow
to -- this module is the other half: a live, per-job memory ceiling
sized off whatever RAM this machine actually has free right now,
divided fairly across the jobs currently allowed to run.
"""

# Imports
from __future__ import annotations
import json
import os
import time
from logging import Logger
from pathlib import Path
from typing import Callable


# GetOrganelle's own documented calibration anchor (--help output,
# confirmed against the installed get_organelle_from_reads.py):
# "for personal computer with 8G memory, we suggest no more than 3E5"
# for -P/--pre-grouping. Used to scale -P proportionally to whatever
# this job's real memory budget turns out to be, instead of either
# hardcoding one machine's number or leaving it unset everywhere.
_GETORGANELLE_PRE_GROUPING_ANCHOR_MB = 8192
_GETORGANELLE_PRE_GROUPING_ANCHOR_VALUE = 300_000

# The minimum real headroom a new GetOrganelle-class job must see
# before it's allowed to start -- real measured p75 peak memory across
# 145 completed primary-assembly runs, 2026-11-06 (same anchor as
# _GETORGANELLE_PRE_GROUPING_ANCHOR_MB). NOT a reservation held on
# behalf of any specific already-running job: /proc/meminfo's own
# MemAvailable already, honestly, reflects every process's REAL
# current use, including a young job that genuinely hasn't allocated
# much yet -- there is nothing to correct for there. A real round of
# extension-stage growth plays out over minutes (confirmed live,
# 2026-11-06 -- completed runs ranged from several minutes to over an
# hour), while this check re-polls every few seconds (see
# wait_for_dispatch_admission's poll_seconds), so a job that starts
# growing shows up in the very next real reading with ample lead time
# to stop admitting further new jobs -- no guess about any specific
# job's future growth needed (Emil's own point, 2026-11-07: an earlier
# version of this reserved a flat placeholder for every young job
# regardless of whether it ever actually grew, taxing throughput for
# protection the live polling already provides for free).
GETORGANELLE_ADMISSION_FLOOR_MB = 8192

def read_meminfo_mb(key: str) -> int | None:
    """
    Reads one real value (in MB) out of /proc/meminfo -- e.g.
    "MemTotal" or "MemAvailable". Linux-only (both of this pipeline's
    real deployment targets -- fishlabworkstn and its WSL2 counterpart
    -- run a real Linux kernel); returns None rather than guessing on
    any platform or kernel where this file doesn't exist or doesn't
    parse, so a caller never mistakes "undetectable" for "zero."

    Args:
        key (str): The /proc/meminfo field name (e.g. "MemAvailable").

    Returns:
        int | None: The value in MB, or None if it couldn't be read.
    """
    try:
        with open("/proc/meminfo", "r", encoding = "utf-8") as handle:
            for line in handle:
                if line.startswith(f"{key}:"):
                    # Real /proc/meminfo lines look like "MemAvailable:   12345678 kB".
                    parts = line.split()
                    if len(parts) >= 2 and parts[1].isdigit():
                        return int(parts[1]) // 1024
    except OSError:
        return None
    return None

def read_total_ram_mb() -> int | None:
    """This machine's real total installed RAM, in MB, or None if undetectable."""
    return read_meminfo_mb("MemTotal")

def read_available_ram_mb() -> int | None:
    """
    This machine's real CURRENTLY available RAM, in MB -- the kernel's
    own "MemAvailable" estimate (free + reclaimable cache/buffers),
    not just "free", so page cache already in use by completed jobs
    doesn't masquerade as unavailable. None if undetectable.
    """
    return read_meminfo_mb("MemAvailable")

def recommended_concurrent_slots(min_mb_per_job: int = 8192, reserved_fraction: float = 0.15,
                                 max_slots: int | None = None) -> int:
    """
    A live, RAM-aware concurrency ceiling for GetOrganelle-class jobs
    (primary assembly, rescue mapping/reassembly) -- how many of them
    this machine can run at once without already being memory-starved
    on day one, before any single job even has a chance to grow.
    Never hardcoded to one machine's spec: reads this machine's real
    total RAM right now and divides by a configurable "comfortable
    minimum per job" anchor.

    min_mb_per_job's default (8192, 8GB) and reserved_fraction (0.15)
    deliberately match compute_job_memory_budget_mb's own defaults --
    confirmed live, 2026-11-06, that an EARLIER mismatch here (a lower
    6GB anchor with no reserved_fraction applied at all) let this
    function recommend more slots than the live per-job budget
    function could actually honor once it divided the same real RAM
    pool, producing a per-job ceiling (~5GB) below real GetOrganelle's
    own measured need (median 6.3GB, p75 8.4GB across 145 real runs) --
    a real, confirmed regression: two samples that assembled cleanly
    under the old, uncapped system failed with a genuine MemoryError
    under the new one. Keeping both functions' math consistent is what
    actually matters, not either number in isolation.

    Falls back to a conservative value of 1 when RAM is undetectable
    -- correctness over throughput when this machine's real headroom
    is simply unknown.

    Args:
        min_mb_per_job (int): How much RAM one job should comfortably be able to claim without the whole pool already being oversubscribed at idle. Defaults to 8192 (8GB) -- real measured p75 peak memory across 145 completed primary-assembly runs, 2026-11-06.
        reserved_fraction (float): The same safety margin compute_job_memory_budget_mb reserves for the OS/other processes -- kept identical here so the two functions agree on how much RAM is actually being divided up.
        max_slots (int | None): An optional hard ceiling (e.g. to leave headroom for other users/processes on a shared machine); None applies no extra ceiling beyond the machine's own real CPU count.

    Returns:
        int: The recommended slot count, at least 1.
    """
    total_ram_mb = read_total_ram_mb()
    cpu_count = os.cpu_count() or 1

    if total_ram_mb is None:
        slots = 1
    else:
        usable_mb = total_ram_mb * (1.0 - reserved_fraction)
        slots = max(1, int(usable_mb // max(1, min_mb_per_job)))

    slots = min(slots, cpu_count)
    if max_slots is not None:
        slots = min(slots, max_slots)
    return max(1, slots)

def compute_job_memory_budget_mb(concurrent_slots: int, reserved_fraction: float = 0.15) -> int | None:
    """
    This specific job's own fair-share memory ceiling, in MB -- read
    fresh, right now, off whatever RAM is actually free on this
    machine at this exact moment, not a number precomputed once at
    pipeline launch. Early in a run, with few jobs active, this
    naturally grants a generous budget (the machine's real headroom is
    used, not left idle); as more jobs pile up and real available RAM
    drops, newly-dispatched jobs correctly see a smaller, still-safe
    share -- the same live number every caller would see running the
    same check by hand.

    Args:
        concurrent_slots (int): How many GetOrganelle-class jobs are allowed to run at once (the same number sizing `getorganelle_slots`) -- this job's share is the current headroom divided across all of them.
        reserved_fraction (float): The fraction of currently-available RAM to leave unclaimed, as a safety margin for the OS and any other process on this machine (e.g. a shared box's other users). Defaults to 0.15 (15%).

    Returns:
        int | None: This job's memory budget in MB, or None if available RAM is undetectable -- callers must treat None as "no limit can be safely computed," never substitute a guessed number.
    """
    available_mb = read_available_ram_mb()
    if available_mb is None:
        return None
    usable_mb = available_mb * (1.0 - reserved_fraction)
    return max(1, int(usable_mb / max(1, concurrent_slots)))

def scale_pre_grouping(budget_mb: int | None) -> int | None:
    """
    Scales GetOrganelle's own -P/--pre-grouping value proportionally
    to a real computed memory budget, anchored to GetOrganelle's own
    documented guidance ("for personal computer with 8G memory, we
    suggest no more than 3E5"). This is a SOFT measure: it reduces how
    much pre-grouped read data GetOrganelle keeps in memory before a
    job has any chance to approach a hard ceiling, rather than relying
    on a hard kill as the only real defense.

    Args:
        budget_mb (int | None): This job's own computed memory budget, or None when undetectable.

    Returns:
        int | None: A scaled -P value, or None when budget_mb is None (meaning: don't override -P at all, defer to the user's own config or GetOrganelle's own default).
    """
    if budget_mb is None:
        return None
    scaled = _GETORGANELLE_PRE_GROUPING_ANCHOR_VALUE * budget_mb / _GETORGANELLE_PRE_GROUPING_ANCHOR_MB
    return max(10_000, int(scaled))

def dispatch_admission_ready(floor_mb: int = GETORGANELLE_ADMISSION_FLOOR_MB, reserved_fraction: float = 0.15) -> bool | None:
    """
    Whether a new GetOrganelle-class job may safely start RIGHT NOW --
    purely a live read of real current headroom, never a guess about
    any specific already-running job's future growth. See
    GETORGANELLE_ADMISSION_FLOOR_MB for why no such guess is needed.

    Args:
        floor_mb (int): The minimum real headroom (after reserved_fraction) a new job must see before starting. Defaults to GETORGANELLE_ADMISSION_FLOOR_MB.
        reserved_fraction (float): The same OS/safety margin every other budget function in this module reserves. Defaults to 0.15.

    Returns:
        bool | None: True if there's real room, False if not, None if available RAM is undetectable -- callers must treat None as "unknown," never substitute a guess (see wait_for_dispatch_admission, which fails open on None).
    """
    available_mb = read_available_ram_mb()
    if available_mb is None:
        return None
    usable_mb = available_mb * (1.0 - reserved_fraction)
    return usable_mb >= floor_mb

def wait_for_dispatch_admission(logger: Logger, floor_mb: int = GETORGANELLE_ADMISSION_FLOOR_MB,
                                reserved_fraction: float = 0.15, poll_seconds: float = 5.0,
                                sleep_fn: Callable[[float], None] = time.sleep) -> None:
    """
    Blocks (re-polling real system state every poll_seconds) until
    dispatch_admission_ready confirms it's genuinely safe for a new
    GetOrganelle-class job to start, or returns immediately if
    available RAM is undetectable (fails open, same as every other
    function in this module -- never blocks forever on an environment
    where the real check simply can't run, e.g. local macOS
    development).

    No hard cap on how many jobs may run concurrently -- throughput is
    bounded only by real, currently-measured headroom, never by an
    arbitrary slot count chosen in advance (Emil's explicit call,
    2026-11-07, after a flat getorganelle_slots ceiling sized for
    worst-case safety measurably starved real throughput on an earlier
    run that same night). poll_seconds is deliberately short (5s, not
    this module's usual 20s dynamic-ceiling cadence) so a newly-tight
    real condition is noticed promptly -- cheap to check this often
    against how slowly real GetOrganelle memory growth actually plays
    out (minutes, not seconds), so the short interval costs nothing in
    practice while keeping the real reaction window tight.

    Args:
        logger (Logger): Where to log once, the first time a call actually has to wait (never spammed every poll).
        floor_mb (int): See dispatch_admission_ready.
        reserved_fraction (float): See dispatch_admission_ready.
        poll_seconds (float): How often to re-check while waiting.
        sleep_fn (Callable[[float], None]): The sleep implementation to call between polls -- overridable for tests. Defaults to time.sleep.

    Returns:
        None
    """
    logged_wait = False
    while True:
        ready = dispatch_admission_ready(floor_mb, reserved_fraction)
        if ready is None or ready:
            return
        if not logged_wait:
            logger.info(f"Dispatch admission: waiting for real headroom before starting (need ~{floor_mb}MB free).")
            logged_wait = True
        sleep_fn(poll_seconds)

# How much real RAM is deliberately left OUTSIDE the pipeline's own
# cgroup, for the OS/kernel itself -- small, since the proactive
# monitor (not this margin) is the real, primary defense; this just
# keeps the box itself responsive even in the worst case.
_OOM_MONITOR_OS_RESERVE_MB = 3072

# The real gap kept between the monitor's own proactive kill
# threshold and the cgroup's hard memory.max -- wide enough that a
# real check every few seconds (see launcher.oom_monitor) always has
# a genuine chance to land a controlled, recorded kill BEFORE the
# kernel's own uncontrolled, unrecorded cgroup OOM kill would trigger
# instead. Emil's own explicit call, 2026-11-07: the ceiling should be
# as high as real safety allows, not conservatively low -- this is
# sized only to win the race against the kernel, not to leave RAM
# idle out of excess caution.
_OOM_MONITOR_SAFETY_GAP_MB = 3072

def read_cgroup_working_set_bytes(cgroup_path: Path) -> int | None:
    """
    This cgroup's real "working set" -- memory.current minus its own
    real inactive_file (the LRU-cold, first-to-be-reclaimed file
    cache) -- the same correction container runtimes/cAdvisor/
    Kubernetes apply to report a container's own real memory pressure,
    and conceptually the same one read_available_ram_mb already
    applies at the whole-machine level via /proc/meminfo's own
    MemAvailable. Raw memory.current includes reclaimable file cache
    (FASTQ/BLAST-database reads, in this pipeline's case) that the
    kernel evicts for free under real pressure, long before any OOM
    risk exists -- confirmed live, 2026-11-07: a fresh run's cgroup
    showed ~95GB memory.current while real anon+kernel usage was only
    ~5GB, nearly all of the rest reclaimable file cache. Treating raw
    memory.current as real pressure made the proactive monitor kill
    small, healthy GetOrganelle jobs that were never actually at risk.

    Args:
        cgroup_path (Path): The real cgroup directory to read (see current_cgroup_path).

    Returns:
        int | None: The real working-set size, in bytes, or None if memory.current or memory.stat can't be read.
    """
    try:
        current_bytes = int((cgroup_path / "memory.current").read_text().strip())
        inactive_file_bytes = 0
        with (cgroup_path / "memory.stat").open("r", encoding = "utf-8") as handle:
            for line in handle:
                key, _, value = line.strip().partition(" ")
                if key == "inactive_file":
                    inactive_file_bytes = int(value)
                    break
    except (OSError, ValueError):
        return None
    return max(0, current_bytes - inactive_file_bytes)

def cgroup_memory_max_mb() -> int | None:
    """
    The real memory.max to set on the pipeline's own containing
    cgroup -- this machine's real total RAM, minus only a small
    reserve for the OS itself (_OOM_MONITOR_OS_RESERVE_MB). None if
    real total RAM is undetectable.

    Returns:
        int | None: The real cgroup ceiling, in MB.
    """
    total_mb = read_total_ram_mb()
    if total_mb is None:
        return None
    return max(1024, total_mb - _OOM_MONITOR_OS_RESERVE_MB)

def oom_monitor_threshold_mb() -> int | None:
    """
    The real aggregate usage, in MB, above which launcher.oom_monitor
    proactively intervenes -- set below cgroup_memory_max_mb() by
    _OOM_MONITOR_SAFETY_GAP_MB, so a real, controlled, recorded kill
    always has room to land before the kernel's own cgroup-level OOM
    killer would otherwise trigger uncontrolled and unrecorded. None
    if real total RAM is undetectable.

    Returns:
        int | None: The real proactive threshold, in MB.
    """
    ceiling_mb = cgroup_memory_max_mb()
    if ceiling_mb is None:
        return None
    return max(1024, ceiling_mb - _OOM_MONITOR_SAFETY_GAP_MB)

def current_cgroup_path() -> Path | None:
    """
    This process's own real cgroup v2 directory, read fresh from
    /proc/self/cgroup -- lets a process placed inside a fresh
    `systemd-run --user --scope` cgroup discover that cgroup's real
    path itself, without the caller already having to know it. None
    on non-Linux, a cgroup v1-only system (no unified "0::" line), or
    any other unreadable/unparseable case -- never a guess.

    Returns:
        Path | None: The real absolute cgroup directory under /sys/fs/cgroup, or None if undetectable.
    """
    try:
        with open("/proc/self/cgroup", "r", encoding = "utf-8") as handle:
            for line in handle:
                parts = line.strip().split(":", 2)
                if len(parts) == 3 and parts[0] == "0" and parts[1] == "":
                    return Path("/sys/fs/cgroup") / parts[2].lstrip("/")
    except OSError:
        return None
    return None

def job_root_from_global_log(global_log_file: Path) -> Path:
    """
    Derives this job's own output root from its global log file path
    (always .../outputs/<job_id>/logs/00-job/<job_id>.global.log,
    already passed to every exec script as --global-log-file -- no new
    CLI plumbing needed). Used to find this job's own mlr_<job_id>
    scratch folder, one of the four places this pipeline may ever
    write a file on a shared machine.

    Args:
        global_log_file (Path): The job's global log file path, as received via --global-log-file.

    Returns:
        Path: This job's own output root directory.
    """
    return Path(global_log_file).parent.parent.parent

def oom_state_path(job_root: Path) -> Path:
    """
    Where this job's real OOM-kill/reservation state is persisted --
    inside the job's own mlr_<job_id> scratch folder (sanctioned write
    location), read by every GetOrganelleRunner dispatch and written
    by the proactive OOM monitor alone.

    Args:
        job_root (Path): This job's own output root (see job_root_from_global_log).

    Returns:
        Path: The state file's path.
    """
    return job_root / f"mrl_{job_root.name}" / "oom_escalation_state.json"

def read_reserved_mb(sample_id: str, state_path: Path, default_mb: int) -> int:
    """
    This sample's own real reserved memory share for its NEXT dispatch,
    based on what it has actually needed before -- default_mb (the
    normal, un-escalated floor) if this sample has never been killed,
    or its own real escalated reservation (see record_oom_kill) if it
    has. Never guesses; a sample that's never been seen gets exactly
    the same share as everyone else.

    Args:
        sample_id (str): The sample id.
        state_path (Path): See oom_state_path.
        default_mb (int): The normal, un-escalated reservation -- GETORGANELLE_ADMISSION_FLOOR_MB in production.

    Returns:
        int: The real reservation to request for this sample's next dispatch, in MB.
    """
    try:
        with state_path.open("r", encoding = "utf-8") as handle:
            state = json.load(handle)
    except (OSError, json.JSONDecodeError):
        return default_mb
    record = state.get(sample_id)
    if record is None or "reserved_mb" not in record:
        return default_mb
    return max(default_mb, int(record["reserved_mb"]))

def record_oom_kill(sample_id: str, peak_mb: float, state_path: Path, total_slots: int,
                    slot_mb: int = GETORGANELLE_ADMISSION_FLOOR_MB) -> dict:
    """
    Records a real observed OOM kill for sample_id -- called ONLY by
    the proactive OOM monitor, immediately after it kills a process,
    using that process's own real last-measured RSS (not a guess).

    Escalation is in whole slots, never a blind multiplier (Emil's own
    explicit call, 2026-11-07): the next reservation is exactly enough
    real slots to have covered the peak just observed, plus ONE extra
    slot of headroom -- nothing more. Once that requirement would
    already need the entire machine (>= total_slots), the next retry
    is capped at running completely alone, as the final escalation
    step. If a sample that already ran alone gets killed again, that
    is real, definitive proof this sample's need exceeds the whole
    machine's real capacity -- recorded as exceeds_machine_capacity.
    This is the cross-dispatch backstop: BaseTool.run()'s own
    max_oom_retries already gives a killed job a bounded number of
    real in-place retries first (same dispatch); this escalation only
    matters once that budget is exhausted and a later, separate
    dispatch (the next rescue attempt) needs its own real reservation
    sized off what was actually observed.

    Args:
        sample_id (str): The sample id.
        peak_mb (float): This sample's own real peak RSS (MB) at the moment it was killed.
        state_path (Path): See oom_state_path.
        total_slots (int): The real total slot count this machine supports right now (see recommended_concurrent_slots) -- defines what "running alone" means.
        slot_mb (int): One slot's worth of memory, in MB. Defaults to GETORGANELLE_ADMISSION_FLOOR_MB.

    Returns:
        dict: The new state record now on file for this sample.
    """
    state_path.parent.mkdir(parents = True, exist_ok = True)
    try:
        with state_path.open("r", encoding = "utf-8") as handle:
            state = json.load(handle)
    except (OSError, json.JSONDecodeError):
        state = {}
    existing = state.get(sample_id, {})

    if existing.get("ran_alone"):
        # Already escalated all the way to the whole machine and STILL
        # killed -- a real, proven ceiling, not escalated further.
        record = {**existing, "kills": existing.get("kills", 0) + 1,
                 "last_observed_peak_mb": peak_mb, "exceeds_machine_capacity": True}
    else:
        slots_needed = max(1, -(-int(peak_mb) // slot_mb))  # ceil division, real units
        newly_required_slots = min(total_slots, slots_needed + 1)
        # Reservations only ever grow -- a smaller observed peak on a
        # later kill (a different real bottleneck that round) must
        # never undo real evidence an earlier kill already established.
        reserved_slots = max(newly_required_slots, existing.get("reserved_slots", 0))
        ran_alone = reserved_slots >= total_slots
        record = {
            "reserved_slots": reserved_slots,
            "reserved_mb": reserved_slots * slot_mb,
            "kills": existing.get("kills", 0) + 1,
            "last_observed_peak_mb": peak_mb,
            "ran_alone": ran_alone,
            "exceeds_machine_capacity": False,
        }

    state[sample_id] = record
    tmp_path = state_path.with_suffix(".json.tmp")
    with tmp_path.open("w", encoding = "utf-8") as handle:
        json.dump(state, handle, indent = 4)
    tmp_path.replace(state_path)
    return record
