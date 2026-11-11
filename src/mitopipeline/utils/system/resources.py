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
                                 max_slots: int | None = None, cpu_count_override: int | None = None,
                                 total_ram_mb_override: int | None = None) -> int:
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
        cpu_count_override (int | None): This machine's real core count, when already known (e.g. the pipeline's own `resources.cores` config knob) -- avoids a second, possibly inconsistent os.cpu_count() read. None reads it live.
        total_ram_mb_override (int | None): This machine's real total RAM in MB, when already known (e.g. the pipeline's own `resources.ram_gb` config knob). None reads it live.

    Returns:
        int: The recommended slot count, at least 1.
    """
    total_ram_mb = total_ram_mb_override if total_ram_mb_override is not None else read_total_ram_mb()
    cpu_count = cpu_count_override if cpu_count_override is not None else (os.cpu_count() or 1)

    if total_ram_mb is None:
        slots = 1
    else:
        usable_mb = total_ram_mb * (1.0 - reserved_fraction)
        slots = max(1, int(usable_mb // max(1, min_mb_per_job)))

    slots = min(slots, cpu_count)
    if max_slots is not None:
        slots = min(slots, max_slots)
    return max(1, slots)

def thread_count_for_reserved_slots(reserved_slots: int, total_slots: int, default_threads: int,
                                    os_reserve_cores: int = 1, cpu_count_override: int | None = None) -> int:
    """
    How many real CPU threads a sample should actually be given, given
    how many of the machine's real memory slots it has earned through
    OOM escalation -- real gap found live 2026-11-09/10: memory
    reservation escalates gradually, one slot at a time, tracking each
    kill's own real observed peak, but threads used to only ever jump
    ONCE, at the very end, when a sample reached a genuine exclusive
    trial (all slots). Every kill before that final one left real idle
    cores on the table (confirmed live: a sample stuck at -t 4 through
    8 real kills and several real hours, on a 24-core machine, before
    finally reaching exclusivity). Threads now track the SAME
    reserved_slots fraction memory already does, so the scale-up
    happens progressively across the escalation climb instead of only
    at its very last step -- which still converges to exactly the
    prior all-or-nothing behavior once reserved_slots reaches
    total_slots.

    Args:
        reserved_slots (int): This sample's own current reserved_slots, from its real OOM escalation record (oom_escalation_state.json). 0 (no record yet) yields default_threads.
        total_slots (int): This machine's real total concurrency ceiling (recommended_concurrent_slots()) -- what reserved_slots is itself capped against.
        default_threads (int): This job's own ordinary, unescalated thread count (whatever the pipeline's config/CLI already specifies) -- the real floor; this function only ever scales UP from it, never down, so an unescalated job's thread count is always left exactly as configured.
        os_reserve_cores (int): How many real cores to leave for the OS/oom_monitor itself, same small-margin philosophy as this module's RAM-side reserved_fraction. Defaults to 1.
        cpu_count_override (int | None): This machine's real core count, when already known. None reads it live.

    Returns:
        int: At least default_threads, never more than this machine's real core count minus os_reserve_cores.
    """
    cpu_count = cpu_count_override if cpu_count_override is not None else (os.cpu_count() or 1)
    max_threads = max(1, cpu_count - os_reserve_cores)
    if total_slots <= 0 or reserved_slots <= 0:
        return default_threads
    proportional = -(-(max_threads * min(reserved_slots, total_slots)) // total_slots)  # ceil division, real units
    return max(default_threads, min(max_threads, proportional))

def recommended_thread_count(total_slots: int, reserved_slots: int = 0, os_reserve_cores: int = 1,
                             cpu_count_override: int | None = None) -> int:
    """
    This machine's own live, intelligently-decided thread count for
    one job in a pool of up to total_slots concurrently-running jobs
    -- replaces a manually hand-tuned, per-tool, per-machine `threads:`
    config value with a number derived purely from this machine's real
    core count and the SAME slot count already sizing memory (never a
    number someone has to notice is wrong and re-tune by hand, e.g.
    the real, documented history in prod_config.yaml of getorganelle/
    assembly_rescue's own threads being lowered by hand after observing
    real CPU utilization on one specific machine).

    An ordinary, unescalated job gets a fair SHARE of the real core
    budget (cores / total_slots) -- enough concurrent jobs can run at
    once to use every real core, without any one of them being able to
    starve the others. A job that has earned extra memory through real
    OOM escalation (reserved_slots > 0) gets a proportionally larger
    share too, via thread_count_for_reserved_slots -- see that
    function's own docstring for why that climbs gradually rather than
    jumping once at the end.

    Args:
        total_slots (int): How many jobs this pool is sized to run concurrently (recommended_concurrent_slots()) -- the same number this job's own fair share is divided against.
        reserved_slots (int): This job's own current reserved_slots from a real OOM escalation record, if any. 0 (the default) means "no escalation yet" -- just the ordinary fair share.
        os_reserve_cores (int): How many real cores to leave for the OS itself. Defaults to 1.
        cpu_count_override (int | None): This machine's real core count, when already known (e.g. the pipeline's own `resources.cores` config knob). None reads it live.

    Returns:
        int: At least 1, never more than this machine's real core count minus os_reserve_cores.
    """
    cpu_count = cpu_count_override if cpu_count_override is not None else (os.cpu_count() or 1)
    max_threads = max(1, cpu_count - os_reserve_cores)
    fair_share = max(1, max_threads // max(1, total_slots))
    if reserved_slots <= 0:
        return fair_share
    return thread_count_for_reserved_slots(reserved_slots, total_slots, fair_share,
                                           os_reserve_cores = os_reserve_cores, cpu_count_override = cpu_count)

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

def dispatch_admission_ready(floor_mb: int = GETORGANELLE_ADMISSION_FLOOR_MB) -> bool | None:
    """
    Whether a new GetOrganelle-class job may safely start RIGHT NOW --
    purely a live read of real current headroom, never a guess about
    any specific already-running job's future growth. See
    GETORGANELLE_ADMISSION_FLOOR_MB for why no such guess is needed.

    Real gap found live 2026-11-11: this used to apply ANOTHER 15%
    safety margin on top of floor_mb (available_mb * 0.85 >= floor_mb)
    -- fine for the ordinary flat 8GB floor, but floor_mb for a sample
    already escalated to the real machine ceiling (read_reserved_mb)
    is itself already ~85% of total RAM (that's exactly how
    recommended_concurrent_slots sized it). Re-applying the same
    margin here compounded into needing ~96% of total RAM free --
    confirmed live: a sample needing 81920MB (of 98.2GB total) never
    actually saw enough "available" RAM to admit, even with nothing
    else running. floor_mb already carries whatever margin it needs;
    comparing raw availability against it directly is correct for
    both the flat ordinary floor and an already-margined escalated
    one.

    Args:
        floor_mb (int): The minimum real headroom a new job must see before starting. Defaults to GETORGANELLE_ADMISSION_FLOOR_MB.

    Returns:
        bool | None: True if there's real room, False if not, None if available RAM is undetectable -- callers must treat None as "unknown," never substitute a guess (see wait_for_dispatch_admission, which fails open on None).
    """
    available_mb = read_available_ram_mb()
    if available_mb is None:
        return None
    return available_mb >= floor_mb

def wait_for_dispatch_admission(logger: Logger, floor_mb: int = GETORGANELLE_ADMISSION_FLOOR_MB,
                                poll_seconds: float = 5.0,
                                sleep_fn: Callable[[float], None] = time.sleep,
                                sample_id: str | None = None, exclusive_lock_path: Path | None = None) -> None:
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

    Also waits out any OTHER sample's genuine exclusive trial (see
    exclusive_trial_blocks_others) when sample_id/exclusive_lock_path
    are given -- a real gap found live 2026-11-09: a sample that
    escalated all the way to needing the whole machine only ever had
    its OWN admission floor raised to match. That floor is just a
    point-in-time RAM reading, never an enforced reservation against
    OTHER callers, so a small, unrelated job could still see "enough
    room for ME" a moment later and start concurrently -- the exact
    same oversubscription the escalation was meant to rule out, now
    one layer up. This has no separate release step to get wrong: the
    other side's check is a live PID-liveness read, so a finished or
    killed exclusive trial stops blocking automatically.

    Args:
        logger (Logger): Where to log once, the first time a call actually has to wait (never spammed every poll).
        floor_mb (int): See dispatch_admission_ready.
        poll_seconds (float): How often to re-check while waiting.
        sleep_fn (Callable[[float], None]): The sleep implementation to call between polls -- overridable for tests. Defaults to time.sleep.
        sample_id (str | None): This dispatch's own sample id, so it never blocks on its own exclusive trial. None skips the exclusive-trial check entirely (e.g. a caller that doesn't track one).
        exclusive_lock_path (Path | None): See exclusive_trial_blocks_others. None skips the check.

    Returns:
        None
    """
    logged_wait = False
    while True:
        blocked_by_other = (sample_id is not None and exclusive_lock_path is not None
                            and exclusive_trial_blocks_others(sample_id, exclusive_lock_path))
        ready = dispatch_admission_ready(floor_mb)
        if not blocked_by_other and (ready is None or ready):
            return
        if not logged_wait:
            reason = "another sample's own genuine solo trial is in progress" if blocked_by_other \
                else f"real headroom before starting (need ~{floor_mb}MB free)"
            logger.info(f"Dispatch admission: waiting for {reason}.")
            logged_wait = True
        sleep_fn(poll_seconds)

def exclusive_trial_lock_path(job_root: Path) -> Path:
    """
    Where a real "this sample has earned the whole machine, let it run
    genuinely alone" marker is persisted -- same sanctioned scratch
    folder as oom_state_path, which this is a companion to.

    Args:
        job_root (Path): This job's own output root (see job_root_from_global_log).

    Returns:
        Path: The lock file's path.
    """
    return job_root / f"mrl_{job_root.name}" / "oom_exclusive_trial.json"

def acquire_exclusive_trial(sample_id: str, lock_path: Path) -> None:
    """
    Records that sample_id (this real process, by its own real pid) is
    about to run its own genuine solo trial -- called only by a
    dispatch whose own escalation record already shows reserved_slots
    at the real machine ceiling (see GetOrganelleRunner), right before
    it actually invokes GetOrganelle. No lock is acquired if another
    sample's own trial is still genuinely alive (exclusive_trial_blocks_others
    already made this dispatch wait for that); overwrites a stale
    entry from a real process that has since exited.

    Deliberately has no matching "release" call: this process's own
    real pid IS the lock's proof of life, checked fresh by
    exclusive_trial_blocks_others every time, so a normal exit, an
    uncaught exception, or even a real SIGKILL all self-release the
    instant this process is actually gone -- nothing to clean up, and
    nothing that can leak into a permanent block.

    Args:
        sample_id (str): The real sample id about to run alone.
        lock_path (Path): See exclusive_trial_lock_path.

    Returns:
        None
    """
    # Real race found live 2026-11-11: this used a FIXED tmp filename
    # shared by every caller -- harmless if only ever one sample earns
    # this at a time, but two samples reaching the ceiling within
    # moments of each other both write to the SAME tmp path, and
    # whichever one's rename runs second finds it already consumed by
    # the first (FileNotFoundError on tmp_path.replace), crashing that
    # dispatch outright instead of a clean retry. Suffixing with this
    # process's own real pid makes every caller's tmp file genuinely
    # unique, so concurrent writers can never collide.
    lock_path.parent.mkdir(parents = True, exist_ok = True)
    tmp_path = lock_path.with_suffix(f".json.{os.getpid()}.tmp")
    with tmp_path.open("w", encoding = "utf-8") as handle:
        json.dump({"sample_id": sample_id, "pid": os.getpid()}, handle, indent = 4)
    tmp_path.replace(lock_path)

def exclusive_trial_blocks_others(calling_sample_id: str, lock_path: Path,
                                  is_pid_alive: Callable[[int], bool] | None = None) -> bool:
    """
    Whether some OTHER sample currently holds a genuinely live
    exclusive trial (see acquire_exclusive_trial) -- true real
    exclusion, not just a generous admission floor: this reads the
    real recorded pid and confirms via /proc that it is still actually
    alive right now, so a trial that already finished (success,
    non-OOM failure, or a real SIGKILL -- none of which run any
    Python cleanup code) is recognized as over automatically, the very
    next time anyone checks. Never blocks a sample on its own trial,
    and a missing/unreadable/stale lock file never blocks anyone.

    Args:
        calling_sample_id (str): The sample about to ask "am I blocked."
        lock_path (Path): See exclusive_trial_lock_path.
        is_pid_alive (Callable[[int], bool] | None): Overridable real-pid-liveness check, for tests -- defaults to a real /proc read (Linux only; always reports "not alive" where /proc doesn't exist, same fail-open convention as the rest of this module).

    Returns:
        bool: True only if a DIFFERENT sample's own trial is genuinely still running right now.
    """
    if is_pid_alive is None:
        is_pid_alive = lambda pid: Path(f"/proc/{pid}").exists()
    try:
        with lock_path.open("r", encoding = "utf-8") as handle:
            record = json.load(handle)
    except (OSError, json.JSONDecodeError):
        return False
    holder = record.get("sample_id")
    pid = record.get("pid")
    if holder is None or holder == calling_sample_id or not isinstance(pid, int):
        return False
    return is_pid_alive(pid)

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

def read_cgroup_oom_kill_count(cgroup_path: Path) -> int | None:
    """
    This cgroup's real, cumulative `oom_kill` counter from its own
    memory.events -- a standard cgroup v2 interface file the kernel
    itself increments for every real OOM kill in this cgroup, whether
    launcher.oom_monitor caused it or the kernel's own OOM killer acted
    on its own. Exists specifically to close a real race: our own
    monitor polls every few seconds, so a real spike fast enough to
    outrun that poll can let the kernel's own killer fire first, with
    nothing of ours ever recording a matching oom_escalation_state.json
    entry. This is the one real, independent signal that still proves
    a genuine OOM kill happened even when our own monitor never saw it
    coming -- see api.getorganelle's own oom_retry_decision, which
    compares this against a real baseline taken at dispatch start
    rather than trusting an absolute count alone (this counter is
    cgroup-wide, not scoped to any one sample).

    Args:
        cgroup_path (Path): The real cgroup directory to read (see current_cgroup_path).

    Returns:
        int | None: The real cumulative oom_kill count, or None if memory.events can't be read or parsed.
    """
    try:
        with open(cgroup_path / "memory.events", "r", encoding = "utf-8") as handle:
            for line in handle:
                key, _, value = line.strip().partition(" ")
                if key == "oom_kill":
                    return int(value)
    except (OSError, ValueError):
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

def read_oom_escalation_record(sample_id: str, state_path: Path) -> dict | None:
    """
    This sample's own full, real OOM-escalation record (reserved_slots,
    reserved_mb, kills, last_observed_peak_mb, ran_alone,
    exceeds_machine_capacity), as written by record_oom_kill after a
    real kill -- unlike read_reserved_mb (which only ever returns the
    mb figure, silently defaulting when unseen), this returns the
    whole record so a caller can check ran_alone/exceeds_machine_capacity
    directly. Used by GetOrganelleRunner's own in-place OOM-retry
    decision (see api.getorganelle) to tell "still escalating, keep
    retrying" apart from "already proven to exceed real machine
    capacity, stop."

    Args:
        sample_id (str): The sample id.
        state_path (Path): See oom_state_path.

    Returns:
        dict | None: The real record, or None if this sample has never been killed, or the state file can't be read.
    """
    try:
        with state_path.open("r", encoding = "utf-8") as handle:
            state = json.load(handle)
    except (OSError, json.JSONDecodeError):
        return None
    return state.get(sample_id)

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
    slot of headroom -- nothing more, UNLESS that would mean a real
    kill made zero real progress (see below). Once that requirement
    would already need the entire machine (>= total_slots), the next
    retry is capped at running completely alone, as the final
    escalation step. If a sample that already ran alone gets killed
    again, that is real, definitive proof this sample's need exceeds
    the whole machine's real capacity -- recorded as
    exceeds_machine_capacity, which is exactly what BaseTool.run()'s
    own oom_retry_decision checks (see api.getorganelle) to finally
    stop requeuing this sample in place, same dispatch, no fixed retry
    count. The same real record also sizes a later, separate
    dispatch's own initial floor (e.g. the sample's next rescue
    attempt) via read_reserved_mb.

    Real gap found and fixed, 2026-11-09: admission has no hard
    concurrency cap (Emil's own explicit call, 2026-11-07, after a
    fixed ceiling measurably starved throughput) -- it only checks
    real free RAM at the moment a NEW job starts, never reserving room
    against jobs already escalated. That means a sample with a
    genuinely modest own peak (nowhere near a full-machine need) can
    still be the real largest consumer at the instant SYSTEM-WIDE
    pressure crosses the threshold, purely because many OTHER
    concurrently-running jobs are collectively oversubscribing the
    same pool -- confirmed live the same night: six real samples each
    killed 11-31 times, every single kill's own peak (4.8-13.6GB)
    nowhere near escalating reserved_slots toward total_slots, so
    ran_alone/exceeds_machine_capacity never had a chance to fire.
    Peak-only escalation has no way to tell "this sample is too big"
    apart from "the pool is just crowded right now" -- it can retry
    forever without ever reaching a real answer either way. Fixed by
    also escalating by at least one real slot per kill regardless of
    peak (still capped at total_slots, still never shrinking), so
    every real kill is guaranteed forward progress toward a genuine
    solo trial -- a repeatedly-killed sample now reaches
    ran_alone/exceeds_machine_capacity within at most total_slots real
    kills, whatever its own peak happens to be.

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
        prior_slots = existing.get("reserved_slots", 0)
        # At least one real slot of forward progress every real kill,
        # regardless of this kill's own peak -- a modest, genuinely
        # unlucky-timing peak must never leave a repeatedly-killed
        # sample stuck forever short of a genuine solo trial (see this
        # function's own docstring, "Real gap found and fixed").
        newly_required_slots = min(total_slots, max(slots_needed + 1, prior_slots + 1))
        # Reservations only ever grow -- a smaller observed peak on a
        # later kill (a different real bottleneck that round) must
        # never undo real evidence an earlier kill already established.
        reserved_slots = max(newly_required_slots, prior_slots)
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
