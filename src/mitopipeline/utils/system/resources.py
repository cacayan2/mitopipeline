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
import os
from pathlib import Path

# GetOrganelle's own documented calibration anchor (--help output,
# confirmed against the installed get_organelle_from_reads.py):
# "for personal computer with 8G memory, we suggest no more than 3E5"
# for -P/--pre-grouping. Used to scale -P proportionally to whatever
# this job's real memory budget turns out to be, instead of either
# hardcoding one machine's number or leaving it unset everywhere.
_GETORGANELLE_PRE_GROUPING_ANCHOR_MB = 8192
_GETORGANELLE_PRE_GROUPING_ANCHOR_VALUE = 300_000

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

def count_active_getorganelle_processes() -> int:
    """
    Real current count of running get_organelle_from_reads.py processes
    on this machine, via a direct /proc scan -- the REAL number of
    GetOrganelle-class jobs active right now, not the configured
    maximum slot count. Early in a run, or whenever fewer than the max
    are actually running, a job's live fair share (see
    compute_live_fair_share_mb) should reflect that real headroom
    rather than a number sized for a fully-loaded machine that isn't
    the current reality.

    Args:
        None

    Returns:
        int: The real count, 0 if undetectable (non-Linux, no /proc) or none are running.
    """
    count = 0
    proc_dir = Path("/proc")
    if not proc_dir.is_dir():
        return 0
    for entry in proc_dir.iterdir():
        if not entry.name.isdigit():
            continue
        try:
            cmdline = (entry / "cmdline").read_bytes()
        except OSError:
            continue
        if b"get_organelle_from_reads" in cmdline:
            count += 1
    return count

def compute_live_fair_share_mb(reserved_fraction: float = 0.15, min_divisor: int = 1) -> int | None:
    """
    This moment's real fair-share memory budget for ONE GetOrganelle
    -class job -- real currently-available RAM, divided across however
    many such jobs are ACTUALLY running right now (count_active_
    getorganelle_processes), not a fixed configured slot count.
    Recomputed fresh on every call, so a job watching this live (see
    BaseTool's dynamic_memory_limit_mb) sees its own fair share grow as
    sibling jobs finish and free up real RAM, and shrink if new ones
    start -- the live equivalent of "if this ran alone, it would get
    the whole machine," generalized to however many are genuinely
    sharing the pool at this exact instant, rather than a number fixed
    once at launch and never revisited (Emil's own explicit point,
    2026-11-06: a statically-divided ceiling can cut off a job that
    would have reached a real, natural conclusion if given the room
    that was genuinely, actually free at the time).

    Args:
        reserved_fraction (float): The fraction of currently-available RAM to leave unclaimed, as a safety margin for the OS and any other process on this machine. Defaults to 0.15 (15%), matching compute_job_memory_budget_mb.
        min_divisor (int): Never divide by fewer than this many jobs, even if the real count reads as 0 (e.g. a race right as this job itself started, before it shows up in its own /proc scan) -- defaults to 1, meaning "assume at least yourself."

    Returns:
        int | None: This moment's real fair-share budget in MB, or None if available RAM is undetectable.
    """
    available_mb = read_available_ram_mb()
    if available_mb is None:
        return None
    active = max(min_divisor, count_active_getorganelle_processes())
    usable_mb = available_mb * (1.0 - reserved_fraction)
    return max(1, int(usable_mb / active))

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
