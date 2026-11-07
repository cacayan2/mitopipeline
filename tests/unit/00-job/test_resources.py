"""
test_resources.py

Unit tests for utils.system.resources -- live, portable
machine-resource detection used to size GetOrganelle-class jobs'
memory budgets. Every real /proc/meminfo read is monkeypatched to a
controlled value here; these tests confirm the arithmetic and the
honest-None-on-unknown contract, not any one real machine's numbers.
"""

# Imports
import builtins
from pathlib import Path
import pytest
from mitopipeline.utils.system import resources

_MEMINFO_TEXT = (
    "MemTotal:       16777216 kB\n"
    "MemFree:         1048576 kB\n"
    "MemAvailable:   10485760 kB\n"
)

class _FakeLogger:
    """A minimal stand-in -- these tests only ever call .info()."""
    def __init__(self) -> None:
        self.info_messages: list[str] = []
    def info(self, message: str) -> None:
        self.info_messages.append(message)

def test_read_meminfo_mb_parses_a_real_meminfo_line(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    meminfo = tmp_path / "meminfo"
    meminfo.write_text(_MEMINFO_TEXT)

    real_open = builtins.open
    def _fake_open(path, *args, **kwargs):
        if path == "/proc/meminfo":
            return real_open(meminfo, *args, **kwargs)
        return real_open(path, *args, **kwargs)
    monkeypatch.setattr(builtins, "open", _fake_open)

    assert resources.read_meminfo_mb("MemTotal") == 16384
    assert resources.read_meminfo_mb("MemAvailable") == 10240

def test_read_meminfo_mb_returns_none_when_unreadable(monkeypatch: pytest.MonkeyPatch) -> None:
    def _raise(*args, **kwargs):
        raise OSError("no such file")
    monkeypatch.setattr(builtins, "open", _raise)
    assert resources.read_meminfo_mb("MemTotal") is None

def test_read_meminfo_mb_returns_none_for_a_missing_key(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    meminfo = tmp_path / "meminfo"
    meminfo.write_text(_MEMINFO_TEXT)
    real_open = builtins.open
    monkeypatch.setattr(builtins, "open", lambda path, *a, **k: real_open(meminfo, *a, **k) if path == "/proc/meminfo" else real_open(path, *a, **k))
    assert resources.read_meminfo_mb("SwapTotal") is None

def test_recommended_concurrent_slots_scales_with_detected_ram(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(resources, "read_total_ram_mb", lambda: 49152)   # 48GB
    monkeypatch.setattr(resources.os, "cpu_count", lambda: 64)
    assert resources.recommended_concurrent_slots(min_mb_per_job = 6144, reserved_fraction = 0.0) == 8

def test_recommended_concurrent_slots_applies_the_same_reserved_fraction_as_the_budget_function(monkeypatch: pytest.MonkeyPatch) -> None:
    """Confirmed live, 2026-11-06: a mismatch here (slots sized off raw
    total RAM while the live per-job budget divided RESERVED RAM) let
    this function recommend more slots than the budget function could
    actually honor, producing a per-job ceiling below what GetOrganelle
    actually needed and causing a real regression. Both must agree."""
    monkeypatch.setattr(resources, "read_total_ram_mb", lambda: 49152)   # 48GB
    monkeypatch.setattr(resources.os, "cpu_count", lambda: 64)
    # 49152 * 0.85 = 41779.2; // 6144 = 6
    assert resources.recommended_concurrent_slots(min_mb_per_job = 6144, reserved_fraction = 0.15) == 6

def test_recommended_concurrent_slots_never_exceeds_cpu_count(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(resources, "read_total_ram_mb", lambda: 1_000_000)   # ~977GB, plenty
    monkeypatch.setattr(resources.os, "cpu_count", lambda: 4)
    assert resources.recommended_concurrent_slots(min_mb_per_job = 6144, reserved_fraction = 0.0) == 4

def test_recommended_concurrent_slots_respects_an_explicit_max(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(resources, "read_total_ram_mb", lambda: 1_000_000)
    monkeypatch.setattr(resources.os, "cpu_count", lambda: 64)
    assert resources.recommended_concurrent_slots(min_mb_per_job = 6144, reserved_fraction = 0.0, max_slots = 3) == 3

def test_recommended_concurrent_slots_falls_back_to_one_when_ram_undetectable(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(resources, "read_total_ram_mb", lambda: None)
    monkeypatch.setattr(resources.os, "cpu_count", lambda: 64)
    assert resources.recommended_concurrent_slots() == 1

def test_compute_job_memory_budget_mb_divides_available_ram_across_slots(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(resources, "read_available_ram_mb", lambda: 10000)
    # 10000 * 0.85 / 5 = 1700
    assert resources.compute_job_memory_budget_mb(concurrent_slots = 5, reserved_fraction = 0.15) == 1700

def test_compute_job_memory_budget_mb_returns_none_when_ram_undetectable(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(resources, "read_available_ram_mb", lambda: None)
    assert resources.compute_job_memory_budget_mb(concurrent_slots = 8) is None

def test_count_active_getorganelle_processes_counts_real_matching_cmdlines(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    proc_dir = tmp_path / "proc"
    for pid, cmd in [(100, b"get_organelle_from_reads.py\x00-1\x00r1.fq"),
                     (200, b"get_organelle_from_reads.py\x00-1\x00r2.fq"),
                     (300, b"fastqc\x00--threads\x002"),
                     (400, b"")]:
        pid_dir = proc_dir / str(pid)
        pid_dir.mkdir(parents = True)
        (pid_dir / "cmdline").write_bytes(cmd)
    (proc_dir / "not_a_pid").mkdir()

    monkeypatch.setattr(resources, "Path", lambda p: proc_dir if p == "/proc" else Path(p))
    assert resources.count_active_getorganelle_processes() == 2

def test_count_active_getorganelle_processes_returns_zero_without_proc(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(resources, "Path", lambda p: Path("/definitely/does/not/exist") if p == "/proc" else Path(p))
    assert resources.count_active_getorganelle_processes() == 0

def test_compute_live_fair_share_mb_divides_by_the_real_active_count(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(resources, "read_available_ram_mb", lambda: 10000)
    monkeypatch.setattr(resources, "count_active_getorganelle_processes", lambda: 4)
    # 10000 * 0.85 / 4 = 2125
    assert resources.compute_live_fair_share_mb() == 2125

def test_compute_live_fair_share_mb_grows_as_fewer_jobs_are_actually_running(monkeypatch: pytest.MonkeyPatch) -> None:
    """The whole point of this function: the SAME real available RAM
    divided across fewer real active jobs gives each one a bigger
    live share -- confirming a job's own ceiling can genuinely grow as
    siblings finish, not just stay fixed at an initial division."""
    monkeypatch.setattr(resources, "read_available_ram_mb", lambda: 10000)
    monkeypatch.setattr(resources, "count_active_getorganelle_processes", lambda: 4)
    busy_share = resources.compute_live_fair_share_mb()
    monkeypatch.setattr(resources, "count_active_getorganelle_processes", lambda: 1)
    quiet_share = resources.compute_live_fair_share_mb()
    assert quiet_share > busy_share

def test_compute_live_fair_share_mb_never_divides_by_fewer_than_min_divisor(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(resources, "read_available_ram_mb", lambda: 10000)
    monkeypatch.setattr(resources, "count_active_getorganelle_processes", lambda: 0)
    # Even reading 0 active processes (e.g. a race right at this job's
    # own startup) must still divide by at least 1, never explode.
    assert resources.compute_live_fair_share_mb(min_divisor = 1) == 8500

def test_compute_live_fair_share_mb_returns_none_when_ram_undetectable(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(resources, "read_available_ram_mb", lambda: None)
    assert resources.compute_live_fair_share_mb() is None

def test_scale_pre_grouping_matches_getorganelles_own_documented_anchor() -> None:
    """GetOrganelle's own --help: 'for personal computer with 8G memory, we suggest no more than 3E5' for -P."""
    assert resources.scale_pre_grouping(8192) == 300_000

def test_scale_pre_grouping_scales_proportionally() -> None:
    assert resources.scale_pre_grouping(4096) == 150_000

def test_scale_pre_grouping_returns_none_when_budget_is_none() -> None:
    assert resources.scale_pre_grouping(None) is None

def test_scale_pre_grouping_never_goes_below_a_sane_floor() -> None:
    assert resources.scale_pre_grouping(1) == 10_000

def test_dispatch_admission_ready_true_when_real_headroom_covers_the_floor(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(resources, "read_available_ram_mb", lambda: 10000)
    # 10000 * 0.85 = 8500 >= floor_mb=8192
    assert resources.dispatch_admission_ready(floor_mb = 8192) is True

def test_dispatch_admission_ready_false_when_real_headroom_is_tight(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(resources, "read_available_ram_mb", lambda: 5000)
    # 5000 * 0.85 = 4250 < floor_mb=8192
    assert resources.dispatch_admission_ready(floor_mb = 8192) is False

def test_dispatch_admission_ready_returns_none_when_ram_undetectable(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(resources, "read_available_ram_mb", lambda: None)
    assert resources.dispatch_admission_ready() is None

def test_wait_for_dispatch_admission_returns_immediately_when_already_ready(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(resources, "read_available_ram_mb", lambda: 10000)
    calls = []
    resources.wait_for_dispatch_admission(logger = _FakeLogger(), floor_mb = 8192,
                                          sleep_fn = lambda seconds: calls.append(seconds))
    assert calls == []

def test_wait_for_dispatch_admission_polls_until_real_room_opens_up(monkeypatch: pytest.MonkeyPatch) -> None:
    """Simulates real headroom that's tight at first, then opens up --
    confirms this actually blocks and re-checks, rather than trusting
    a one-time snapshot."""
    readings = iter([5000, 5000, 10000])
    monkeypatch.setattr(resources, "read_available_ram_mb", lambda: next(readings))
    sleeps = []
    resources.wait_for_dispatch_admission(logger = _FakeLogger(), floor_mb = 8192,
                                          sleep_fn = lambda seconds: sleeps.append(seconds))
    assert len(sleeps) == 2

def test_wait_for_dispatch_admission_fails_open_when_ram_undetectable(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(resources, "read_available_ram_mb", lambda: None)
    calls = []
    resources.wait_for_dispatch_admission(logger = _FakeLogger(), sleep_fn = lambda seconds: calls.append(seconds))
    assert calls == []

def test_wait_for_dispatch_admission_logs_only_once_while_waiting(monkeypatch: pytest.MonkeyPatch) -> None:
    readings = iter([5000, 5000, 10000])
    monkeypatch.setattr(resources, "read_available_ram_mb", lambda: next(readings))
    logger = _FakeLogger()
    resources.wait_for_dispatch_admission(logger = logger, floor_mb = 8192, sleep_fn = lambda seconds: None)
    assert len(logger.info_messages) == 1
