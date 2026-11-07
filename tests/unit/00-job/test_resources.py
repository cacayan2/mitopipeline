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

def test_job_root_from_global_log_strips_the_logs_00_job_suffix() -> None:
    global_log = Path("/home/fishl/mitogenome_project/outputs/final_output_v020/logs/00-job/final_output_v020.global.log")
    assert resources.job_root_from_global_log(global_log) == Path("/home/fishl/mitogenome_project/outputs/final_output_v020")

def test_oom_state_path_lives_inside_the_jobs_own_mlr_folder() -> None:
    job_root = Path("/home/fishl/mitogenome_project/outputs/final_output_v020")
    assert resources.oom_state_path(job_root) == job_root / "mrl_final_output_v020" / "oom_escalation_state.json"

def test_read_reserved_mb_returns_default_for_an_unseen_sample(tmp_path: Path) -> None:
    state_path = tmp_path / "oom_escalation_state.json"
    assert resources.read_reserved_mb("never_killed_sample", state_path, default_mb = 8192) == 8192

def test_read_reserved_mb_returns_default_when_state_file_missing(tmp_path: Path) -> None:
    state_path = tmp_path / "does_not_exist" / "oom_escalation_state.json"
    assert resources.read_reserved_mb("any_sample", state_path, default_mb = 8192) == 8192

def test_record_oom_kill_reserves_exactly_the_needed_slots_plus_one(tmp_path: Path) -> None:
    """A real peak of 20000MB with an 8192MB slot needs ceil(20000/8192)=3
    slots to cover it -- plus one slot of headroom = 4, never a blind
    doubling of whatever it had before."""
    state_path = tmp_path / "oom_escalation_state.json"
    record = resources.record_oom_kill("big_sample", peak_mb = 20000.0, state_path = state_path, total_slots = 10, slot_mb = 8192)
    assert record["reserved_slots"] == 4
    assert record["reserved_mb"] == 4 * 8192
    assert record["ran_alone"] is False
    assert record["exceeds_machine_capacity"] is False
    assert record["kills"] == 1

def test_record_oom_kill_persists_and_is_read_back_by_read_reserved_mb(tmp_path: Path) -> None:
    # peak 9000MB needs ceil(9000/8192)=2 slots just to cover it, plus
    # 1 slot of headroom = 3 total.
    state_path = tmp_path / "oom_escalation_state.json"
    resources.record_oom_kill("persisted_sample", peak_mb = 9000.0, state_path = state_path, total_slots = 10, slot_mb = 8192)
    assert resources.read_reserved_mb("persisted_sample", state_path, default_mb = 8192) == 3 * 8192

def test_record_oom_kill_caps_at_running_alone_when_requirement_exceeds_total_slots(tmp_path: Path) -> None:
    """A real peak so large that needed+1 slots would exceed the whole
    machine's real total_slots is capped at total_slots exactly (running
    alone), not an impossible number of slots -- and marks ran_alone so
    the next kill (if any) is recognized as a proven, final ceiling."""
    state_path = tmp_path / "oom_escalation_state.json"
    record = resources.record_oom_kill("huge_sample", peak_mb = 90000.0, state_path = state_path, total_slots = 10, slot_mb = 8192)
    assert record["reserved_slots"] == 10
    assert record["ran_alone"] is True
    assert record["exceeds_machine_capacity"] is False

def test_record_oom_kill_marks_exceeds_capacity_when_killed_again_after_running_alone(tmp_path: Path) -> None:
    """The real, definitive terminal case: already escalated to the
    whole machine and STILL killed -- proven to exceed real capacity,
    not escalated further (Emil's explicit call, 2026-11-07)."""
    state_path = tmp_path / "oom_escalation_state.json"
    resources.record_oom_kill("doomed_sample", peak_mb = 90000.0, state_path = state_path, total_slots = 10, slot_mb = 8192)
    record = resources.record_oom_kill("doomed_sample", peak_mb = 85000.0, state_path = state_path, total_slots = 10, slot_mb = 8192)
    assert record["exceeds_machine_capacity"] is True
    assert record["kills"] == 2

def test_record_oom_kill_never_shrinks_a_samples_own_prior_reservation(tmp_path: Path) -> None:
    """A second, smaller observed peak (e.g. a different real
    bottleneck this time) must never reduce an already-escalated
    reservation below what real evidence already proved it needs --
    reservations only ever grow."""
    state_path = tmp_path / "oom_escalation_state.json"
    resources.record_oom_kill("sample", peak_mb = 30000.0, state_path = state_path, total_slots = 10, slot_mb = 8192)
    first_reserved = resources.read_reserved_mb("sample", state_path, default_mb = 8192)
    record = resources.record_oom_kill("sample", peak_mb = 9000.0, state_path = state_path, total_slots = 10, slot_mb = 8192)
    assert record["reserved_mb"] >= first_reserved
