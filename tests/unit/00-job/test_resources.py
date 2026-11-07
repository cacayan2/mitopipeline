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
    assert resources.recommended_concurrent_slots(min_mb_per_job = 6144) == 8

def test_recommended_concurrent_slots_never_exceeds_cpu_count(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(resources, "read_total_ram_mb", lambda: 1_000_000)   # ~977GB, plenty
    monkeypatch.setattr(resources.os, "cpu_count", lambda: 4)
    assert resources.recommended_concurrent_slots(min_mb_per_job = 6144) == 4

def test_recommended_concurrent_slots_respects_an_explicit_max(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(resources, "read_total_ram_mb", lambda: 1_000_000)
    monkeypatch.setattr(resources.os, "cpu_count", lambda: 64)
    assert resources.recommended_concurrent_slots(min_mb_per_job = 6144, max_slots = 3) == 3

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

def test_scale_pre_grouping_matches_getorganelles_own_documented_anchor() -> None:
    """GetOrganelle's own --help: 'for personal computer with 8G memory, we suggest no more than 3E5' for -P."""
    assert resources.scale_pre_grouping(8192) == 300_000

def test_scale_pre_grouping_scales_proportionally() -> None:
    assert resources.scale_pre_grouping(4096) == 150_000

def test_scale_pre_grouping_returns_none_when_budget_is_none() -> None:
    assert resources.scale_pre_grouping(None) is None

def test_scale_pre_grouping_never_goes_below_a_sane_floor() -> None:
    assert resources.scale_pre_grouping(1) == 10_000
