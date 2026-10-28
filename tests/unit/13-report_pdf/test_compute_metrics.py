"""
test_compute_metrics.py

Unit tests for aggregating the real per-invocation compute metrics
BaseTool already persists into a per-stage report summary, and for
collecting the real (non-identifying) execution-host specs.
"""

# Imports
import json
from pathlib import Path
from mitopipeline.utils.reporting.compute_metrics import collect_compute_summary, collect_machine_info

def _write_metric(path: Path, tool_name: str, runtime: float, user_cpu: float,
                   system_cpu: float, peak_mb: float) -> None:
    path.parent.mkdir(parents = True, exist_ok = True)
    path.write_text(json.dumps({
        "tool_name": tool_name, "runtime_seconds": runtime,
        "user_cpu_seconds": user_cpu, "system_cpu_seconds": system_cpu, "peak_memory_mb": peak_mb,
    }))

def test_collect_compute_summary_aggregates_by_stage(tmp_path: Path) -> None:
    metrics_dir = tmp_path / "metrics"
    _write_metric(metrics_dir / "08b-annotation" / "sample1" / "mitos2.sample1.1.metrics.json",
                  "mitos2", runtime = 5.0, user_cpu = 4.0, system_cpu = 1.0, peak_mb = 100.0)
    _write_metric(metrics_dir / "08b-annotation" / "sample1" / "mitos2.sample1.2.metrics.json",
                  "mitos2", runtime = 3.0, user_cpu = 2.0, system_cpu = 0.5, peak_mb = 150.0)

    summaries = collect_compute_summary(metrics_dir)

    assert len(summaries) == 1
    s = summaries[0]
    assert s.stage == "08b-annotation"
    assert s.tools == "mitos2"
    assert s.invocation_count == 2
    assert s.total_runtime_seconds == 8.0
    assert s.total_cpu_seconds == 7.5
    assert s.peak_memory_mb == 150.0

def test_collect_compute_summary_separates_stages_and_isolates_samples(tmp_path: Path) -> None:
    metrics_dir = tmp_path / "metrics"
    _write_metric(metrics_dir / "04d-assembly" / "sample1" / "getorganelle.sample1.1.metrics.json",
                  "getorganelle", runtime = 10.0, user_cpu = 9.0, system_cpu = 1.0, peak_mb = 500.0)
    _write_metric(metrics_dir / "04d-assembly" / "sample2" / "getorganelle.sample2.1.metrics.json",
                  "getorganelle", runtime = 20.0, user_cpu = 18.0, system_cpu = 2.0, peak_mb = 600.0)

    all_samples = collect_compute_summary(metrics_dir)
    assert len(all_samples) == 1
    assert all_samples[0].invocation_count == 2
    assert all_samples[0].total_runtime_seconds == 30.0

    sample1_only = collect_compute_summary(metrics_dir, sample_id = "sample1")
    assert len(sample1_only) == 1
    assert sample1_only[0].invocation_count == 1
    assert sample1_only[0].total_runtime_seconds == 10.0

def test_collect_compute_summary_lists_multiple_real_tools_for_a_shared_stage(tmp_path: Path) -> None:
    metrics_dir = tmp_path / "metrics"
    _write_metric(metrics_dir / "07-assembly_rescue" / "sample1" / "minimap2.sample1.1.metrics.json",
                  "minimap2", runtime = 2.0, user_cpu = 1.0, system_cpu = 0.0, peak_mb = 50.0)
    _write_metric(metrics_dir / "07-assembly_rescue" / "sample1" / "samtools.sample1.1.metrics.json",
                  "samtools", runtime = 1.0, user_cpu = 0.5, system_cpu = 0.0, peak_mb = 40.0)

    summaries = collect_compute_summary(metrics_dir)
    assert summaries[0].tools in ("minimap2, samtools", "samtools, minimap2")
    assert summaries[0].invocation_count == 2

def test_collect_compute_summary_empty_when_metrics_dir_missing(tmp_path: Path) -> None:
    assert collect_compute_summary(tmp_path / "does_not_exist") == ()

def test_collect_compute_summary_skips_a_malformed_record(tmp_path: Path) -> None:
    metrics_dir = tmp_path / "metrics"
    bad = metrics_dir / "08b-annotation" / "sample1" / "mitos2.sample1.1.metrics.json"
    bad.parent.mkdir(parents = True)
    bad.write_text("not real json")
    _write_metric(metrics_dir / "08b-annotation" / "sample1" / "mitos2.sample1.2.metrics.json",
                  "mitos2", runtime = 5.0, user_cpu = 4.0, system_cpu = 1.0, peak_mb = 100.0)

    summaries = collect_compute_summary(metrics_dir)
    assert summaries[0].invocation_count == 1

def test_collect_machine_info_returns_real_non_crashing_values() -> None:
    info = collect_machine_info()
    assert info.os_name
    assert info.architecture
    assert info.python_version
    assert info.cpu_logical_threads is None or info.cpu_logical_threads > 0
    assert info.cpu_physical_cores is None or info.cpu_physical_cores > 0
    # Best-effort real detection, every one of these is legitimately
    # allowed to be None (RAM speed/type needs root; GPU needs
    # nvidia-smi/lspci present; cpu model needs /proc/cpuinfo) --
    # the real assertion here is simply "never crashes."
    assert info.memory_speed_mhz is None or info.memory_speed_mhz > 0
    assert info.gpu_model is None or isinstance(info.gpu_model, str)

def test_cpu_model_and_clock_extracts_embedded_ghz_from_model_string() -> None:
    from mitopipeline.utils.reporting.compute_metrics import _cpu_model_and_clock
    cpuinfo = [{"model name": "Intel(R) Xeon(R) Platinum 8259CL CPU @ 2.50GHz", "cpu MHz": "2494.140"}]
    model, clock = _cpu_model_and_clock(cpuinfo)
    assert model == "Intel(R) Xeon(R) Platinum 8259CL CPU @ 2.50GHz"
    assert clock == 2500.0

def test_cpu_model_and_clock_falls_back_to_current_mhz_when_no_embedded_ghz() -> None:
    from mitopipeline.utils.reporting.compute_metrics import _cpu_model_and_clock
    cpuinfo = [{"model name": "Some Generic CPU", "cpu MHz": "1800.0"}]
    model, clock = _cpu_model_and_clock(cpuinfo)
    assert model == "Some Generic CPU"
    assert clock == 1800.0

def test_cpu_model_and_clock_none_for_empty_cpuinfo() -> None:
    from mitopipeline.utils.reporting.compute_metrics import _cpu_model_and_clock
    assert _cpu_model_and_clock([]) == (None, None)

def test_cpu_core_thread_counts_handles_hyperthreading() -> None:
    from mitopipeline.utils.reporting.compute_metrics import _cpu_core_thread_counts
    cpuinfo = [
        {"processor": "0", "physical id": "0", "core id": "0"},
        {"processor": "1", "physical id": "0", "core id": "0"},
        {"processor": "2", "physical id": "0", "core id": "1"},
        {"processor": "3", "physical id": "0", "core id": "1"},
    ]
    physical, logical = _cpu_core_thread_counts(cpuinfo)
    assert physical == 2
    assert logical == 4

def test_cpu_core_thread_counts_handles_multi_socket() -> None:
    from mitopipeline.utils.reporting.compute_metrics import _cpu_core_thread_counts
    cpuinfo = [
        {"processor": "0", "physical id": "0", "core id": "0"},
        {"processor": "1", "physical id": "1", "core id": "0"},
    ]
    physical, logical = _cpu_core_thread_counts(cpuinfo)
    assert physical == 2
    assert logical == 2

def test_cpu_core_thread_counts_none_for_empty_cpuinfo() -> None:
    from mitopipeline.utils.reporting.compute_metrics import _cpu_core_thread_counts
    assert _cpu_core_thread_counts([]) == (None, None)
