"""
compute_metrics.py

Aggregates the real per-invocation compute metrics every BaseTool
call already persists (runtime, CPU, peak memory -- see
BaseTool._write_metrics()) into a per-stage summary for the run-level
and per-sample reports, so the resource cost already being recorded
on every real run is actually visible somewhere (Emil's own direct
ask) instead of sitting unread under metrics/.
"""

# Imports
import json
import os
import platform
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path

@dataclass(frozen = True, slots = True)
class StageComputeSummary:
    """
    One pipeline stage's own real aggregated compute cost.

    Attributes:
        stage (str): The stage directory name (e.g. "08b-annotation").
        tools (str): The real tool name(s) invoked under this stage, comma-joined when more than one.
        invocation_count (int): How many real tool invocations this covers.
        total_runtime_seconds (float): Summed wall-clock runtime across every invocation.
        total_cpu_seconds (float): Summed user+system CPU time across every invocation.
        peak_memory_mb (float): The single highest peak memory observed across every invocation.
    """
    stage: str
    tools: str
    invocation_count: int
    total_runtime_seconds: float
    total_cpu_seconds: float
    peak_memory_mb: float

@dataclass(frozen = True, slots = True)
class MachineInfo:
    """
    The real execution host's own hardware/software environment --
    deliberately only non-identifying specs (no hostname, username,
    IP, or any other machine-identifying detail), for reproducibility
    context in a report, not for tracing a report back to a specific
    machine. Any individual spec that can't be real-detected (varies
    by host, virtualization, and installed tooling) is None rather
    than guessed.

    Attributes:
        os_name (str): The real OS name (e.g. "Linux").
        os_release (str): The real OS kernel/release string.
        architecture (str): The real CPU architecture (e.g. "x86_64").
        cpu_model (str | None): The real CPU model string (e.g. "Intel(R) Xeon(R) Platinum 8259CL CPU @ 2.50GHz"), or None if undetectable.
        cpu_physical_cores (int | None): The real physical core count, or None if undetectable.
        cpu_logical_threads (int | None): The real logical thread count (what os.cpu_count() reports), or None if undetectable.
        cpu_clock_mhz (float | None): The real CPU clock speed in MHz (from the model string, or the max/current frequency the kernel reports -- whichever real value was actually found), or None if undetectable.
        total_memory_gb (float | None): The real total system memory in GiB, or None if undetectable.
        memory_speed_mhz (int | None): The real RAM speed in MHz, or None if undetectable (needs dmidecode, usually root -- commonly unavailable, especially under virtualization).
        memory_type (str | None): The real RAM type (e.g. "DDR4"), or None if undetectable (same dmidecode dependency as memory_speed_mhz).
        gpu_model (str | None): The real GPU model string, or None if no GPU was detected. Never implies the pipeline uses it -- see the report's own fixed disclaimer text.
        python_version (str): The real Python interpreter version that ran the pipeline.
    """
    os_name: str
    os_release: str
    architecture: str
    cpu_model: str | None
    cpu_physical_cores: int | None
    cpu_logical_threads: int | None
    cpu_clock_mhz: float | None
    total_memory_gb: float | None
    memory_speed_mhz: int | None
    memory_type: str | None
    gpu_model: str | None
    python_version: str

def _total_memory_gb() -> float | None:
    """Reads real total system memory from /proc/meminfo (Linux only -- the pipeline's real execution environment); None elsewhere or if unreadable."""
    try:
        with open("/proc/meminfo", "r", encoding = "utf-8") as handle:
            for line in handle:
                if line.startswith("MemTotal:"):
                    kib = int(line.split()[1])
                    return round(kib / (1024 * 1024), 1)
    except (OSError, ValueError, IndexError):
        pass
    return None

def _read_proc_cpuinfo() -> list[dict[str, str]]:
    """Parses /proc/cpuinfo into one dict per real logical processor block; empty list outside Linux or if unreadable."""
    try:
        text = Path("/proc/cpuinfo").read_text(encoding = "utf-8")
    except OSError:
        return []
    blocks, current = [], {}
    for line in text.splitlines():
        if not line.strip():
            if current:
                blocks.append(current)
                current = {}
            continue
        if ":" in line:
            key, _, value = line.partition(":")
            current[key.strip()] = value.strip()
    if current:
        blocks.append(current)
    return blocks

def _cpu_model_and_clock(cpuinfo: list[dict[str, str]]) -> tuple[str | None, float | None]:
    """
    Reads the real CPU model string from /proc/cpuinfo's own "model
    name" field, and the clock speed it carries embedded (e.g. "@
    2.50GHz") when present -- the most reliable real source under
    virtualization, where /sys's own cpufreq files are often absent
    or report the virtualized host's own placeholder values. Falls
    back to /proc/cpuinfo's own real "cpu MHz" (current, not
    necessarily max) reading when no embedded value is found.
    """
    if not cpuinfo:
        return None, None
    model = cpuinfo[0].get("model name") or None
    clock_mhz = None
    if model:
        match = re.search(r"@\s*([\d.]+)\s*GHz", model, re.IGNORECASE)
        if match:
            clock_mhz = round(float(match.group(1)) * 1000, 1)
    if clock_mhz is None:
        raw = cpuinfo[0].get("cpu MHz")
        if raw:
            try:
                clock_mhz = round(float(raw), 1)
            except ValueError:
                pass
    return model, clock_mhz

def _cpu_core_thread_counts(cpuinfo: list[dict[str, str]]) -> tuple[int | None, int | None]:
    """Derives real physical core and logical thread counts from /proc/cpuinfo -- physical cores as the count of unique (physical id, core id) pairs (correct across multi-socket and hyperthreaded hosts alike), logical threads as the real processor-block count."""
    if not cpuinfo:
        return None, None
    logical_threads = len(cpuinfo)
    pairs = {(block.get("physical id", "0"), block.get("core id")) for block in cpuinfo if block.get("core id") is not None}
    physical_cores = len(pairs) if pairs else None
    return physical_cores, logical_threads

def _memory_speed_and_type() -> tuple[int | None, str | None]:
    """
    Reads real RAM speed/type via `dmidecode --type 17`, without
    sudo -- this genuinely requires root on most real hosts to read
    the DMI/SMBIOS tables, so this commonly (honestly) returns
    (None, None) rather than a guess; never escalates privileges
    itself. A 5s timeout guards against an unexpectedly slow/hanging
    call blocking report generation.
    """
    try:
        result = subprocess.run(
            ["dmidecode", "--type", "17"], capture_output = True, text = True, timeout = 5,
        )
    except (OSError, subprocess.TimeoutExpired):
        return None, None
    if result.returncode != 0:
        return None, None

    speed_mhz, memory_type = None, None
    for line in result.stdout.splitlines():
        line = line.strip()
        if line.startswith("Speed:") and speed_mhz is None:
            match = re.search(r"(\d+)\s*MT/s|\s*(\d+)\s*MHz", line)
            if match:
                speed_mhz = int(match.group(1) or match.group(2))
        elif line.startswith("Type:") and memory_type is None:
            value = line.split(":", 1)[1].strip()
            if value and value not in ("Unknown", "<OUT OF SPEC>"):
                memory_type = value
    return speed_mhz, memory_type

def _gpu_model() -> str | None:
    """
    Detects the real GPU model, trying `nvidia-smi` first (the only
    real path that works both on bare Linux and under WSL2's own GPU
    passthrough, when configured), falling back to `lspci` (bare
    Linux only -- WSL2's own virtualized PCI bus never exposes a real
    GPU this way). None when neither finds a real GPU, or neither
    tool is available -- never guessed.
    """
    try:
        result = subprocess.run(
            ["nvidia-smi", "--query-gpu=name", "--format=csv,noheader"],
            capture_output = True, text = True, timeout = 5,
        )
        if result.returncode == 0 and result.stdout.strip():
            return result.stdout.strip().splitlines()[0].strip()
    except (OSError, subprocess.TimeoutExpired):
        pass

    try:
        result = subprocess.run(["lspci"], capture_output = True, text = True, timeout = 5)
        if result.returncode == 0:
            for line in result.stdout.splitlines():
                if re.search(r"VGA compatible controller|3D controller|Display controller", line, re.IGNORECASE):
                    return line.split(":", 2)[-1].strip()
    except (OSError, subprocess.TimeoutExpired):
        pass

    return None

def collect_machine_info() -> MachineInfo:
    """
    Collects the real execution host's own non-identifying hardware/
    software specs for a report's reproducibility context. Every real
    external-command-dependent spec (RAM speed/type, GPU model) is
    best-effort and gracefully None when undetectable -- never
    fabricated or guessed.

    Returns:
        MachineInfo: The real collected specs.
    """
    cpuinfo = _read_proc_cpuinfo()
    cpu_model, cpu_clock_mhz = _cpu_model_and_clock(cpuinfo)
    cpu_physical_cores, cpu_logical_threads = _cpu_core_thread_counts(cpuinfo)
    memory_speed_mhz, memory_type = _memory_speed_and_type()

    return MachineInfo(
        os_name = platform.system(),
        os_release = platform.release(),
        architecture = platform.machine(),
        cpu_model = cpu_model,
        cpu_physical_cores = cpu_physical_cores,
        cpu_logical_threads = cpu_logical_threads if cpu_logical_threads is not None else os.cpu_count(),
        cpu_clock_mhz = cpu_clock_mhz,
        total_memory_gb = _total_memory_gb(),
        memory_speed_mhz = memory_speed_mhz,
        memory_type = memory_type,
        gpu_model = _gpu_model(),
        python_version = platform.python_version(),
    )

def _iter_metrics_records(metrics_dir: Path, sample_id: str | None = None) -> list[tuple[str, dict]]:
    """
    Reads every real *.metrics.json under metrics_dir (layout:
    metrics/{stage}/{sample}/{tool}.{sample}.{timestamp}.metrics.json),
    returning (stage, record) pairs. Filters to one sample's own
    records when sample_id is given. A single malformed/unreadable
    record is skipped, not fatal to the rest.
    """
    results: list[tuple[str, dict]] = []
    if not metrics_dir.is_dir():
        return results
    for path in metrics_dir.glob("*/*/*.metrics.json"):
        stage = path.parent.parent.name
        this_sample = path.parent.name
        if sample_id is not None and this_sample != sample_id:
            continue
        try:
            record = json.loads(path.read_text(encoding = "utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        results.append((stage, record))
    return results

def collect_compute_summary(metrics_dir: str | Path, sample_id: str | None = None) -> tuple[StageComputeSummary, ...]:
    """
    Aggregates every real compute-metrics record under metrics_dir
    into one StageComputeSummary per stage, sorted by stage name.

    Args:
        metrics_dir (str | Path): The job's own metrics/ directory.
        sample_id (str | None): Restrict to one sample's own records, or None for every sample (a run-level, batch-wide summary).

    Returns:
        tuple[StageComputeSummary, ...]: One summary per real stage with at least one recorded invocation. Empty when metrics_dir has nothing recorded.
    """
    by_stage: dict[str, list[dict]] = {}
    for stage, record in _iter_metrics_records(Path(metrics_dir), sample_id):
        by_stage.setdefault(stage, []).append(record)

    summaries = []
    for stage in sorted(by_stage):
        records = by_stage[stage]
        tools = ", ".join(dict.fromkeys(r.get("tool_name", "") for r in records if r.get("tool_name")))
        summaries.append(StageComputeSummary(
            stage = stage,
            tools = tools,
            invocation_count = len(records),
            total_runtime_seconds = sum(float(r.get("runtime_seconds") or 0.0) for r in records),
            total_cpu_seconds = sum(
                float(r.get("user_cpu_seconds") or 0.0) + float(r.get("system_cpu_seconds") or 0.0) for r in records
            ),
            peak_memory_mb = max((float(r.get("peak_memory_mb") or 0.0) for r in records), default = 0.0),
        ))
    return tuple(summaries)
