"""
run_getorganelle_assembly.py

Execution layer for GetOrganelle assembly: runs the tool through
GetOrganelleRunner (BaseTool), then classifies the result with
discover_assembly_artifacts() and records a status JSON alongside the
standardized assembly output -- the runtime/failure record the next
stage (and eventual rescue decision) reads.

Reads the resolved configuration the config stage wrote (organelle type, seed/
no-seed, output dir) rather than re-deriving it, so the two stages
never disagree about what this sample's assembly run was configured to do.
"""

# Imports
import argparse
import json
from pathlib import Path
from mitopipeline.api.getorganelle import GetOrganelleRunner, discover_assembly_artifacts
from mitopipeline.logging.logger_factory import make_logger
from mitopipeline.models.sample.sample import Sample

def parse_args() -> argparse.Namespace:
    """
    Parses the command line arguments for this script.

    Returns:
        argparse.Namespace: The parsed command line arguments.
    """
    parser = argparse.ArgumentParser(description = "Runs GetOrganelle assembly for one sample and classifies the result.")

    parser.add_argument("--sample-id", type = str, required = True, help = "The sample id.")
    parser.add_argument("--r1", type = str, required = True, help = "The path to the trimmed R1 fastq file.")
    parser.add_argument("--r2", type = str, required = True, help = "The path to the trimmed R2 fastq file.")
    parser.add_argument("--config-json", type = Path, required = True, help = "The resolved GetOrganelle configuration JSON written by run_getorganelle_config.py.")
    parser.add_argument("--working-dir", type = Path, required = True, help = "The path to the working directory.")
    parser.add_argument("--log-file", type = Path, required = True, help = "The path to the log file.")
    parser.add_argument("--global-log-file", type = Path, required = True, help = "The path to the global log file.")
    parser.add_argument("--threads", type = int, required = True, help = "The number of threads to use.")
    parser.add_argument("--metrics-dir", type = Path, default = None, help = "Where to persist this run's compute metrics (CPU/memory/runtime) as JSON. Omit to skip persistence.")

    return parser.parse_args()

def main() -> int:
    """
    Runs GetOrganelle assembly for one sample and records its
    classified outcome.

    Returns:
        int: The return code for the process.
    """
    args = parse_args()

    logger = make_logger(
        name = "getorganelle_assembly",
        log_file_path = args.log_file,
        global_log_file_path = args.global_log_file,
    )
    context = f"{{{args.sample_id}}}"

    try:
        config = json.loads(args.config_json.read_text(encoding = "utf-8"))
    except (OSError, json.JSONDecodeError) as e:
        logger.error(f"{context} Could not read resolved GetOrganelle config {args.config_json}: {e}")
        return 1

    sample = Sample(sample_id = args.sample_id, r1 = Path(args.r1), r2 = Path(args.r2))
    tool_options = {"seed_file": config["seed_file"]} if config.get("seeded") else {}
    output_dir = Path(config["assembly_output_dir"])

    # A re-invocation (Snakemake's --rerun-incomplete re-running a job
    # it can no longer trust, a manual retry, etc.) can land here after
    # GetOrganelle already finished this sample successfully. By then
    # postprocess_outputs() has already standardized the directory to
    # data.fasta/graph.gfa and deleted GetOrganelle's native output, so
    # the directory is neither fresh (GetOrganelle refuses to reuse an
    # existing -o directory) nor does it still hold the native filenames
    # discover_assembly_artifacts() looks for -- a genuine prior success
    # would otherwise be reclassified as "missing". Trusting the
    # standardized output plus its own recorded status avoids re-running
    # GetOrganelle at all in that case.
    prior_status = _prior_success(output_dir, args.sample_id)
    if prior_status is not None:
        logger.info(f"{context} GetOrganelle assembly already completed for this sample (classified as '{prior_status}'); skipping re-run.")
        return 0

    runner = GetOrganelleRunner(
        working_dir = Path(args.working_dir),
        output_dir = output_dir,
        sample = sample,
        organelle_type = config["organelle_type"],
        tool_options = tool_options,
        threads = args.threads,
        logger = logger,
        metrics_dir = args.metrics_dir,
    )

    logger.info(f"{context} Initiating GetOrganelle assembly ({'seeded' if config.get('seeded') else 'unseeded'}, organelle_type={config['organelle_type']}).")

    # Recording runtime and failures: a validate_outputs() failure (no
    # native candidate at all) exits via sys.exit inside run() -- caught
    # here as SystemExit so a failed assembly is reported and classified
    # as "missing", not left as an uncaught crash. Any other exception
    # is a genuine, unexpected failure, not a normal "no assembly" outcome.
    try:
        result = runner.run()
    except SystemExit as e:
        logger.warning(f"{context} GetOrganelle assembly did not produce a usable output (exit {e.code}); classifying as missing.")
        status = discover_assembly_artifacts(output_dir, logger)
        _write_status(output_dir, args.sample_id, status.status, success = False)
        return 1
    except Exception as e:
        logger.error(f"{context} GetOrganelle assembly failed unexpectedly: {e}.")
        return 1

    if not result.success:
        logger.error(f"{context} GetOrganelle assembly exited with return code {result.return_code} after {result.runtime_seconds:.1f}s.")
        status = discover_assembly_artifacts(output_dir, logger)
        _write_status(output_dir, args.sample_id, status.status, success = False)
        return result.return_code or 1

    # Not re-deriving the classification here: postprocess_outputs()
    # (run inside runner.run(), above) already standardized/cleaned up
    # GetOrganelle's native output filenames, so discover_assembly_artifacts()
    # would find nothing if called again now. runner.assembly_status is
    # the classification captured before that cleanup happened.
    status = runner.assembly_status
    logger.info(f"{context} GetOrganelle assembly complete after {result.runtime_seconds:.1f}s; classified as '{status.status}' ({status.candidate_count} native candidate(s)).")
    _write_status(output_dir, args.sample_id, status.status, success = True)
    return 0

def _prior_success(output_dir: Path, sample_id: str) -> str | None:
    """
    Checks whether this sample's output_dir already holds a validated,
    successful assembly from an earlier invocation -- the standardized
    data.fasta/graph.gfa pair plus a status record reporting success.
    All three must be present and non-empty; a partial or failed prior
    attempt (missing file, empty file, or a recorded failure) is never
    treated as a prior success, so a genuinely incomplete directory
    still gets a real re-run.

    Args:
        output_dir (Path): The sample's assembly output directory.
        sample_id (str): The sample id.

    Returns:
        str | None: The prior classification ("circular" or "fragmented") if a validated success is already present, otherwise None.
    """
    status_path = output_dir / f"{sample_id}.assembly_status.json"
    standardized_fasta = output_dir / "data.fasta"
    standardized_gfa = output_dir / "graph.gfa"

    for path in (status_path, standardized_fasta, standardized_gfa):
        if not path.is_file() or path.stat().st_size == 0:
            return None

    try:
        record = json.loads(status_path.read_text(encoding = "utf-8"))
    except (OSError, json.JSONDecodeError):
        return None

    if not record.get("success"):
        return None
    return record.get("status")

def _write_status(output_dir: Path, sample_id: str, status: str, success: bool) -> None:
    """
    Records this sample's assembly classification alongside its
    standardized output, for the rescue trigger and any later
    stage to read without re-deriving it.

    Args:
        output_dir (Path): The sample's assembly output directory.
        sample_id (str): The sample id.
        status (str): The classification from discover_assembly_artifacts() ("circular", "fragmented", or "missing").
        success (bool): Whether the GetOrganelle command itself completed successfully.

    Returns:
        None
    """
    output_dir.mkdir(parents = True, exist_ok = True)
    status_path = output_dir / f"{sample_id}.assembly_status.json"
    with status_path.open("w", encoding = "utf-8") as handle:
        json.dump({"sample_id": sample_id, "status": status, "success": success}, handle, indent = 4)

if __name__ == "__main__":
    raise SystemExit(main())
