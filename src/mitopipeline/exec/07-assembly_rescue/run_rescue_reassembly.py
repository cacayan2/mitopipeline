"""
run_rescue_reassembly.py

Execution layer for assembly rescue's reassembly step (#406): takes
the rescued reads #426's mapping chain produced and re-runs
GetOrganelle through the SAME GetOrganelleRunner (BaseTool) the
primary attempt uses (#421), but with this attempt's own
getorganelle_profiles entry (config.yaml's documented retry
mechanism -- a different kmer/target_coverage/max_rounds each
attempt) and its own output directory, so repeated attempts never
collide and #402's contig selection can rank every attempt afterward.

Mirrors run_getorganelle_assembly.py's own run/classify/record
lifecycle deliberately -- this IS that same lifecycle, just fed
rescued reads and a named profile's tuning instead of the originally
trimmed reads and no tuning at all.
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
    parser = argparse.ArgumentParser(description = "Re-runs GetOrganelle on rescued reads with one getorganelle_profiles entry's tuning, and classifies the result.")

    parser.add_argument("--sample-id", type = str, required = True, help = "The sample id.")
    parser.add_argument("--r1", type = Path, required = True, help = "The rescued R1 reads (#426's output).")
    parser.add_argument("--r2", type = Path, required = True, help = "The rescued R2 reads (#426's output).")
    parser.add_argument("--organelle-type", type = str, default = "animal_mt", help = "The GetOrganelle organelle type.")
    parser.add_argument("--profile-name", type = str, required = True, help = "This attempt's getorganelle_profiles entry name (config.yaml), for logging/provenance only.")
    parser.add_argument("--kmer", type = str, default = None, help = "This profile's -k value (e.g. '21,33,55,85,115').")
    parser.add_argument("--target-coverage", type = float, default = None, help = "This profile's --reduce-reads-for-coverage value.")
    parser.add_argument("--max-rounds", type = int, default = None, help = "This profile's -R value.")
    parser.add_argument("--seed-file", type = Path, default = None, help = "A resolved seed FASTA for this attempt, if any (omit for unseeded).")
    parser.add_argument("--output-dir", type = Path, required = True, help = "This attempt's OWN output directory -- must differ from the primary attempt's and every other rescue attempt's.")
    parser.add_argument("--working-dir", type = Path, required = True, help = "The working directory.")
    parser.add_argument("--threads", type = int, required = True, help = "The number of threads to use.")
    parser.add_argument("--log-file", type = Path, required = True, help = "The path to the log file.")
    parser.add_argument("--global-log-file", type = Path, required = True, help = "The path to the global log file.")
    parser.add_argument("--metrics-dir", type = Path, default = None, help = "Where to persist this run's compute metrics as JSON. Omit to skip persistence.")

    return parser.parse_args()

def main() -> int:
    """
    Re-runs GetOrganelle on rescued reads with this attempt's profile
    tuning, and records its classified outcome.

    Returns:
        int: The return code for the process.
    """
    args = parse_args()

    logger = make_logger(
        name = "rescue_reassembly",
        log_file_path = args.log_file,
        global_log_file_path = args.global_log_file,
    )
    context = f"{{{args.sample_id}}}"

    sample = Sample(sample_id = args.sample_id, r1 = args.r1, r2 = args.r2)
    tool_options = {}
    if args.seed_file is not None:
        tool_options["seed_file"] = str(args.seed_file)
    if args.kmer is not None:
        tool_options["kmer"] = args.kmer
    if args.target_coverage is not None:
        tool_options["target_coverage"] = args.target_coverage
    if args.max_rounds is not None:
        tool_options["max_rounds"] = args.max_rounds

    runner = GetOrganelleRunner(
        working_dir = args.working_dir,
        output_dir = args.output_dir,
        sample = sample,
        organelle_type = args.organelle_type,
        tool_options = tool_options,
        threads = args.threads,
        logger = logger,
        metrics_dir = args.metrics_dir,
    )

    logger.info(f"{context} Initiating rescue reassembly (profile '{args.profile_name}': "
                f"kmer={args.kmer}, target_coverage={args.target_coverage}, max_rounds={args.max_rounds}).")

    try:
        result = runner.run()
    except SystemExit as e:
        logger.warning(f"{context} Rescue reassembly (profile '{args.profile_name}') did not produce a usable output (exit {e.code}); classifying as missing.")
        status = discover_assembly_artifacts(args.output_dir, logger)
        _write_status(args.output_dir, args.sample_id, status.status, success = False, profile_name = args.profile_name)
        return 1
    except Exception as e:
        logger.error(f"{context} Rescue reassembly (profile '{args.profile_name}') failed unexpectedly: {e}.")
        return 1

    if not result.success:
        logger.error(f"{context} Rescue reassembly exited with return code {result.return_code} after {result.runtime_seconds:.1f}s.")
        status = discover_assembly_artifacts(args.output_dir, logger)
        _write_status(args.output_dir, args.sample_id, status.status, success = False, profile_name = args.profile_name)
        return result.return_code or 1

    # Same reasoning as run_getorganelle_assembly.py: postprocess_outputs()
    # (inside run(), above) already standardized/cleaned up GetOrganelle's
    # native filenames, so runner.assembly_status (captured BEFORE that
    # cleanup) is the real classification -- re-deriving it now would
    # always find nothing.
    status = runner.assembly_status
    logger.info(f"{context} Rescue reassembly (profile '{args.profile_name}') complete after "
                f"{result.runtime_seconds:.1f}s; classified as '{status.status}' ({status.candidate_count} native candidate(s)).")
    _write_status(args.output_dir, args.sample_id, status.status, success = True, profile_name = args.profile_name)
    return 0

def _write_status(output_dir: Path, sample_id: str, status: str, success: bool, profile_name: str) -> None:
    """
    Records this rescue attempt's classification, same shape as the
    primary attempt's (plus which profile produced it), for #402's
    cross-attempt contig selection to read without re-deriving it.

    Args:
        output_dir (Path): This attempt's own output directory.
        sample_id (str): The sample id.
        status (str): The classification from discover_assembly_artifacts().
        success (bool): Whether the GetOrganelle command itself completed successfully.
        profile_name (str): Which getorganelle_profiles entry produced this attempt.

    Returns:
        None
    """
    output_dir.mkdir(parents = True, exist_ok = True)
    status_path = output_dir / f"{sample_id}.assembly_status.json"
    with status_path.open("w", encoding = "utf-8") as handle:
        json.dump({"sample_id": sample_id, "status": status, "success": success,
                  "profile_name": profile_name}, handle, indent = 4)

if __name__ == "__main__":
    raise SystemExit(main())
