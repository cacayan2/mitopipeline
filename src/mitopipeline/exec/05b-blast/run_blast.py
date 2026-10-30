"""
run_blast.py

Execution layer for BLAST identification: runs BLAST through
BlastRunner (BaseTool) against the resolved configuration the config stage wrote,
then records a result JSON (hit count, runtime, success) -- the
record the hit-selection and retrieval stages read next.

Zero hits is a valid, non-failure outcome (BlastRunner.validate_outputs
itself only rejects a missing output file, never an empty one) -- this
step's own result record distinguishes "ran successfully, found
nothing" from "the tool itself failed," so nothing downstream has to
guess which one happened from an empty TSV alone.
"""

# Imports
import argparse
import json
from pathlib import Path
from mitopipeline.api.blast import BlastRunner
from mitopipeline.logging.logger_factory import make_logger
from mitopipeline.models.sample.sample import Sample
from mitopipeline.utils.blast.blast_hits import parse_blast_tsv

def parse_args() -> argparse.Namespace:
    """
    Parses the command line arguments for this script.

    Returns:
        argparse.Namespace: The parsed command line arguments.
    """
    parser = argparse.ArgumentParser(description = "Runs BLAST identification for one sample and records the result.")

    parser.add_argument("--sample-id", type = str, required = True, help = "The sample id.")
    parser.add_argument("--config-json", type = Path, required = True, help = "The resolved BLAST configuration JSON written by run_blast_config.py.")
    parser.add_argument("--output-dir", type = Path, required = True, help = "The directory to write this step's own result record to.")
    parser.add_argument("--working-dir", type = Path, required = True, help = "The working directory.")
    parser.add_argument("--threads", type = int, required = True, help = "The number of threads to use.")
    parser.add_argument("--log-file", type = Path, required = True, help = "The path to the log file.")
    parser.add_argument("--global-log-file", type = Path, required = True, help = "The path to the global log file.")
    parser.add_argument("--metrics-dir", type = Path, default = None, help = "Where to persist this run's compute metrics as JSON. Omit to skip persistence.")

    return parser.parse_args()

def _write_result(output_dir: Path, sample_id: str, success: bool, hit_count: int | None,
                   output_file: Path | None, reason: str) -> None:
    """
    Records this sample's BLAST outcome for the hit-selection and retrieval stages to read, without
    re-deriving it (and without re-parsing the TSV just to know whether
    it ran).

    Args:
        output_dir (Path): Where to write the result record.
        sample_id (str): The sample id.
        success (bool): Whether BLAST itself completed successfully.
        hit_count (int | None): How many hits were found, or None if the run failed outright.
        output_file (Path | None): The tabular output path, or None if the run failed outright.
        reason (str): A short, human-readable explanation.

    Returns:
        None
    """
    output_dir.mkdir(parents = True, exist_ok = True)
    path = output_dir / f"{sample_id}.blast_result.json"
    with path.open("w", encoding = "utf-8") as handle:
        json.dump({
            "sample_id": sample_id,
            "success": success,
            "hit_count": hit_count,
            "output_file": str(output_file) if output_file else None,
            "reason": reason,
        }, handle, indent = 4)

def main() -> int:
    """
    Runs BLAST for one sample and records the result.

    Returns:
        int: The return code for the process.
    """
    args = parse_args()

    logger = make_logger(
        name = "blast",
        log_file_path = args.log_file,
        global_log_file_path = args.global_log_file,
    )
    context = f"{{{args.sample_id}}}"

    try:
        config = json.loads(args.config_json.read_text(encoding = "utf-8"))
    except (OSError, json.JSONDecodeError) as e:
        logger.error(f"{context} Could not read resolved BLAST config {args.config_json}: {e}")
        return 1

    output_file = Path(config["output_file"])
    query_fasta = Path(config["query_fasta"])

    # BlastRunner's own API takes a Sample object (shared BaseTool
    # convention for _log_context), even though BLAST has no paired
    # reads of its own to carry -- the query FASTA stands in for both,
    # since it's never actually read as r1/r2 by anything.
    sample = Sample(sample_id = args.sample_id, r1 = query_fasta, r2 = query_fasta)

    runner = BlastRunner(
        working_dir = args.working_dir,
        query_fasta = query_fasta,
        database = config["database"],
        output_file = output_file,
        sample = sample,
        mode = config["mode"],
        threads = args.threads,
        output_format = config["output_format"],
        logger = logger,
        metrics_dir = args.metrics_dir,
    )

    logger.info(f"{context} Initiating BLAST ({config['mode']} mode, database={config['database']}).")

    try:
        result = runner.run()
    except SystemExit as e:
        logger.error(f"{context} BLAST did not produce a usable output (exit {e.code}).")
        _write_result(args.output_dir, args.sample_id, success = False, hit_count = None,
                     output_file = None, reason = f"blast validation failed (exit {e.code})")
        return 1

    if not result.success:
        logger.error(f"{context} BLAST exited with return code {result.return_code} after {result.runtime_seconds:.1f}s.")
        _write_result(args.output_dir, args.sample_id, success = False, hit_count = None,
                     output_file = None, reason = f"blast exited {result.return_code}")
        return result.return_code or 1

    hits = parse_blast_tsv(output_file, logger = logger)
    logger.info(f"{context} BLAST complete after {result.runtime_seconds:.1f}s: {len(hits)} hit(s) found.")
    _write_result(args.output_dir, args.sample_id, success = True, hit_count = len(hits),
                 output_file = output_file, reason = "completed successfully")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
