"""
run_circular_coverage.py

Execution layer for circular-genome BLAST coverage correction
(the biological validation gate's own "blast_identity" check): for a circular assembly only,
re-runs BLAST against a doubled copy of the selected genome so an
origin-spanning match reads as one unbroken HSP, then records the
real coverage against the true, undoubled genome length --
run_biological_validation.py reads this result next, in place of the
raw top-hit qcovs, for circular samples.
"""

# Imports
import argparse
import json
from pathlib import Path
from mitopipeline.api.blast import BlastRunner
from mitopipeline.logging.logger_factory import make_logger
from mitopipeline.models.sample.sample import Sample
from mitopipeline.utils.blast.blast_hits import parse_blast_tsv
from mitopipeline.utils.blast.circular_coverage import compute_circular_query_coverage, double_fasta_sequence

def parse_args() -> argparse.Namespace:
    """
    Parses the command line arguments for this script.

    Returns:
        argparse.Namespace: The parsed command line arguments.
    """
    parser = argparse.ArgumentParser(description = "Corrects BLAST query coverage for a circular assembly.")

    parser.add_argument("--sample-id", type = str, required = True, help = "The sample id.")
    parser.add_argument("--contig-selection-json", type = Path, required = True, help = "The selection JSON written by run_contig_selection.py -- this step only ever runs for a circular selected attempt.")
    parser.add_argument("--database", type = str, required = True, help = "The BLAST database name (same one the BLAST config stage resolved).")
    parser.add_argument("--mode", type = str, default = "local", choices = ["local", "remote"], help = "BLAST mode. Defaults to local.")
    parser.add_argument("--working-dir", type = Path, required = True, help = "The working directory.")
    parser.add_argument("--threads", type = int, default = 4, help = "The number of threads to use.")
    parser.add_argument("--output-dir", type = Path, required = True, help = "The directory to write this step's own output to.")
    parser.add_argument("--log-file", type = Path, required = True, help = "The path to the log file.")
    parser.add_argument("--global-log-file", type = Path, required = True, help = "The path to the global log file.")

    return parser.parse_args()

def _write_result(output_dir: Path, sample_id: str, success: bool, corrected_coverage: float | None,
                   alignment_length: int | None, genome_length: int | None, reason: str) -> None:
    """Records this sample's circular-coverage correction for run_biological_validation.py to read next."""
    output_dir.mkdir(parents = True, exist_ok = True)
    path = output_dir / f"{sample_id}.circular_coverage_result.json"
    with path.open("w", encoding = "utf-8") as handle:
        json.dump({
            "sample_id": sample_id, "success": success,
            "corrected_query_coverage": corrected_coverage,
            "alignment_length": alignment_length, "genome_length": genome_length,
            "reason": reason,
        }, handle, indent = 4)

def main() -> int:
    """
    Corrects BLAST query coverage for one circular assembly.

    Returns:
        int: The return code for the process.
    """
    args = parse_args()

    logger = make_logger(
        name = "circular_coverage",
        log_file_path = args.log_file,
        global_log_file_path = args.global_log_file,
    )
    context = f"{{{args.sample_id}}}"

    try:
        selection = json.loads(args.contig_selection_json.read_text(encoding = "utf-8"))
    except (OSError, json.JSONDecodeError) as e:
        logger.error(f"{context} Could not read contig selection {args.contig_selection_json}: {e}")
        return 1

    if not selection.get("selected") or selection.get("status") != "circular":
        logger.info(f"{context} Selected attempt is not circular; circular-coverage correction does not apply.")
        _write_result(args.output_dir, args.sample_id, success = True, corrected_coverage = None,
                      alignment_length = None, genome_length = None, reason = "selected attempt is not circular")
        return 0

    args.output_dir.mkdir(parents = True, exist_ok = True)
    doubled_fasta = args.output_dir / f"{args.sample_id}.doubled.fasta"

    try:
        genome_length = double_fasta_sequence(selection["fasta_path"], doubled_fasta)
    except (FileNotFoundError, ValueError) as e:
        logger.error(f"{context} Could not build doubled query: {e}")
        _write_result(args.output_dir, args.sample_id, success = False, corrected_coverage = None,
                      alignment_length = None, genome_length = None, reason = str(e))
        return 1

    output_file = args.output_dir / f"{args.sample_id}.circular_blast.tsv"
    sample = Sample(sample_id = args.sample_id, r1 = doubled_fasta, r2 = doubled_fasta)
    runner = BlastRunner(
        working_dir = args.working_dir, query_fasta = doubled_fasta, database = args.database,
        output_file = output_file, sample = sample, mode = args.mode, threads = args.threads, logger = logger,
    )

    logger.info(f"{context} Initiating circular-coverage BLAST ({args.mode} mode, database={args.database}, genome_length={genome_length}bp).")

    try:
        result = runner.run()
    except SystemExit as e:
        logger.error(f"{context} Circular-coverage BLAST did not produce a usable output (exit {e.code}).")
        _write_result(args.output_dir, args.sample_id, success = False, corrected_coverage = None,
                      alignment_length = None, genome_length = genome_length, reason = f"blast validation failed (exit {e.code})")
        return 1

    if not result.success:
        logger.error(f"{context} Circular-coverage BLAST exited with return code {result.return_code}.")
        _write_result(args.output_dir, args.sample_id, success = False, corrected_coverage = None,
                      alignment_length = None, genome_length = genome_length, reason = f"blast exited {result.return_code}")
        return result.return_code or 1

    hits = parse_blast_tsv(output_file, logger = logger)
    corrected = compute_circular_query_coverage(hits, genome_length)
    if corrected is None:
        logger.warning(f"{context} Circular-coverage BLAST found no hits.")
        _write_result(args.output_dir, args.sample_id, success = True, corrected_coverage = None,
                      alignment_length = None, genome_length = genome_length, reason = "no hits found")
        return 0

    coverage, alignment_length = corrected
    logger.info(f"{context} Corrected circular query coverage: {coverage:.2f}% (alignment {alignment_length}bp / genome {genome_length}bp).")
    _write_result(args.output_dir, args.sample_id, success = True, corrected_coverage = coverage,
                 alignment_length = alignment_length, genome_length = genome_length, reason = "completed successfully")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
