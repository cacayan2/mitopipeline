"""
run_blast_hit_selection.py

Execution layer for BLAST hit ranking/selection (the hit-ranking module's own logic,
never wired to a CLI until now): ranks a sample's raw tabular BLAST
output, selects the top hits, and writes them as the top-hit TSV
run_blast_reference_retrieval.py and run_biological_validation.py
both read next.
"""

# Imports
import argparse
from pathlib import Path
from mitopipeline.logging.logger_factory import make_logger
from mitopipeline.utils.blast.blast_hits import parse_blast_tsv, select_top_blast_hits, write_top_blast_hits

def parse_args() -> argparse.Namespace:
    """
    Parses the command line arguments for this script.

    Returns:
        argparse.Namespace: The parsed command line arguments.
    """
    parser = argparse.ArgumentParser(description = "Ranks and selects the top BLAST hits for one sample.")

    parser.add_argument("--sample-id", type = str, required = True, help = "The sample id.")
    parser.add_argument("--blast-tsv", type = Path, required = True, help = "The raw tabular BLAST output written by run_blast.py.")
    parser.add_argument("--maximum-hits", type = int, default = 6, help = "The maximum number of top hits to select.")
    parser.add_argument("--output-dir", type = Path, required = True, help = "The directory to write the selected top-hit TSV to.")
    parser.add_argument("--log-file", type = Path, required = True, help = "The path to the log file.")
    parser.add_argument("--global-log-file", type = Path, required = True, help = "The path to the global log file.")

    return parser.parse_args()

def main() -> int:
    """
    Ranks and selects one sample's top BLAST hits.

    Returns:
        int: The return code for the process.
    """
    args = parse_args()

    logger = make_logger(
        name = "blast_hit_selection",
        log_file_path = args.log_file,
        global_log_file_path = args.global_log_file,
    )
    context = f"{{{args.sample_id}}}"

    if not args.blast_tsv.is_file():
        logger.error(f"{context} BLAST output not found: {args.blast_tsv}.")
        return 1

    hits = parse_blast_tsv(args.blast_tsv, logger = logger)
    # Zero hits is a valid outcome (BLAST ran successfully, found
    # nothing) -- an empty top-hit TSV is still written, not an error.
    selected = select_top_blast_hits(hits, maximum_hits = args.maximum_hits, logger = logger) if hits else []

    output_path = args.output_dir / f"{args.sample_id}.top_hits.tsv"
    write_top_blast_hits(selected, output_path, sample_id = args.sample_id, logger = logger)

    logger.info(f"{context} Selected {len(selected)} of {len(hits)} hit(s) to {output_path}.")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
