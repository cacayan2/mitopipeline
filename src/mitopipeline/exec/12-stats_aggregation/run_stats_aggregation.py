"""
run_stats_aggregation.py

Execution layer for per-sample statistics aggregation (#353-356):
extracts assembly/BLAST/annotation statistics from the real result
JSON each earlier stage already wrote, merges them by sample_id, and
records the result -- the reporting stage reads this next rather
than re-deriving statistics from raw tool output itself.
"""

# Imports
import argparse
import json
from pathlib import Path
from mitopipeline.logging.logger_factory import make_logger
from mitopipeline.utils.statistics.extract import (
    extract_annotation_statistics,
    extract_assembly_statistics,
    extract_blast_statistics,
    merge_sample_statistics,
)

def parse_args() -> argparse.Namespace:
    """
    Parses the command line arguments for this script.

    Returns:
        argparse.Namespace: The parsed command line arguments.
    """
    parser = argparse.ArgumentParser(description = "Aggregates one sample's assembly/BLAST/annotation statistics.")

    parser.add_argument("--sample-id", type = str, required = True, help = "The sample id.")
    parser.add_argument("--contig-selection-json", type = Path, required = True, help = "The selection JSON written by run_contig_selection.py.")
    parser.add_argument("--recruited-reads-qc-json", type = Path, default = None, help = "The selected attempt's recruited-reads QC record, if available.")
    parser.add_argument("--top-hits-tsv", type = Path, default = None, help = "The ranked top-hit TSV written by run_blast_hit_selection.py, if available.")
    parser.add_argument("--biological-validation-json", type = Path, default = None, help = "The validation result JSON written by run_biological_validation.py, if available.")
    parser.add_argument("--circular-coverage-json", type = Path, default = None, help = "The result JSON written by run_circular_coverage.py, for a circular selected attempt only, if available.")
    parser.add_argument("--annotation-result-json", type = Path, default = None, help = "The result JSON written by run_mitos2_annotation.py, if available.")
    parser.add_argument("--output-dir", type = Path, required = True, help = "The directory to write this step's own output to.")
    parser.add_argument("--log-file", type = Path, required = True, help = "The path to the log file.")
    parser.add_argument("--global-log-file", type = Path, required = True, help = "The path to the global log file.")

    return parser.parse_args()

def main() -> int:
    """
    Aggregates one sample's statistics and records the result.

    Returns:
        int: The return code for the process (always 0 -- a missing earlier stage is recorded as None fields, never a hard failure here).
    """
    args = parse_args()

    logger = make_logger(
        name = "stats_aggregation",
        log_file_path = args.log_file,
        global_log_file_path = args.global_log_file,
    )
    context = f"{{{args.sample_id}}}"

    assembly = extract_assembly_statistics(
        args.sample_id, args.contig_selection_json,
        recruited_reads_qc_json = args.recruited_reads_qc_json, logger = logger,
    )
    blast = None
    if args.top_hits_tsv is not None:
        blast = extract_blast_statistics(
            args.sample_id, args.top_hits_tsv,
            biological_validation_json = args.biological_validation_json,
            circular_coverage_json = args.circular_coverage_json, logger = logger,
        )
    annotation = None
    if args.annotation_result_json is not None:
        annotation = extract_annotation_statistics(args.sample_id, args.annotation_result_json, logger = logger)

    statistics = merge_sample_statistics(args.sample_id, assembly, blast, annotation)

    args.output_dir.mkdir(parents = True, exist_ok = True)
    output_path = args.output_dir / f"{args.sample_id}.statistics.json"
    with output_path.open("w", encoding = "utf-8") as handle:
        json.dump(statistics.to_dict(), handle, indent = 4)

    logger.info(f"{context} Statistics aggregated: "
               f"assembly={'yes' if assembly else 'no'}, blast={'yes' if blast else 'no'}, annotation={'yes' if annotation else 'no'}.")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
