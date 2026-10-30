"""
run_mafft_alignment.py

Execution layer for MAFFT alignment: runs the combined,
batch-level FASTA through MafftRunner (BaseTool), then
validates the aligned output actually carries the same sequences
(by count and id) as the input before the tree stage ever sees it --
MAFFT reordering or losing a record would otherwise go unnoticed
until IQ-TREE failed on something harder to diagnose.
"""

# Imports
import argparse
import json
from pathlib import Path
from Bio import SeqIO
from mitopipeline.api.mafft import MafftRunner
from mitopipeline.logging.logger_factory import make_logger

def parse_args() -> argparse.Namespace:
    """
    Parses the command line arguments for this script.

    Returns:
        argparse.Namespace: The parsed command line arguments.
    """
    parser = argparse.ArgumentParser(description = "Aligns the combined phylogeny FASTA with MAFFT and validates the result.")

    parser.add_argument("--run-id", type = str, required = True, help = "A label for this batch-level run, for logging context.")
    parser.add_argument("--input-fasta", type = Path, required = True, help = "The combined, unaligned FASTA (the dataset-combination step's own output).")
    parser.add_argument("--output-fasta", type = Path, required = True, help = "The path to write the aligned FASTA to.")
    parser.add_argument("--strategy", type = str, default = "auto", help = "MAFFT's alignment strategy.")
    parser.add_argument("--threads", type = int, required = True, help = "The number of threads to use.")
    parser.add_argument("--working-dir", type = Path, required = True, help = "The working directory for the command (must already exist).")
    parser.add_argument("--output-dir", type = Path, required = True, help = "The directory to write this step's own result record to.")
    parser.add_argument("--log-file", type = Path, required = True, help = "The path to the log file.")
    parser.add_argument("--global-log-file", type = Path, required = True, help = "The path to the global log file.")
    parser.add_argument("--metrics-dir", type = Path, default = None, help = "Where to persist this run's compute metrics as JSON. Omit to skip persistence.")

    return parser.parse_args()

def _write_result(output_dir: Path, success: bool, sequence_count: int | None, reason: str) -> None:
    output_dir.mkdir(parents = True, exist_ok = True)
    with (output_dir / "mafft_alignment_result.json").open("w", encoding = "utf-8") as handle:
        json.dump({"success": success, "sequence_count": sequence_count, "reason": reason}, handle, indent = 4)

def main() -> int:
    """
    Aligns the combined phylogeny FASTA and records the result.

    Returns:
        int: The return code for the process.
    """
    args = parse_args()

    logger = make_logger(
        name = "mafft_alignment",
        log_file_path = args.log_file,
        global_log_file_path = args.global_log_file,
    )

    try:
        input_ids = {record.id for record in SeqIO.parse(args.input_fasta, "fasta")}
    except Exception as e:
        logger.error(f"Could not read input FASTA {args.input_fasta}: {e}")
        _write_result(args.output_dir, success = False, sequence_count = None, reason = f"could not read input FASTA: {e}")
        return 1

    runner = MafftRunner(
        working_dir = args.working_dir, input_fasta = args.input_fasta, output_fasta = args.output_fasta,
        sample_id = args.run_id, strategy = args.strategy, threads = args.threads, logger = logger,
        metrics_dir = args.metrics_dir,
    )

    try:
        result = runner.run()
    except SystemExit as e:
        logger.error(f"MAFFT did not produce a usable alignment (exit {e.code}).")
        _write_result(args.output_dir, success = False, sequence_count = None, reason = f"mafft validation failed (exit {e.code})")
        return 1

    if not result.success:
        logger.error(f"MAFFT exited with return code {result.return_code} after {result.runtime_seconds:.1f}s.")
        _write_result(args.output_dir, success = False, sequence_count = None, reason = f"mafft exited {result.return_code}")
        return result.return_code or 1

    aligned_records = list(SeqIO.parse(args.output_fasta, "fasta"))
    aligned_ids = {record.id for record in aligned_records}
    if aligned_ids != input_ids:
        logger.error(f"Aligned output doesn't match the input sequence set. Input: {sorted(input_ids)}, aligned: {sorted(aligned_ids)}.")
        _write_result(args.output_dir, success = False, sequence_count = len(aligned_records),
                     reason = "aligned sequence set doesn't match the input")
        return 1

    alignment_lengths = {len(record.seq) for record in aligned_records}
    if len(alignment_lengths) > 1:
        logger.error(f"Aligned sequences have inconsistent lengths: {alignment_lengths} -- not a valid alignment.")
        _write_result(args.output_dir, success = False, sequence_count = len(aligned_records),
                     reason = "aligned sequences have inconsistent lengths")
        return 1

    logger.info(f"MAFFT alignment complete after {result.runtime_seconds:.1f}s: {len(aligned_records)} sequence(s), "
               f"{alignment_lengths.pop()}bp aligned length.")
    _write_result(args.output_dir, success = True, sequence_count = len(aligned_records), reason = "completed successfully")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
