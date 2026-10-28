"""
run_iqtree.py

Execution layer for IQ-TREE phylogenetic inference (#436): runs
MAFFT's aligned output through IqtreeRunner (BaseTool), then
validates the real Newick tree it produces -- it parses, and its tip
count matches the alignment's own sequence count -- before anything
downstream (reporting, visualization) treats it as a real tree.
"""

# Imports
import argparse
import json
from pathlib import Path
from Bio import Phylo, SeqIO
from mitopipeline.api.iqtree import IqtreeRunner
from mitopipeline.logging.logger_factory import make_logger

def parse_args() -> argparse.Namespace:
    """
    Parses the command line arguments for this script.

    Returns:
        argparse.Namespace: The parsed command line arguments.
    """
    parser = argparse.ArgumentParser(description = "Infers a maximum-likelihood phylogeny with IQ-TREE and validates the result.")

    parser.add_argument("--run-id", type = str, required = True, help = "A label for this batch-level run, for logging context.")
    parser.add_argument("--alignment", type = Path, required = True, help = "The aligned FASTA (MAFFT's own output, #348).")
    parser.add_argument("--output-prefix", type = Path, required = True, help = "The prefix IQ-TREE writes its own output under.")
    parser.add_argument("--output-newick", type = Path, required = True, help = "A stable, normalized copy of IQ-TREE's own real treefile.")
    parser.add_argument("--sequence-type", type = str, default = "DNA", help = "IQ-TREE's -st value.")
    parser.add_argument("--model", type = str, default = "MFP", help = "IQ-TREE's -m value.")
    parser.add_argument("--ultrafast-bootstrap", type = int, default = 1000, help = "Ultrafast bootstrap replicate count.")
    parser.add_argument("--sh-alrt", type = int, default = 1000, help = "SH-aLRT replicate count.")
    parser.add_argument("--random-seed", type = int, default = 12345, help = "A fixed random seed.")
    parser.add_argument("--threads", type = int, required = True, help = "The number of threads to use.")
    parser.add_argument("--working-dir", type = Path, required = True, help = "The working directory for the command (must already exist).")
    parser.add_argument("--output-dir", type = Path, required = True, help = "The directory to write this step's own result record to.")
    parser.add_argument("--log-file", type = Path, required = True, help = "The path to the log file.")
    parser.add_argument("--global-log-file", type = Path, required = True, help = "The path to the global log file.")
    parser.add_argument("--metrics-dir", type = Path, default = None, help = "Where to persist this run's compute metrics as JSON. Omit to skip persistence.")

    return parser.parse_args()

def _write_result(output_dir: Path, success: bool, tip_count: int | None, newick_path: Path | None, reason: str) -> None:
    output_dir.mkdir(parents = True, exist_ok = True)
    with (output_dir / "iqtree_result.json").open("w", encoding = "utf-8") as handle:
        json.dump({
            "success": success, "tip_count": tip_count,
            "newick_path": str(newick_path) if newick_path else None, "reason": reason,
        }, handle, indent = 4)

def main() -> int:
    """
    Infers a phylogeny with IQ-TREE and records the result.

    Returns:
        int: The return code for the process.
    """
    args = parse_args()

    logger = make_logger(
        name = "iqtree",
        log_file_path = args.log_file,
        global_log_file_path = args.global_log_file,
    )

    try:
        expected_tip_count = len(list(SeqIO.parse(args.alignment, "fasta")))
    except Exception as e:
        logger.error(f"Could not read alignment {args.alignment}: {e}")
        _write_result(args.output_dir, success = False, tip_count = None, newick_path = None, reason = f"could not read alignment: {e}")
        return 1

    # IQ-TREE's own real minimum is 3 sequences -- confirmed directly
    # against its installed error message. Checking this before
    # invoking the tool turns a cryptic "iqtree exited 2" (the exact
    # real exit code observed when this check is skipped) into a
    # readable, specific reason the report can show directly, with
    # the raw tool-level detail staying in the log where it belongs.
    MINIMUM_TAXA = 3
    if expected_tip_count < MINIMUM_TAXA:
        reason = f"Skipped: at least {MINIMUM_TAXA} assemblies are required to build a tree, but only {expected_tip_count} passed validation."
        logger.info(reason)
        _write_result(args.output_dir, success = False, tip_count = expected_tip_count, newick_path = None, reason = reason)
        return 0

    runner = IqtreeRunner(
        working_dir = args.working_dir, alignment_fasta = args.alignment, output_prefix = args.output_prefix,
        sample_id = args.run_id, sequence_type = args.sequence_type, model = args.model,
        ultrafast_bootstrap = args.ultrafast_bootstrap, sh_alrt = args.sh_alrt, random_seed = args.random_seed,
        threads = args.threads, logger = logger, metrics_dir = args.metrics_dir,
    )

    try:
        result = runner.run()
    except SystemExit as e:
        logger.error(f"IQ-TREE did not produce a usable tree (exit {e.code}).")
        _write_result(args.output_dir, success = False, tip_count = None, newick_path = None, reason = f"iqtree validation failed (exit {e.code})")
        return 1

    if not result.success:
        logger.error(f"IQ-TREE exited with return code {result.return_code} after {result.runtime_seconds:.1f}s.")
        _write_result(args.output_dir, success = False, tip_count = None, newick_path = None, reason = f"iqtree exited {result.return_code}")
        return result.return_code or 1

    try:
        tree = Phylo.read(runner.treefile, "newick")
    except Exception as e:
        logger.error(f"IQ-TREE's own treefile could not be parsed as valid Newick: {e}")
        _write_result(args.output_dir, success = False, tip_count = None, newick_path = None, reason = f"malformed newick tree: {e}")
        return 1

    tip_count = tree.count_terminals()
    if tip_count != expected_tip_count:
        logger.error(f"Tree has {tip_count} tips, expected {expected_tip_count} (one per aligned sequence).")
        _write_result(args.output_dir, success = False, tip_count = tip_count, newick_path = None,
                     reason = f"tip count {tip_count} does not match alignment sequence count {expected_tip_count}")
        return 1

    args.output_newick.parent.mkdir(parents = True, exist_ok = True)
    args.output_newick.write_text(runner.treefile.read_text())

    logger.info(f"IQ-TREE complete after {result.runtime_seconds:.1f}s: {tip_count} tips, "
               f"written to {args.output_newick}.")
    _write_result(args.output_dir, success = True, tip_count = tip_count, newick_path = args.output_newick, reason = "completed successfully")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
