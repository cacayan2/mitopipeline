"""
run_tree_figure.py

Execution layer for rendering IQ-TREE's own real tree output
as a labeled PNG figure for the run-level report. Job-level, batch-
scoped -- the same scope as run_iqtree.py itself, not per-sample.

A missing or unsuccessful IQ-TREE result (phylogeny didn't run, or
found too few taxa to build a tree) is a normal, non-error outcome
this script must still complete for, not a failure of its own.
"""

# Imports
import argparse
import json
from pathlib import Path
from mitopipeline.logging.logger_factory import make_logger
from mitopipeline.utils.visualization.tree_figure import render_tree_figure

def parse_args() -> argparse.Namespace:
    """
    Parses the command line arguments for this script.

    Returns:
        argparse.Namespace: The parsed command line arguments.
    """
    parser = argparse.ArgumentParser(description = "Renders IQ-TREE's own real Newick tree as a PNG figure.")

    parser.add_argument("--iqtree-result-json", type = Path, required = True, help = "The result JSON written by run_iqtree.py.")
    parser.add_argument("--dpi", type = int, default = 600, help = "The rendered PNG's resolution. Defaults to 600.")
    parser.add_argument("--midpoint-root", action = "store_true", help = "Whether to midpoint-root the tree before drawing.")
    parser.add_argument("--output-dir", type = Path, required = True, help = "The directory to write the rendered PNG and this step's own result JSON to.")
    parser.add_argument("--log-file", type = Path, required = True, help = "The path to the log file.")
    parser.add_argument("--global-log-file", type = Path, required = True, help = "The path to the global log file.")

    return parser.parse_args()

def _write_result(output_dir: Path, success: bool, tree_figure_png: Path | None, reason: str) -> None:
    output_dir.mkdir(parents = True, exist_ok = True)
    with (output_dir / "tree_figure_result.json").open("w", encoding = "utf-8") as handle:
        json.dump({
            "success": success,
            "tree_figure_png": str(tree_figure_png) if tree_figure_png else None,
            "reason": reason,
        }, handle, indent = 4)

def main() -> int:
    """
    Renders IQ-TREE's own real tree as a PNG figure, when one exists.

    Returns:
        int: Always 0 -- a missing/unsuccessful upstream tree is a normal outcome, not a failure of this step.
    """
    args = parse_args()

    logger = make_logger(
        name = "tree_figure",
        log_file_path = args.log_file,
        global_log_file_path = args.global_log_file,
    )

    try:
        iqtree_result = json.loads(args.iqtree_result_json.read_text(encoding = "utf-8"))
    except (OSError, json.JSONDecodeError) as e:
        logger.info(f"No usable IQ-TREE result to render ({e}); skipping tree figure.")
        _write_result(args.output_dir, success = False, tree_figure_png = None, reason = f"could not read IQ-TREE result: {e}")
        return 0

    if not iqtree_result.get("success") or not iqtree_result.get("newick_path"):
        logger.info("IQ-TREE did not produce a usable tree this run; skipping tree figure.")
        _write_result(args.output_dir, success = False, tree_figure_png = None, reason = "no usable tree to render")
        return 0

    output_png = args.output_dir / "tree_figure.png"
    try:
        render_tree_figure(
            newick_path = iqtree_result["newick_path"], output_png = output_png,
            dpi = args.dpi, midpoint_root = args.midpoint_root, logger = logger,
        )
    except Exception as e:
        logger.warning(f"Tree figure rendering failed: {e}.")
        _write_result(args.output_dir, success = False, tree_figure_png = None, reason = f"rendering failed: {e}")
        return 0

    _write_result(args.output_dir, success = True, tree_figure_png = output_png, reason = "completed successfully")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
