"""
run_circular_map.py

Execution layer for circular mitochondrial genome map rendering
: parses the selected assembly attempt's own GFF/FASTA and
renders its circular genome map PNG. A rendering failure
for one sample -- a malformed GFF, an unreadable FASTA, an empty
annotation -- is recorded as this sample's own failed result, never
raised uncaught, so it never takes down any other sample's run.
"""

# Imports
import argparse
import json
from pathlib import Path
from mitopipeline.logging.logger_factory import make_logger
from mitopipeline.utils.visualization.circular_map import render_circular_genome_map

def parse_args() -> argparse.Namespace:
    """
    Parses the command line arguments for this script.

    Returns:
        argparse.Namespace: The parsed command line arguments.
    """
    parser = argparse.ArgumentParser(description = "Renders one sample's circular mitochondrial genome map.")

    parser.add_argument("--sample-id", type = str, required = True, help = "The sample id.")
    parser.add_argument("--gff-path", type = Path, required = True, help = "The selected attempt's MITOS2 result.gff.")
    parser.add_argument("--fasta-path", type = Path, required = True, help = "The selected attempt's genome FASTA.")
    parser.add_argument("--output-dir", type = Path, required = True, help = "The directory to write this step's own output to.")
    parser.add_argument("--dpi", type = int, default = 600, help = "The rendered PNG's resolution. Defaults to 600.")
    parser.add_argument("--log-file", type = Path, required = True, help = "The path to the log file.")
    parser.add_argument("--global-log-file", type = Path, required = True, help = "The path to the global log file.")

    return parser.parse_args()

def _write_result(output_dir: Path, sample_id: str, success: bool, output_png: Path | None,
                   feature_count: int | None, reason: str) -> None:
    """
    Records this sample's visualization outcome for downstream
    reporting to read, without re-deriving whether the render
    succeeded from the filesystem alone.

    Args:
        output_dir (Path): Where to write the result record.
        sample_id (str): The sample id.
        success (bool): Whether the render completed successfully.
        output_png (Path | None): The rendered PNG path, or None if the render failed.
        feature_count (int | None): How many drawable features were rendered, or None if the render failed.
        reason (str): A short, human-readable explanation.

    Returns:
        None
    """
    output_dir.mkdir(parents = True, exist_ok = True)
    path = output_dir / f"{sample_id}.circular_map_result.json"
    with path.open("w", encoding = "utf-8") as handle:
        json.dump({
            "sample_id": sample_id,
            "success": success,
            "output_png": str(output_png) if output_png else None,
            "feature_count": feature_count,
            "reason": reason,
        }, handle, indent = 4)

def main() -> int:
    """
    Renders one sample's circular mitochondrial genome map.

    Returns:
        int: The return code for the process.
    """
    args = parse_args()

    logger = make_logger(
        name = "circular_map",
        log_file_path = args.log_file,
        global_log_file_path = args.global_log_file,
    )
    context = f"{{{args.sample_id}}}"

    output_png = args.output_dir / f"{args.sample_id}.circular_map.png"

    try:
        map_data = render_circular_genome_map(
            sample_id = args.sample_id, gff_path = args.gff_path, fasta_path = args.fasta_path,
            output_png = output_png, dpi = args.dpi, logger = logger,
        )
    except (FileNotFoundError, ValueError) as e:
        logger.error(f"{context} Could not render circular genome map: {e}")
        _write_result(args.output_dir, args.sample_id, success = False, output_png = None,
                      feature_count = None, reason = str(e))
        return 1

    logger.info(f"{context} Circular genome map rendered: {len(map_data.features)} feature(s).")
    _write_result(args.output_dir, args.sample_id, success = True, output_png = output_png,
                 feature_count = len(map_data.features), reason = "completed successfully")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
