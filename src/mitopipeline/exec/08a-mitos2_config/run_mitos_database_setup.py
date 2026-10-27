"""
run_mitos_database_setup.py

Execution layer for the job-level MITOS2 reference data setup:
downloads one refseqver's reference data from Zenodo (MITOS2's own
published source) before any sample's own annotation runs, so a
fresh install never requires a user to fetch and extract this
archive by hand.

A download failure is recorded, never raised past this script --
annotation with no usable reference data is a real failure for that
sample (MITOS2 cannot run at all), but it is still only that sample's
own failure, not a reason to abort the whole batch over a Zenodo
network hiccup.
"""

# Imports
import argparse
import json
from pathlib import Path
from mitopipeline.logging.logger_factory import make_logger
from mitopipeline.utils.annotation.mitos_database import MitosReferenceData

def parse_args() -> argparse.Namespace:
    """
    Parses the command line arguments for this script.

    Returns:
        argparse.Namespace: The parsed command line arguments.
    """
    parser = argparse.ArgumentParser(description = "Downloads and validates the local MITOS2 reference data.")

    parser.add_argument("--refdir", type = Path, required = True, help = "The directory MITOS2's own -R flag points at.")
    parser.add_argument("--refseqver", type = str, required = True, help = "The reference version to ensure is available (e.g. 'refseq89m').")
    parser.add_argument("--zenodo-record", type = str, required = True, help = "The Zenodo record id to download from.")
    parser.add_argument("--output-dir", type = Path, required = True, help = "The directory to write this step's own result JSON to.")
    parser.add_argument("--log-file", type = Path, required = True, help = "The path to the log file.")
    parser.add_argument("--global-log-file", type = Path, required = True, help = "The path to the global log file.")

    return parser.parse_args()

def _write_result(output_dir: Path, success: bool, version_dir: Path, reason: str) -> None:
    """Records this setup step's outcome for the per-sample annotation rule to read."""
    output_dir.mkdir(parents = True, exist_ok = True)
    path = output_dir / "mitos_database.json"
    with path.open("w", encoding = "utf-8") as handle:
        json.dump({
            "success": success, "refdir": str(version_dir), "reason": reason,
        }, handle, indent = 4)

def main() -> int:
    """
    Ensures the local MITOS2 reference data is downloaded and
    validated, recording the outcome.

    Returns:
        int: Always 0 -- a download failure is recorded, not raised; see module docstring.
    """
    args = parse_args()

    logger = make_logger(
        name = "mitos_database_setup",
        log_file_path = args.log_file,
        global_log_file_path = args.global_log_file,
    )

    reference_data = MitosReferenceData(
        refdir = args.refdir, refseqver = args.refseqver,
        record_id = args.zenodo_record, logger = logger,
    )

    try:
        resolved_dir = reference_data.prepare()
    except SystemExit as e:
        logger.warning(f"MITOS2 reference data setup failed (exit {e.code}). Annotation will fail for every sample this run.")
        _write_result(args.output_dir, success = False, version_dir = reference_data.version_dir,
                     reason = "download or validation failed")
        return 0
    except Exception as e:
        logger.warning(f"MITOS2 reference data setup failed: {e}. Annotation will fail for every sample this run.")
        _write_result(args.output_dir, success = False, version_dir = reference_data.version_dir,
                     reason = f"download or validation failed: {e}")
        return 0

    _write_result(args.output_dir, success = True, version_dir = resolved_dir, reason = "reference data ready")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
