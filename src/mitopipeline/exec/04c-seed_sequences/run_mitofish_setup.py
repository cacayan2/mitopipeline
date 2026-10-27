"""
run_mitofish_setup.py

Execution layer for the job-level MitoFish database setup: ensures
the local MitoFish reference database is downloaded and validated
before any sample's own seed_selection runs, so a fresh install never
requires a user to fetch and place these files by hand. Resolves the
latest MitoFish release automatically unless a version is pinned.

A download failure is recorded, never raised past this script --
seed selection already treats a missing/unavailable database as a
normal fallback-to-unseeded outcome, so there is no correctness
reason to fail the whole run over a network hiccup fetching an
optional seed reference.
"""

# Imports
import argparse
import json
from pathlib import Path
from mitopipeline.logging.logger_factory import make_logger
from mitopipeline.utils.seed.mitofish_database import MitoFishDatabase, mitofish_urls_for_version, resolve_latest_mitofish_version

def parse_args() -> argparse.Namespace:
    """
    Parses the command line arguments for this script.

    Returns:
        argparse.Namespace: The parsed command line arguments.
    """
    parser = argparse.ArgumentParser(description = "Downloads and validates the local MitoFish reference database.")

    parser.add_argument("--database-root", type = Path, required = True, help = "The directory to store the (version-subdirectoried) MitoFish database under.")
    parser.add_argument("--version", type = str, default = None, help = "A specific MitoFish release to pin to, formatted 'YYYY.MM'. Omit to resolve the latest release automatically.")
    parser.add_argument("--output-dir", type = Path, required = True, help = "The directory to write this step's own result JSON to.")
    parser.add_argument("--log-file", type = Path, required = True, help = "The path to the log file.")
    parser.add_argument("--global-log-file", type = Path, required = True, help = "The path to the global log file.")

    return parser.parse_args()

def _write_result(output_dir: Path, success: bool, database_dir: Path, version: str | None, reason: str) -> None:
    """Records this setup step's outcome for seed_selection to read next."""
    output_dir.mkdir(parents = True, exist_ok = True)
    path = output_dir / "mitofish_database.json"
    with path.open("w", encoding = "utf-8") as handle:
        json.dump({
            "success": success, "database_dir": str(database_dir),
            "version": version, "reason": reason,
        }, handle, indent = 4)

def main() -> int:
    """
    Ensures the local MitoFish reference database is downloaded and
    validated, recording the outcome for seed_selection to read.

    Returns:
        int: Always 0 -- a download failure is recorded, not raised; see module docstring.
    """
    args = parse_args()

    logger = make_logger(
        name = "mitofish_setup",
        log_file_path = args.log_file,
        global_log_file_path = args.global_log_file,
    )

    try:
        version = args.version or resolve_latest_mitofish_version(logger)
        tables_url, sequences_url = mitofish_urls_for_version(version)
    except Exception as e:
        logger.warning(f"Could not resolve a MitoFish release to download: {e}. Seed selection will fall back to unseeded assembly.")
        _write_result(args.output_dir, success = False, database_dir = args.database_root, version = args.version,
                     reason = f"version resolution failed: {e}")
        return 0

    database = MitoFishDatabase(
        database_dir = args.database_root, version = version,
        tables_url = tables_url, sequences_url = sequences_url, logger = logger,
    )

    try:
        resolved_dir = database.prepare()
    except SystemExit:
        logger.warning("MitoFish database preparation failed; seed selection will fall back to unseeded assembly.")
        _write_result(args.output_dir, success = False, database_dir = database.database_dir, version = version,
                     reason = "download or validation failed")
        return 0

    _write_result(args.output_dir, success = True, database_dir = resolved_dir, version = version,
                 reason = "database ready")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
