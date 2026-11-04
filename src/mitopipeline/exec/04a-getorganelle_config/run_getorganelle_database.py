"""
run_getorganelle_database.py

Execution layer for the job-level GetOrganelle database check: makes
sure GetOrganelle's default LabelDatabase (BLAST) and SeedDatabase
(bowtie2) for the configured organelle type exist and are usable by
THIS environment's BLAST/bowtie2, downloading or rebuilding them as
needed, before any sample's assembly runs.

Unlike the optional MitoFish seed database, these are mandatory: with
an unusable LabelDatabase index GetOrganelle logs "slimming ... failed"
and every sample ends in "No valid assembly graph found!". So a failure
here exits non-zero and stops the run at once with the reason, instead
of letting a whole batch of assemblies fail one by one.
"""

# Imports
import argparse
import json
from pathlib import Path

from mitopipeline.logging.logger_factory import make_logger
from mitopipeline.utils.assembly.getorganelle_database import GetOrganelleDatabase

def parse_args() -> argparse.Namespace:
    """
    Parses the command line arguments for this script.

    Returns:
        argparse.Namespace: The parsed command line arguments.
    """
    parser = argparse.ArgumentParser(description = "Verifies (and repairs) GetOrganelle's default databases for this environment.")
    parser.add_argument("--organelle-type", type = str, required = True, help = "The GetOrganelle organelle type (e.g. animal_mt).")
    parser.add_argument("--config-dir", type = Path, default = None, help = "GetOrganelle's config/database directory. Defaults to $GETORG_PATH or ~/.GetOrganelle, as GetOrganelle itself does.")
    parser.add_argument("--output-dir", type = Path, required = True, help = "The directory to write this step's own result JSON to.")
    parser.add_argument("--log-file", type = Path, required = True, help = "The path to the log file.")
    parser.add_argument("--global-log-file", type = Path, required = True, help = "The path to the global log file.")
    return parser.parse_args()

def main() -> int:
    """
    Runs the database check and records the outcome.

    Returns:
        int: 0 when both databases are usable; 1 otherwise (the run must stop -- no assembly can succeed).
    """
    args = parse_args()
    logger = make_logger(
        name = "getorganelle_database",
        log_file_path = args.log_file,
        global_log_file_path = args.global_log_file,
    )

    result = GetOrganelleDatabase(organelle_type = args.organelle_type, logger = logger, config_dir = args.config_dir).ensure()

    args.output_dir.mkdir(parents = True, exist_ok = True)
    record_path = args.output_dir / "getorganelle_database.json"
    with record_path.open("w", encoding = "utf-8") as handle:
        json.dump({
            "success": result.success,
            "config_dir": str(result.config_dir),
            "organelle_type": result.organelle_type,
            "blast_version": result.blast_version,
            "downloaded": result.downloaded,
            "label_index_rebuilt": result.label_index_rebuilt,
            "seed_index_rebuilt": result.seed_index_rebuilt,
            "reason": result.reason,
            "actions": result.actions,
        }, handle, indent = 4)

    if not result.success:
        logger.error(f"GetOrganelle databases are not usable in this environment ({result.reason}); "
                     f"no assembly can succeed until this is fixed. See {record_path}.")
        return 1
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
