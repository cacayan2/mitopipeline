"""
run_sra_metadata.py

Execution layer for SRA metadata preparation: builds one
reviewable SRA submission row per raw sequencing sample from the
runtime manifest and the configured sequencing defaults. Job-level --
one table covers every sample in the run.
"""

# Imports
import argparse
import json
from pathlib import Path
from mitopipeline.logging.logger_factory import make_logger
from mitopipeline.utils.archival.sra_metadata import write_sra_metadata

def parse_args() -> argparse.Namespace:
    """
    Parses the command line arguments for this script.

    Returns:
        argparse.Namespace: The parsed command line arguments.
    """
    parser = argparse.ArgumentParser(description = "Builds the run's SRA metadata table.")

    parser.add_argument("--runtime-manifest", type = Path, required = True, help = "The runtime manifest TSV.")
    parser.add_argument("--defaults-json", type = Path, default = None, help = "The configured SRA defaults (stage_options.archival.sra), as a JSON file.")
    parser.add_argument("--output-dir", type = Path, required = True, help = "The directory to write the SRA metadata table to.")
    parser.add_argument("--log-file", type = Path, required = True, help = "The path to the log file.")
    parser.add_argument("--global-log-file", type = Path, required = True, help = "The path to the global log file.")

    return parser.parse_args()

def main() -> int:
    """
    Builds and writes the run's SRA metadata table.

    Returns:
        int: The return code for the process.
    """
    args = parse_args()

    logger = make_logger(
        name = "sra_metadata",
        log_file_path = args.log_file,
        global_log_file_path = args.global_log_file,
    )

    defaults = {}
    if args.defaults_json is not None and args.defaults_json.is_file():
        defaults = json.loads(args.defaults_json.read_text(encoding = "utf-8"))

    output_path = args.output_dir / "sra_metadata.tsv"
    try:
        write_sra_metadata(args.runtime_manifest, output_path, defaults = defaults, logger = logger)
    except (FileNotFoundError, ValueError) as e:
        logger.error(f"Could not build SRA metadata: {e}")
        return 1

    return 0

if __name__ == "__main__":
    raise SystemExit(main())
