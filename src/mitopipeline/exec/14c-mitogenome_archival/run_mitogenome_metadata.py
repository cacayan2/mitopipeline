"""
run_mitogenome_metadata.py

Execution layer for GenBank submission-folder preparation (#370/
#373): collects one sample's real selected assembly and MITOS2
annotation into its own accession-independent folder. Per-sample --
each folder depends only on that sample's own real artifacts.
"""

# Imports
import argparse
import json
from pathlib import Path
from mitopipeline.logging.logger_factory import make_logger
from mitopipeline.utils.archival.genbank_metadata import build_genbank_submission
from mitopipeline.utils.reporting.collect import read_sample_metadata

def parse_args() -> argparse.Namespace:
    """
    Parses the command line arguments for this script.

    Returns:
        argparse.Namespace: The parsed command line arguments.
    """
    parser = argparse.ArgumentParser(description = "Prepares one sample's GenBank submission folder.")

    parser.add_argument("--sample-id", type = str, required = True, help = "The sample id.")
    parser.add_argument("--runtime-manifest", type = Path, required = True, help = "The runtime manifest TSV.")
    parser.add_argument("--contig-selection-json", type = Path, required = True, help = "The selection JSON written by run_contig_selection.py.")
    parser.add_argument("--annotation-result-json", type = Path, required = True, help = "The result JSON written by run_mitos2_annotation.py.")
    parser.add_argument("--defaults-json", type = Path, default = None, help = "The configured GenBank defaults (stage_options.archival.biosample), as a JSON file.")
    parser.add_argument("--mode", type = str, default = "copy", choices = ["copy", "symlink"], help = "How to collect the real assembly/annotation files. Defaults to copy.")
    parser.add_argument("--output-dir", type = Path, required = True, help = "The directory to write this sample's submission folder and record under.")
    parser.add_argument("--log-file", type = Path, required = True, help = "The path to the log file.")
    parser.add_argument("--global-log-file", type = Path, required = True, help = "The path to the global log file.")

    return parser.parse_args()

def main() -> int:
    """
    Prepares one sample's GenBank submission folder and records the result.

    Returns:
        int: The return code for the process.
    """
    args = parse_args()

    logger = make_logger(
        name = "mitogenome_metadata",
        log_file_path = args.log_file,
        global_log_file_path = args.global_log_file,
    )

    metadata = read_sample_metadata(args.runtime_manifest, args.sample_id, logger = logger)
    organism = metadata.get("organism") or " ".join(
        part for part in (metadata.get("genus", ""), metadata.get("species", "")) if part
    )

    defaults = {}
    if args.defaults_json is not None and args.defaults_json.is_file():
        defaults = json.loads(args.defaults_json.read_text(encoding = "utf-8"))

    record = build_genbank_submission(
        args.sample_id, organism, args.contig_selection_json, args.annotation_result_json,
        args.output_dir, defaults = defaults, mode = args.mode, logger = logger,
    )

    record_path = args.output_dir / f"{args.sample_id}.genbank_record.json"
    record_path.write_text(json.dumps(record, indent = 4), encoding = "utf-8")

    return 0

if __name__ == "__main__":
    raise SystemExit(main())
