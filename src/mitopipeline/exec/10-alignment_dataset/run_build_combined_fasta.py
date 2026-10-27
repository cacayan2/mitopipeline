"""
run_build_combined_fasta.py

Execution layer for building the combined, batch-level phylogeny
FASTA (#344): reads every sample's own contig selection (#402) and
biological validation (#414) records, decides which are accepted,
and builds the combined dataset MAFFT/IQ-TREE align and compare next.

A sample's own biological validation FAIL verdict excludes it from
the combined dataset -- REVIEW does not, since that's a borderline
case for a person to look at, not something this step should itself
treat as rejected.
"""

# Imports
import argparse
import json
from pathlib import Path
from mitopipeline.logging.logger_factory import make_logger
from mitopipeline.utils.phylogeny.combined_fasta import CombinedFastaEntry, build_combined_fasta

def parse_args() -> argparse.Namespace:
    """
    Parses the command line arguments for this script.

    Returns:
        argparse.Namespace: The parsed command line arguments.
    """
    parser = argparse.ArgumentParser(description = "Builds the combined, batch-level phylogeny FASTA from every sample's own selection/validation.")

    parser.add_argument("--sample-id", type = str, required = True, action = "append", help = "A sample id considered for this batch. Repeat once per sample.")
    parser.add_argument("--contig-selection-dir", type = Path, required = True, help = "The root directory holding every sample's own contig_selection/biological_validation JSON (e.g. JOB_DIR/06-contig_selection).")
    parser.add_argument("--output-dir", type = Path, required = True, help = "The directory to write the combined FASTA to.")
    parser.add_argument("--log-file", type = Path, required = True, help = "The path to the log file.")
    parser.add_argument("--global-log-file", type = Path, required = True, help = "The path to the global log file.")

    return parser.parse_args()

def _load_entry(sample_id: str, contig_selection_dir: Path, logger) -> CombinedFastaEntry:
    """
    Reads one sample's own contig_selection.json and
    biological_validation.json into a CombinedFastaEntry. A missing
    or unreadable record is treated as "not accepted, no FASTA" --
    the same real outcome as a sample that genuinely had nothing
    selected -- rather than crashing the whole batch over one sample.

    Args:
        sample_id (str): The sample id.
        contig_selection_dir (Path): The root directory holding this sample's own records.
        logger: The logger to use.

    Returns:
        CombinedFastaEntry: The entry read for this sample.
    """
    sample_dir = contig_selection_dir / sample_id
    selection_path = sample_dir / f"{sample_id}.contig_selection.json"
    validation_path = sample_dir / f"{sample_id}.biological_validation.json"

    try:
        selection = json.loads(selection_path.read_text(encoding = "utf-8"))
    except (OSError, json.JSONDecodeError) as e:
        logger.warning(f"{{{sample_id}}} Could not read contig selection {selection_path} ({e!r}); treating as not accepted.")
        return CombinedFastaEntry(sample_id = sample_id, fasta_path = None, accepted = False)

    if not selection.get("selected"):
        return CombinedFastaEntry(sample_id = sample_id, fasta_path = None, accepted = False)

    fasta_path = Path(selection["fasta_path"])

    try:
        validation = json.loads(validation_path.read_text(encoding = "utf-8"))
        selected_attempt = validation.get("selected_attempt")
        verdicts = {r["attempt"]: r["overall_verdict"] for r in validation.get("results", [])}
        overall_verdict = verdicts.get(selected_attempt)
    except (OSError, json.JSONDecodeError, KeyError) as e:
        logger.warning(f"{{{sample_id}}} Could not read biological validation {validation_path} ({e!r}); treating as not accepted.")
        return CombinedFastaEntry(sample_id = sample_id, fasta_path = fasta_path, accepted = False)

    accepted = overall_verdict != "FAIL"
    return CombinedFastaEntry(sample_id = sample_id, fasta_path = fasta_path, accepted = accepted)

def main() -> int:
    """
    Builds and records the combined, batch-level phylogeny FASTA.

    Returns:
        int: The return code for the process.
    """
    args = parse_args()

    logger = make_logger(
        name = "build_combined_fasta",
        log_file_path = args.log_file,
        global_log_file_path = args.global_log_file,
    )

    entries = [_load_entry(sample_id, args.contig_selection_dir, logger) for sample_id in args.sample_id]

    output_path = args.output_dir / "combined.fasta"
    try:
        result = build_combined_fasta(entries, output_path, logger = logger)
    except ValueError as e:
        logger.error(f"Could not build combined FASTA: {e}")
        args.output_dir.mkdir(parents = True, exist_ok = True)
        with (args.output_dir / "combined_fasta_result.json").open("w", encoding = "utf-8") as handle:
            json.dump({"success": False, "reason": str(e), "included_sample_ids": [], "excluded_sample_ids": [e.sample_id for e in entries]}, handle, indent = 4)
        return 1

    with (args.output_dir / "combined_fasta_result.json").open("w", encoding = "utf-8") as handle:
        json.dump({
            "success": True,
            "output_path": str(result.output_path),
            "included_sample_ids": list(result.included_sample_ids),
            "excluded_sample_ids": list(result.excluded_sample_ids),
            "exclusion_reasons": result.exclusion_reasons,
        }, handle, indent = 4)

    return 0

if __name__ == "__main__":
    raise SystemExit(main())
