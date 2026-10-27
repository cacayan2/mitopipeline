"""
run_mitos2_annotation.py

Execution layer for MITOS2 annotation (#432): runs MITOS2 through
Mitos2Runner (BaseTool) against the resolved configuration #431 wrote,
parses the real gene features out of its own result.gff (#410), and
records a result JSON -- the feature count and gene list downstream
statistics/reporting stages read next.
"""

# Imports
import argparse
import json
from pathlib import Path
from mitopipeline.api.mitos2 import Mitos2Runner
from mitopipeline.logging.logger_factory import make_logger
from mitopipeline.utils.annotation.mitos_features import parse_mitos_gff

def parse_args() -> argparse.Namespace:
    """
    Parses the command line arguments for this script.

    Returns:
        argparse.Namespace: The parsed command line arguments.
    """
    parser = argparse.ArgumentParser(description = "Runs MITOS2 annotation for one sample and records the result.")

    parser.add_argument("--sample-id", type = str, required = True, help = "The sample id.")
    parser.add_argument("--config-json", type = Path, required = True, help = "The resolved MITOS2 configuration JSON written by run_mitos2_config.py.")
    parser.add_argument("--working-dir", type = Path, required = True, help = "The working directory.")
    parser.add_argument("--log-file", type = Path, required = True, help = "The path to the log file.")
    parser.add_argument("--global-log-file", type = Path, required = True, help = "The path to the global log file.")
    parser.add_argument("--metrics-dir", type = Path, default = None, help = "Where to persist this run's compute metrics as JSON. Omit to skip persistence.")

    return parser.parse_args()

def _write_result(output_dir: Path, sample_id: str, success: bool, feature_count: int | None,
                   gene_names: list[str] | None, gff_path: Path | None, reason: str) -> None:
    """
    Records this sample's annotation outcome for downstream statistics
    and reporting to read, without re-parsing the GFF or re-running
    MITOS2 just to know what it found.

    Args:
        output_dir (Path): Where to write the result record.
        sample_id (str): The sample id.
        success (bool): Whether MITOS2 itself completed successfully.
        feature_count (int | None): How many gene features were parsed, or None if the run failed outright.
        gene_names (list[str] | None): The parsed gene names, or None if the run failed outright.
        gff_path (Path | None): The real result.gff path, or None if the run failed outright.
        reason (str): A short, human-readable explanation.

    Returns:
        None
    """
    output_dir.mkdir(parents = True, exist_ok = True)
    path = output_dir / f"{sample_id}.mitos2_annotation_result.json"
    with path.open("w", encoding = "utf-8") as handle:
        json.dump({
            "sample_id": sample_id,
            "success": success,
            "feature_count": feature_count,
            "gene_names": gene_names,
            "gff_path": str(gff_path) if gff_path else None,
            "reason": reason,
        }, handle, indent = 4)

def main() -> int:
    """
    Runs MITOS2 annotation for one sample and records the result.

    Returns:
        int: The return code for the process.
    """
    args = parse_args()

    logger = make_logger(
        name = "mitos2_annotation",
        log_file_path = args.log_file,
        global_log_file_path = args.global_log_file,
    )
    context = f"{{{args.sample_id}}}"

    try:
        config = json.loads(args.config_json.read_text(encoding = "utf-8"))
    except (OSError, json.JSONDecodeError) as e:
        logger.error(f"{context} Could not read resolved MITOS2 config {args.config_json}: {e}")
        return 1

    output_dir = Path(config["annotation_output_dir"])

    runner = Mitos2Runner(
        working_dir = args.working_dir,
        output_dir = output_dir,
        sample_id = args.sample_id,
        input_fasta = Path(config["input_fasta"]),
        genetic_code = config["genetic_code"],
        refdir = Path(config["refdir"]),
        refseqver = config["refseqver"],
        circular = config["circular"],
        tool_options = config["tool_options"],
        logger = logger,
        metrics_dir = args.metrics_dir,
    )

    logger.info(f"{context} Initiating MITOS2 annotation (genetic_code={config['genetic_code']}, refseqver={config['refseqver']}).")

    try:
        result = runner.run()
    except SystemExit as e:
        logger.error(f"{context} MITOS2 did not produce a usable output (exit {e.code}).")
        _write_result(output_dir, args.sample_id, success = False, feature_count = None,
                      gene_names = None, gff_path = None, reason = f"mitos2 validation failed (exit {e.code})")
        return 1

    if not result.success:
        logger.error(f"{context} MITOS2 exited with return code {result.return_code} after {result.runtime_seconds:.1f}s.")
        _write_result(output_dir, args.sample_id, success = False, feature_count = None,
                      gene_names = None, gff_path = None, reason = f"mitos2 exited {result.return_code}")
        return result.return_code or 1

    try:
        features = parse_mitos_gff(runner.result_gff, logger = logger)
    except ValueError as e:
        logger.error(f"{context} MITOS2 ran successfully but its GFF output could not be parsed: {e}")
        _write_result(output_dir, args.sample_id, success = False, feature_count = None,
                      gene_names = None, gff_path = runner.result_gff, reason = f"gff parse error: {e}")
        return 1

    gene_names = [feature.gene_id for feature in features]
    logger.info(f"{context} MITOS2 complete after {result.runtime_seconds:.1f}s: {len(features)} gene feature(s) found.")
    _write_result(output_dir, args.sample_id, success = True, feature_count = len(features),
                 gene_names = gene_names, gff_path = runner.result_gff, reason = "completed successfully")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
