"""
run_table2asn_validation.py

Execution layer for local table2asn validation: runs NCBI's own real
validation tool against one sample's real GenBank submission folder
(the FASTA + 5-column feature table 14c already prepared), recording
a typed PASS/FAIL result. Entirely local -- never contacts or
submits anything to NCBI; see the module docstring in
table2asn_validation.py.
"""

# Imports
import argparse
import json
import shutil
from pathlib import Path
from mitopipeline.api.table2asn import Table2AsnRunner
from mitopipeline.logging.logger_factory import make_logger
from mitopipeline.utils.archival.table2asn_binary import Table2AsnBinary
from mitopipeline.utils.archival.table2asn_validation import build_validation_result, tag_fasta_with_ncbi_defline

def parse_args() -> argparse.Namespace:
    """
    Parses the command line arguments for this script.

    Returns:
        argparse.Namespace: The parsed command line arguments.
    """
    parser = argparse.ArgumentParser(description = "Validates one sample's GenBank submission folder with table2asn.")

    parser.add_argument("--sample-id", type = str, required = True, help = "The sample id.")
    parser.add_argument("--genbank-record-json", type = Path, required = True, help = "The record JSON written by run_mitogenome_metadata.py.")
    parser.add_argument("--binary-path", type = Path, required = True, help = "Where the real table2asn binary should live once prepared (auto-installed if missing).")
    parser.add_argument("--submission-template", type = Path, required = True, help = "The real .sbt submission template.")
    parser.add_argument("--topology", type = str, default = "circular", help = "Defaults to 'circular'.")
    parser.add_argument("--genome-location", type = str, default = "mitochondrion", help = "Defaults to 'mitochondrion'.")
    parser.add_argument("--working-dir", type = Path, required = True, help = "The working directory for the command.")
    parser.add_argument("--output-dir", type = Path, required = True, help = "Where to write this sample's real table2asn inputs/outputs and result JSON.")
    parser.add_argument("--log-file", type = Path, required = True, help = "The path to the log file.")
    parser.add_argument("--global-log-file", type = Path, required = True, help = "The path to the global log file.")
    parser.add_argument("--metrics-dir", type = Path, default = None, help = "Where to persist this run's compute metrics as JSON. Omit to skip persistence.")

    return parser.parse_args()

def _write_skipped_result(output_dir: Path, sample_id: str, reason: str) -> None:
    """Records a real, honest skip -- never a fabricated PASS -- when the prerequisite files genuinely aren't there."""
    output_dir.mkdir(parents = True, exist_ok = True)
    path = output_dir / f"{sample_id}.table2asn_validation_result.json"
    path.write_text(json.dumps({
        "sample_id": sample_id, "status": "SKIPPED", "return_code": None, "sqn_file": "", "validation_file": "",
        "stats_file": "", "discrepancy_file": "", "genbank_flatfile": "", "reject_count": 0, "fatal_count": 0,
        "error_count": 0, "warning_count": 0, "info_count": 0, "blockers": reason,
    }, indent = 4), encoding = "utf-8")

def main() -> int:
    """
    Validates one sample's real GenBank submission folder with
    table2asn and records the result.

    Returns:
        int: Always 0 -- a validation FAIL is a real, recorded outcome, never a reason to abort the batch.
    """
    args = parse_args()

    logger = make_logger(name = "table2asn_validation", log_file_path = args.log_file, global_log_file_path = args.global_log_file)
    context = f"{{{args.sample_id}}}"

    try:
        record = json.loads(args.genbank_record_json.read_text(encoding = "utf-8"))
    except (OSError, json.JSONDecodeError) as e:
        logger.error(f"{context} Could not read GenBank record {args.genbank_record_json}: {e}")
        _write_skipped_result(args.output_dir, args.sample_id, f"could not read GenBank record: {e}")
        return 0

    assembly_fasta, annotation_tbl, organism = record.get("assembly_fasta", ""), record.get("annotation_tbl", ""), record.get("organism", "")
    if not assembly_fasta or not Path(assembly_fasta).is_file() or not annotation_tbl or not Path(annotation_tbl).is_file():
        logger.info(f"{context} No real assembly FASTA and/or feature table available; skipping table2asn validation.")
        _write_skipped_result(args.output_dir, args.sample_id, "missing assembly_fasta and/or annotation_tbl")
        return 0

    try:
        binary = Table2AsnBinary(binary_path = args.binary_path, logger = logger).prepare()
    except Exception as e:
        logger.error(f"{context} Could not prepare the table2asn binary: {e}")
        _write_skipped_result(args.output_dir, args.sample_id, f"table2asn binary unavailable: {e}")
        return 0

    sample_root = args.output_dir / args.sample_id
    input_dir = sample_root / "input"
    output_subdir = sample_root / "output"
    input_dir.mkdir(parents = True, exist_ok = True)

    fasta_input = input_dir / f"{args.sample_id}.fsa"
    feature_table_input = input_dir / f"{args.sample_id}.tbl"
    tag_fasta_with_ncbi_defline(Path(assembly_fasta), fasta_input, args.sample_id, organism, args.topology, args.genome_location)
    shutil.copy2(annotation_tbl, feature_table_input)

    runner = Table2AsnRunner(
        working_dir = sample_root, fasta_path = fasta_input, feature_table_path = feature_table_input,
        submission_template_path = args.submission_template, output_basename = output_subdir / args.sample_id,
        logger = logger, binary_path = str(binary), metrics_dir = args.metrics_dir,
    )

    logger.info(f"{context} Running table2asn validation.")
    try:
        result = runner.run()
    except SystemExit as e:
        logger.error(f"{context} table2asn ran but produced no usable output (exit {e.code}).")
        _write_skipped_result(args.output_dir, args.sample_id, f"table2asn produced no usable output (exit {e.code})")
        return 0

    stdout_path, stderr_path = output_subdir / "table2asn.stdout.txt", output_subdir / "table2asn.stderr.txt"
    stdout_path.write_text(result.stdout, encoding = "utf-8")
    stderr_path.write_text(result.stderr, encoding = "utf-8")

    validation = build_validation_result(args.sample_id, result.return_code, output_subdir / args.sample_id, stdout_path, stderr_path)
    logger.info(f"{context} table2asn validation: {validation.status} after {result.runtime_seconds:.1f}s"
               f"{' (' + validation.blockers + ')' if validation.blockers else ''}.")

    args.output_dir.mkdir(parents = True, exist_ok = True)
    (args.output_dir / f"{args.sample_id}.table2asn_validation_result.json").write_text(
        json.dumps(validation.to_dict(), indent = 4), encoding = "utf-8")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
