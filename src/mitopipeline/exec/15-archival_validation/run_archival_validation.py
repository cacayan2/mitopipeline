"""
run_archival_validation.py

Execution layer for cross-record archival validation and packaging
(#374/#375/#450): combines every sample's own GenBank record into one
table, cross-checks sample identifiers across SRA/BioSample/GenBank,
builds the per-sample submission-readiness report, computes real MD5
checksums over every real archival file, and writes the submission
checklist -- a genuine "what's left to do" report, never a silent pass.
"""

# Imports
import argparse
import json
from pathlib import Path
from mitopipeline.logging.logger_factory import make_logger
from mitopipeline.utils.archival.genbank_metadata import write_genbank_metadata
from mitopipeline.utils.archival.models import GENBANK_COLUMNS
from mitopipeline.utils.archival.package import build_checksum_manifest, write_checksum_manifest, write_submission_checklist
from mitopipeline.utils.archival.validation import build_archival_validation, write_archival_validation

def parse_args() -> argparse.Namespace:
    """
    Parses the command line arguments for this script.

    Returns:
        argparse.Namespace: The parsed command line arguments.
    """
    parser = argparse.ArgumentParser(description = "Validates and packages the run's archival outputs.")

    parser.add_argument("--job-id", type = str, required = True, help = "The pipeline job id.")
    parser.add_argument("--runtime-manifest", type = Path, required = True, help = "The runtime manifest TSV.")
    parser.add_argument("--sra-metadata-tsv", type = Path, required = True, help = "The SRA metadata table written by run_sra_metadata.py.")
    parser.add_argument("--biosample-metadata-tsv", type = Path, required = True, help = "The BioSample metadata table written by run_biosample_metadata.py.")
    parser.add_argument("--genbank-records-dir", type = Path, required = True, help = "The directory containing every sample's own *.genbank_record.json.")
    parser.add_argument("--job-dir", type = Path, required = True, help = "The pipeline job's own output directory, used as the checksum manifest's base directory.")
    parser.add_argument("--output-dir", type = Path, required = True, help = "The directory to write the validation/packaging outputs to.")
    parser.add_argument("--log-file", type = Path, required = True, help = "The path to the log file.")
    parser.add_argument("--global-log-file", type = Path, required = True, help = "The path to the global log file.")

    return parser.parse_args()

def main() -> int:
    """
    Validates and packages the run's archival outputs.

    Returns:
        int: The return code for the process.
    """
    args = parse_args()

    logger = make_logger(
        name = "archival_validation",
        log_file_path = args.log_file,
        global_log_file_path = args.global_log_file,
    )

    genbank_records = []
    if args.genbank_records_dir.is_dir():
        for record_path in sorted(args.genbank_records_dir.glob("*.genbank_record.json")):
            record = json.loads(record_path.read_text(encoding = "utf-8"))
            genbank_records.append({column: record.get(column, "") for column in GENBANK_COLUMNS})

    args.output_dir.mkdir(parents = True, exist_ok = True)
    genbank_metadata_tsv = args.output_dir / "genbank_metadata.tsv"
    write_genbank_metadata(genbank_records, genbank_metadata_tsv, logger = logger)

    validation_records = build_archival_validation(
        args.runtime_manifest, args.sra_metadata_tsv, args.biosample_metadata_tsv,
        genbank_metadata_tsv, args.job_dir, logger = logger,
    )
    write_archival_validation(validation_records, args.output_dir / "archival_validation.tsv", logger = logger)

    real_files = [args.sra_metadata_tsv, args.biosample_metadata_tsv, genbank_metadata_tsv]
    for record in genbank_records:
        if record.get("assembly_fasta"): real_files.append(Path(record["assembly_fasta"]))
        if record.get("annotation_gff"): real_files.append(Path(record["annotation_gff"]))

    checksum_records = build_checksum_manifest(real_files, base_directory = args.job_dir)
    write_checksum_manifest(checksum_records, args.output_dir / "md5_checksums.tsv", logger = logger)
    write_submission_checklist(validation_records, args.job_id, args.output_dir / "submission_checklist.md", logger = logger)

    return 0

if __name__ == "__main__":
    raise SystemExit(main())
