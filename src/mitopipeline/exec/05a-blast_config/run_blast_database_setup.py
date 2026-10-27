"""
run_blast_database_setup.py

Execution layer for the job-level local BLAST reference database
setup: downloads the configured taxonomic scope's complete RefSeq
mitochondrial genomes from NCBI Entrez and builds the local BLAST
database before any sample's own BLAST stages run, so a fresh install
never requires a user to build this database by hand.

A download/build failure is recorded, never raised past this script.
The biological validation gate's own "blast_identity" check already
treats missing BLAST evidence as a normal REVIEW outcome, not a
crash, so there is no correctness reason to fail the whole run over
an NCBI network hiccup.
"""

# Imports
import argparse
import json
from pathlib import Path
from mitopipeline.logging.logger_factory import make_logger
from mitopipeline.utils.blast.reference_database import setup_reference_database

def parse_args() -> argparse.Namespace:
    """
    Parses the command line arguments for this script.

    Returns:
        argparse.Namespace: The parsed command line arguments.
    """
    parser = argparse.ArgumentParser(description = "Downloads and builds the local BLAST reference database.")

    parser.add_argument("--reference-fasta", type = Path, required = True, help = "The path to write/read the reference FASTA.")
    parser.add_argument("--database-prefix", type = Path, required = True, help = "The path prefix for the built BLAST database.")
    parser.add_argument("--database-title", type = str, required = True, help = "A human-readable title recorded in the database.")
    parser.add_argument("--taxon-scope", type = str, required = True, help = "The NCBI organism-level scope to search within (e.g. 'Actinopterygii').")
    parser.add_argument("--entrez-email", type = str, required = True, help = "The contact email NCBI requires for Entrez requests.")
    parser.add_argument("--entrez-api-key", type = str, default = None, help = "An optional NCBI API key.")
    parser.add_argument("--max-records", type = int, required = True, help = "The ceiling on records to retrieve.")
    parser.add_argument("--batch-size", type = int, required = True, help = "How many records to fetch per efetch request.")
    parser.add_argument("--parse-seqids", action = "store_true", help = "Whether to preserve original sequence IDs in the built database.")
    parser.add_argument("--overwrite", action = "store_true", help = "Whether to rebuild even if a valid database already exists.")
    parser.add_argument("--output-dir", type = Path, required = True, help = "The directory to write this step's own result JSON and setup metadata to.")
    parser.add_argument("--log-file", type = Path, required = True, help = "The path to the log file.")
    parser.add_argument("--global-log-file", type = Path, required = True, help = "The path to the global log file.")

    return parser.parse_args()

def _write_result(output_dir: Path, success: bool, database_prefix: Path, reason: str) -> None:
    """Records this setup step's outcome for every downstream BLAST stage to read."""
    output_dir.mkdir(parents = True, exist_ok = True)
    path = output_dir / "blast_database.json"
    with path.open("w", encoding = "utf-8") as handle:
        json.dump({
            "success": success, "database_prefix": str(database_prefix), "reason": reason,
        }, handle, indent = 4)

def main() -> int:
    """
    Ensures the local BLAST reference database is downloaded, built,
    and validated, recording the outcome for downstream BLAST stages
    to read.

    Returns:
        int: Always 0 -- a download/build failure is recorded, not raised; see module docstring.
    """
    args = parse_args()

    logger = make_logger(
        name = "blast_database_setup",
        log_file_path = args.log_file,
        global_log_file_path = args.global_log_file,
    )

    metadata_file = args.output_dir / "blast_database_metadata.json"

    try:
        return_code = setup_reference_database(
            reference_fasta = args.reference_fasta, database_prefix = args.database_prefix,
            database_title = args.database_title, metadata_file = metadata_file,
            email = args.entrez_email, logger = logger, taxon_scope = args.taxon_scope,
            api_key = args.entrez_api_key, max_records = args.max_records, batch_size = args.batch_size,
            parse_seqids = args.parse_seqids, overwrite = args.overwrite,
        )
    except Exception as e:
        logger.warning(f"BLAST database setup failed unexpectedly: {e}. BLAST evidence will be unavailable this run.")
        _write_result(args.output_dir, success = False, database_prefix = args.database_prefix, reason = f"unexpected error: {e}")
        return 0

    if return_code != 0:
        logger.warning("BLAST database setup failed; BLAST evidence will be unavailable this run.")
        _write_result(args.output_dir, success = False, database_prefix = args.database_prefix, reason = "download or build failed")
        return 0

    _write_result(args.output_dir, success = True, database_prefix = args.database_prefix, reason = "database ready")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
