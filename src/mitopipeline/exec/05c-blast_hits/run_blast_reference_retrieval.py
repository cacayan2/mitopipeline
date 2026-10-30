"""
run_blast_reference_retrieval.py

Execution layer for retrieving the selected BLAST hit's reference
FASTA: reads the top-hit TSV blast_hits.write_top_blast_hits()
already wrote, retrieves rank 1's accession from the local BLAST
database with BlastDbCmdRunner (BaseTool), and records the result
(or an explicit, non-crashing failure for an accession the database
doesn't actually have) with provenance for later stages to read.
"""

# Imports
import argparse
import csv
import json
from datetime import datetime, timezone
from pathlib import Path
from mitopipeline.api.blastdbcmd import BlastDbCmdRunner
from mitopipeline.logging.logger_factory import make_logger

def parse_args() -> argparse.Namespace:
    """
    Parses the command line arguments for this script.

    Returns:
        argparse.Namespace: The parsed command line arguments.
    """
    parser = argparse.ArgumentParser(description = "Retrieves the top-ranked selected BLAST hit's reference FASTA from a local database.")

    parser.add_argument("--sample-id", type = str, required = True, help = "The sample id.")
    parser.add_argument("--top-hits-tsv", type = Path, required = True, help = "The top-hit TSV written by blast_hits.write_top_blast_hits().")
    parser.add_argument("--database", type = str, required = True, help = "The local BLAST database name to retrieve from.")
    parser.add_argument("--output-dir", type = Path, required = True, help = "The directory to write the retrieved FASTA and provenance record to.")
    parser.add_argument("--working-dir", type = Path, required = True, help = "The working directory.")
    parser.add_argument("--log-file", type = Path, required = True, help = "The path to the log file.")
    parser.add_argument("--global-log-file", type = Path, required = True, help = "The path to the global log file.")
    parser.add_argument("--metrics-dir", type = Path, default = None, help = "Where to persist this run's compute metrics as JSON. Omit to skip persistence.")

    return parser.parse_args()

def _read_top_hit_accession(top_hits_tsv: Path, sample_id: str) -> str | None:
    """
    Reads the rank-1 (strongest) hit's accession for this sample from
    a top-hit TSV.

    Args:
        top_hits_tsv (Path): The TSV to read.
        sample_id (str): The sample whose rank-1 hit to select (a shared TSV may hold multiple samples).

    Returns:
        str | None: The rank-1 hit's accession (sseqid), or None if this sample has no rank-1 row.
    """
    with top_hits_tsv.open("r", encoding = "utf-8", newline = "") as handle:
        for row in csv.DictReader(handle, delimiter = "\t"):
            if row.get("sample_id") == sample_id and row.get("rank") == "1":
                return row["sseqid"]
    return None

def _write_record(output_dir: Path, sample_id: str, retrieved: bool, database: str,
                   accession: str | None = None, fasta_path: Path | None = None, reason: str = "") -> None:
    """
    Records this retrieval's outcome and provenance for later stages
    to read without re-deriving it.

    Args:
        output_dir (Path): Where to write the record.
        sample_id (str): The sample id.
        retrieved (bool): Whether a FASTA was actually retrieved.
        database (str): The local BLAST database retrieved from.
        accession (str | None): The accession retrieved, if any.
        fasta_path (Path | None): The retrieved FASTA's path, if any.
        reason (str): A short, human-readable explanation.

    Returns:
        None
    """
    output_dir.mkdir(parents = True, exist_ok = True)
    path = output_dir / f"{sample_id}.blast_reference_retrieval.json"
    with path.open("w", encoding = "utf-8") as handle:
        json.dump({
            "sample_id": sample_id,
            "retrieved": retrieved,
            "accession": accession,
            "database": database,
            "fasta_path": str(fasta_path) if fasta_path else None,
            "retrieved_at": datetime.now(timezone.utc).isoformat(),
            "reason": reason,
        }, handle, indent = 4)

def main() -> int:
    """
    Retrieves the selected BLAST hit's reference FASTA for one sample.

    Returns:
        int: The return code for the process.
    """
    args = parse_args()

    logger = make_logger(
        name = "blast_reference_retrieval",
        log_file_path = args.log_file,
        global_log_file_path = args.global_log_file,
    )
    context = f"{{{args.sample_id}}}"

    if not args.top_hits_tsv.is_file():
        logger.error(f"{context} Top-hit TSV not found: {args.top_hits_tsv}.")
        _write_record(args.output_dir, args.sample_id, retrieved = False, database = args.database,
                     reason = f"top-hit TSV not found: {args.top_hits_tsv}")
        return 1

    accession = _read_top_hit_accession(args.top_hits_tsv, args.sample_id)
    if accession is None:
        logger.warning(f"{context} No rank-1 BLAST hit found in {args.top_hits_tsv} for this sample; nothing to retrieve.")
        _write_record(args.output_dir, args.sample_id, retrieved = False, database = args.database,
                     reason = "no rank-1 hit for this sample")
        return 1

    runner = BlastDbCmdRunner(
        working_dir = args.working_dir,
        output_dir = args.output_dir,
        sample_id = args.sample_id,
        accession = accession,
        database = args.database,
        logger = logger,
        metrics_dir = args.metrics_dir,
    )

    try:
        result = runner.run()
    except SystemExit as e:
        logger.warning(f"{context} Could not retrieve accession '{accession}' (exit {e.code}) -- "
                       "likely not present in this database.")
        _write_record(args.output_dir, args.sample_id, retrieved = False, database = args.database,
                     accession = accession, reason = f"accession not found in database '{args.database}'")
        return 1

    if not result.success:
        logger.error(f"{context} blastdbcmd exited with return code {result.return_code}.")
        _write_record(args.output_dir, args.sample_id, retrieved = False, database = args.database,
                     accession = accession, reason = f"blastdbcmd exited {result.return_code}")
        return result.return_code or 1

    logger.info(f"{context} Retrieved accession '{accession}' from '{args.database}' to {runner.output_fasta}.")
    _write_record(args.output_dir, args.sample_id, retrieved = True, database = args.database,
                 accession = accession, fasta_path = runner.output_fasta, reason = "retrieved successfully")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
