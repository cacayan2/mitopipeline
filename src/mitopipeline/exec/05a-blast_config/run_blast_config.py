"""
run_blast_config.py

Execution layer that resolves one sample's BLAST configuration
deterministically -- query FASTA, database, output format, and the
output path the BLAST exec stage will write to -- and records
it as a JSON manifest that stage reads next.

Doesn't run BLAST itself; that's the exec stage's job. This step
exists so a bad configuration (an unsupported mode, a missing query
FASTA) is caught and logged before a real BLAST invocation starts.
"""

# Imports
import argparse
import json
from pathlib import Path
from mitopipeline.api.blast import BLAST_OUTPUT_COLUMNS, DEFAULT_OUTPUT_FORMAT, VALID_BLAST_MODES
from mitopipeline.logging.logger_factory import make_logger

def parse_args() -> argparse.Namespace:
    """
    Parses the command line arguments for this script.

    Returns:
        argparse.Namespace: The parsed command line arguments.
    """
    parser = argparse.ArgumentParser(description = "Resolves one sample's BLAST configuration.")

    parser.add_argument("--sample-id", type = str, required = True, help = "The sample id.")
    parser.add_argument("--query-fasta", type = Path, required = True, help = "The assembled mitogenome FASTA to query.")
    parser.add_argument("--database", type = str, required = True, help = "The BLAST database name (local db name, or a remote NCBI db).")
    parser.add_argument("--mode", type = str, default = "local", help = "Either 'local' or 'remote'.")
    parser.add_argument("--blast-output-dir", type = Path, required = True, help = "The directory the BLAST exec stage will write this sample's tabular output to.")
    parser.add_argument("--output-dir", type = Path, required = True, help = "The directory to write this config step's own output to.")
    parser.add_argument("--log-file", type = Path, required = True, help = "The path to the log file.")
    parser.add_argument("--global-log-file", type = Path, required = True, help = "The path to the global log file.")

    return parser.parse_args()

def resolve_config(sample_id: str, query_fasta: Path, database: str, mode: str,
                    blast_output_dir: Path) -> dict:
    """
    Resolves one sample's deterministic BLAST configuration. Pure
    function (no I/O beyond what's already in its arguments) -- the
    same inputs always resolve to the same configuration.

    Args:
        sample_id (str): The sample id.
        query_fasta (Path): The assembled mitogenome FASTA to query.
        database (str): The requested BLAST database.
        mode (str): The requested BLAST mode.
        blast_output_dir (Path): Where the BLAST exec stage will write this sample's output.

    Returns:
        dict: The resolved configuration, JSON-serializable.

    Raises:
        ValueError: If mode isn't one BlastRunner actually supports, or query_fasta doesn't exist.
        FileNotFoundError: If query_fasta doesn't exist.
    """
    if mode not in VALID_BLAST_MODES:
        raise ValueError(f"Unsupported BLAST mode: {mode!r}. Must be one of {sorted(VALID_BLAST_MODES)}.")

    if not query_fasta.is_file():
        raise FileNotFoundError(f"Query FASTA does not exist: {query_fasta}.")

    if not database or not database.strip():
        raise ValueError("No BLAST database given.")

    output_file = blast_output_dir / f"{sample_id}.blast.tsv"

    return {
        "sample_id": sample_id,
        "query_fasta": str(query_fasta),
        "database": database,
        "mode": mode,
        "output_format": DEFAULT_OUTPUT_FORMAT,
        "output_columns": BLAST_OUTPUT_COLUMNS,
        "output_file": str(output_file),
    }

def main() -> int:
    """
    Resolves and records one sample's BLAST configuration.

    Returns:
        int: The return code for the process.
    """
    args = parse_args()

    logger = make_logger(
        name = "blast_config",
        log_file_path = args.log_file,
        global_log_file_path = args.global_log_file,
    )
    context = f"{{{args.sample_id}}}"

    try:
        config = resolve_config(
            sample_id = args.sample_id,
            query_fasta = args.query_fasta,
            database = args.database,
            mode = args.mode,
            blast_output_dir = args.blast_output_dir,
        )
    except (ValueError, FileNotFoundError) as e:
        logger.error(f"{context} Could not resolve BLAST configuration: {e}")
        return 1

    args.output_dir.mkdir(parents = True, exist_ok = True)
    config_path = args.output_dir / f"{args.sample_id}.blast_config.json"
    with config_path.open("w", encoding = "utf-8") as handle:
        json.dump(config, handle, indent = 4)

    logger.info(f"{context} Resolved BLAST config: database={config['database']}, mode={config['mode']}. Written to {config_path}.")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
