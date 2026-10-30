"""
run_seed_selection.py

Execution layer for seed selection (the seed-selection module's own logic,
never wired to a CLI until now): selects one sample's mitochondrial
seed sequence via select_seed() and records the outcome -- the
provenance (locus, accession, fallback reason) the assembly stage's
own config step reads next.
"""

# Imports
import argparse
import json
from pathlib import Path
from mitopipeline.logging.logger_factory import make_logger
from mitopipeline.utils.seed.seed_selector import select_seed

def parse_args() -> argparse.Namespace:
    """
    Parses the command line arguments for this script.

    Returns:
        argparse.Namespace: The parsed command line arguments.
    """
    parser = argparse.ArgumentParser(description = "Selects one sample's mitochondrial seed sequence.")

    parser.add_argument("--sample-id", type = str, required = True, help = "The sample id.")
    parser.add_argument("--scientific-name", type = str, required = True, help = "The sample's scientific name.")
    parser.add_argument("--database-dir", type = Path, required = True, help = "The local MitoFish reference database directory.")
    parser.add_argument("--output-dir", type = Path, required = True, help = "The directory to write the selected seed FASTA (and this step's own record) to.")
    parser.add_argument("--log-file", type = Path, required = True, help = "The path to the log file.")
    parser.add_argument("--global-log-file", type = Path, required = True, help = "The path to the global log file.")

    return parser.parse_args()

def main() -> int:
    """
    Selects and records one sample's mitochondrial seed sequence.

    Returns:
        int: The return code for the process (always 0 -- select_seed() never raises; a fallback result is a valid, non-error outcome).
    """
    args = parse_args()

    logger = make_logger(
        name = "seed_selection",
        log_file_path = args.log_file,
        global_log_file_path = args.global_log_file,
    )
    context = f"{{{args.sample_id}}}"

    result = select_seed(
        sample_id = args.sample_id,
        scientific_name = args.scientific_name,
        database_dir = args.database_dir,
        output_dir = args.output_dir,
        logger = logger,
    )

    args.output_dir.mkdir(parents = True, exist_ok = True)
    record_path = args.output_dir / f"{args.sample_id}.seed_selection.json"
    with record_path.open("w", encoding = "utf-8") as handle:
        json.dump({
            "sample_id": result.sample_id,
            "scientific_name": result.scientific_name,
            "gene": result.gene,
            "accession": result.accession,
            "source": result.source,
            "fasta_path": str(result.fasta_path) if result.fasta_path else None,
            "fallback_used": result.fallback_used,
            "fallback_reason": result.fallback_reason,
            "match_level": result.match_level,
            "checksum": result.checksum,
            "ncbi_taxon_id": result.ncbi_taxon_id,
            "sequence_length": result.sequence_length,
        }, handle, indent = 4)

    if result.fallback_used:
        logger.info(f"{context} No seed selected ({result.fallback_reason}); GetOrganelle will use its own default database.")
    else:
        logger.info(f"{context} Selected seed: locus={result.gene}, accession={result.accession}, {result.sequence_length}bp.")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
