"""
run_build_sample_dataset.py

Execution layer for building one sample's own real phylogenetic
dataset (restoring v0.1.0's per-sample design): its assembled genome
plus the real sequences for its own top-ranked BLAST reference hits,
fetched live from NCBI Entrez. Unlike the earlier batch-level combined
FASTA this replaces, a sample's own tree never depends on how many
other samples are in the same run.

A sample with no selected assembly, or no BLAST hits, is a normal,
non-error outcome this script must still complete for -- the report
simply shows no phylogeny section for that sample.
"""

# Imports
import argparse
import json
from pathlib import Path
from mitopipeline.logging.logger_factory import make_logger
from mitopipeline.utils.phylogeny.sample_dataset import build_sample_phylogeny_dataset

def parse_args() -> argparse.Namespace:
    """
    Parses the command line arguments for this script.

    Returns:
        argparse.Namespace: The parsed command line arguments.
    """
    parser = argparse.ArgumentParser(description = "Builds one sample's own real phylogenetic dataset (assembled genome + live-fetched BLAST references).")

    parser.add_argument("--sample-id", type = str, required = True, help = "The sample id.")
    parser.add_argument("--assembly-fasta", type = Path, required = True, help = "The sample's own selected/resolved assembly FASTA.")
    parser.add_argument("--top-hits-tsv", type = Path, required = True, help = "The sample's own ranked BLAST top-hit TSV.")
    parser.add_argument("--output-fasta", type = Path, required = True, help = "Where to write the combined FASTA.")
    parser.add_argument("--entrez-email", type = str, required = True, help = "The contact email NCBI requires for Entrez requests.")
    parser.add_argument("--entrez-api-key", type = str, default = None, help = "An optional NCBI API key.")
    parser.add_argument("--output-dir", type = Path, required = True, help = "The directory to write this step's own result JSON to.")
    parser.add_argument("--log-file", type = Path, required = True, help = "The path to the log file.")
    parser.add_argument("--global-log-file", type = Path, required = True, help = "The path to the global log file.")

    return parser.parse_args()

def _write_result(output_dir: Path, sample_id: str, success: bool, output_fasta: Path | None,
                   reference_count: int | None, reason: str) -> None:
    output_dir.mkdir(parents = True, exist_ok = True)
    with (output_dir / f"{sample_id}.phylogeny_dataset_result.json").open("w", encoding = "utf-8") as handle:
        json.dump({
            "sample_id": sample_id, "success": success,
            "output_fasta": str(output_fasta) if (output_fasta and success) else None,
            "reference_count": reference_count, "reason": reason,
        }, handle, indent = 4)

def main() -> int:
    """
    Builds one sample's own real phylogenetic dataset and records the
    outcome.

    Returns:
        int: Always 0 -- a missing assembly/hits or a live NCBI failure is recorded, not raised; see module docstring.
    """
    args = parse_args()

    logger = make_logger(
        name = "phylogeny_dataset",
        log_file_path = args.log_file,
        global_log_file_path = args.global_log_file,
    )
    context = f"{{{args.sample_id}}}"

    if not args.assembly_fasta.is_file():
        logger.info(f"{context} No assembly FASTA at {args.assembly_fasta}; skipping phylogeny dataset.")
        _write_result(args.output_dir, args.sample_id, success = False, output_fasta = None,
                      reference_count = None, reason = "no selected assembly")
        return 0

    if not args.top_hits_tsv.is_file():
        logger.info(f"{context} No BLAST top-hit TSV at {args.top_hits_tsv}; skipping phylogeny dataset.")
        _write_result(args.output_dir, args.sample_id, success = False, output_fasta = None,
                      reference_count = None, reason = "no BLAST hits available")
        return 0

    try:
        reference_count = build_sample_phylogeny_dataset(
            sample_id = args.sample_id, assembly_fasta = args.assembly_fasta, top_hits_tsv = args.top_hits_tsv,
            output_fasta = args.output_fasta, entrez_email = args.entrez_email, entrez_api_key = args.entrez_api_key,
            logger = logger,
        )
    except Exception as e:
        logger.warning(f"{context} Could not build phylogeny dataset: {e}.")
        _write_result(args.output_dir, args.sample_id, success = False, output_fasta = None,
                      reference_count = None, reason = f"could not build dataset: {e}")
        return 0

    logger.info(f"{context} Phylogeny dataset ready: assembled genome + {reference_count} reference sequence(s).")
    _write_result(args.output_dir, args.sample_id, success = True, output_fasta = args.output_fasta,
                 reference_count = reference_count, reason = "completed successfully")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
