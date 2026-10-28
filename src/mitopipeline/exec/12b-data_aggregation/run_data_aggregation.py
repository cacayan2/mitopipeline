"""
run_data_aggregation.py

Execution layer for consolidating every sample's own final assembly/
annotation/phylogeny outputs into one flat aggregated_data/ tree
(Emil's own direct ask, 2026-10-28) -- assemblies/, annotations/
(+ annotations/gff/), and phylogenies/, so a human reviewer never has
to go hunting across a dozen per-stage directories for the files that
actually matter.

Also consolidates the run's real NCBI archival materials into
archival/ (Emil's own direct follow-up ask) -- SRA/BioSample tables,
each sample's own real GenBank submission folder, and the validation
stage's own submission checklist/readiness report/checksums.
"""

# Imports
import argparse
import json
from pathlib import Path
from mitopipeline.logging.logger_factory import make_logger
from mitopipeline.utils.aggregation.data_aggregation import aggregate_archival_data, aggregate_sample_data

def parse_args() -> argparse.Namespace:
    """
    Parses the command line arguments for this script.

    Returns:
        argparse.Namespace: The parsed command line arguments.
    """
    parser = argparse.ArgumentParser(description = "Consolidates every sample's own final outputs into aggregated_data/.")

    parser.add_argument("--job-id", type = str, required = True, help = "The pipeline job id.")
    parser.add_argument("--sample-id", type = str, required = True, action = "append", help = "A sample id to include. Repeat for every sample in the run.")
    parser.add_argument("--job-dir", type = Path, required = True, help = "The pipeline job's own output directory, used to resolve every per-sample source file by its known stage-relative path.")
    parser.add_argument("--output-dir", type = Path, required = True, help = "The aggregated_data/ directory itself.")
    parser.add_argument("--sra-metadata-tsv", type = Path, default = None, help = "The run's own real SRA metadata table, written by run_sra_metadata.py.")
    parser.add_argument("--biosample-metadata-tsv", type = Path, default = None, help = "The run's own real BioSample metadata table, written by run_biosample_metadata.py.")
    parser.add_argument("--genbank-metadata-tsv", type = Path, default = None, help = "The run's own real combined GenBank metadata table, written by run_archival_validation.py.")
    parser.add_argument("--archival-validation-tsv", type = Path, default = None, help = "The run's own real per-sample submission-readiness report, written by run_archival_validation.py.")
    parser.add_argument("--submission-checklist-md", type = Path, default = None, help = "The run's own real human-readable submission checklist, written by run_archival_validation.py.")
    parser.add_argument("--checksums-tsv", type = Path, default = None, help = "The run's own real MD5 checksum manifest, written by run_archival_validation.py.")
    parser.add_argument("--log-file", type = Path, required = True, help = "The path to the log file.")
    parser.add_argument("--global-log-file", type = Path, required = True, help = "The path to the global log file.")

    return parser.parse_args()

def _read_json(path: Path) -> dict | None:
    """Reads a JSON file, returning None when it's absent or unreadable."""
    if not path.is_file():
        return None
    try:
        return json.loads(path.read_text(encoding = "utf-8"))
    except (OSError, json.JSONDecodeError):
        return None

def _selected_assembly_fasta(job_dir: Path, sample_id: str) -> Path | None:
    """Reads a sample's own selected assembly FASTA path from its real contig_selection.json -- the same selected attempt (primary or any rescue) every other stage's own report/stats derivation already uses."""
    selection = _read_json(job_dir / "06-contig_selection" / sample_id / f"{sample_id}.contig_selection.json")
    if selection is None or not selection.get("selected") or not selection.get("fasta_path"):
        return None
    return Path(selection["fasta_path"])

def _annotated_map_png(job_dir: Path, sample_id: str) -> Path | None:
    """Reads a sample's own real annotated genome map PNG path from its real circular_map_result.json."""
    result = _read_json(job_dir / "09-visualization" / sample_id / f"{sample_id}.circular_map_result.json")
    if result is None or not result.get("success") or not result.get("output_png"):
        return None
    return Path(result["output_png"])

def _annotation_gff(job_dir: Path, sample_id: str) -> Path | None:
    """Reads a sample's own real MITOS2 annotation GFF path from its real mitos2_annotation_result.json -- never a hardcoded filename, since Mitos2Runner.postprocess_outputs() renames runmitos.py's own fixed "result.gff" to a sample-prefixed name."""
    result = _read_json(job_dir / "08b-annotation" / sample_id / f"{sample_id}.mitos2_annotation_result.json")
    if result is None or not result.get("success") or not result.get("gff_path"):
        return None
    return Path(result["gff_path"])

def _genbank_submission_dir(job_dir: Path, sample_id: str) -> Path | None:
    """Reads a sample's own real GenBank submission folder path from its real genbank_record.json (written flat under 14c-mitogenome_archival/, not nested under the per-sample done-marker subdirectory)."""
    record = _read_json(job_dir / "14c-mitogenome_archival" / f"{sample_id}.genbank_record.json")
    if record is None or not record.get("submission_folder"):
        return None
    return Path(record["submission_folder"])

def main() -> int:
    """
    Consolidates every sample's own final outputs into aggregated_data/
    and records the outcome.

    Returns:
        int: Always 0 -- a sample with nothing to aggregate (no assembly) is a normal outcome, not a failure.
    """
    args = parse_args()

    logger = make_logger(
        name = "data_aggregation",
        log_file_path = args.log_file,
        global_log_file_path = args.global_log_file,
    )

    results = []
    for sample_id in args.sample_id:
        context = f"{{{sample_id}}}"
        result = aggregate_sample_data(
            sample_id, args.output_dir,
            assembly_fasta = _selected_assembly_fasta(args.job_dir, sample_id),
            annotation_image_png = _annotated_map_png(args.job_dir, sample_id),
            annotation_gff = _annotation_gff(args.job_dir, sample_id),
            phylogeny_newick = args.job_dir / "11-phylogeny" / sample_id / f"{sample_id}.treefile.nwk",
            phylogeny_image_png = args.job_dir / "11-phylogeny" / sample_id / "tree_figure.png",
            genbank_submission_dir = _genbank_submission_dir(args.job_dir, sample_id),
            logger = logger,
        )
        logger.info(
            f"{context} Aggregated: assembly={'yes' if result.assembly_path else 'no'}, "
            f"annotation_image={'yes' if result.annotation_image_path else 'no'}, "
            f"annotation_gff={'yes' if result.annotation_gff_path else 'no'}, "
            f"phylogeny_newick={'yes' if result.phylogeny_newick_path else 'no'}, "
            f"phylogeny_image={'yes' if result.phylogeny_image_path else 'no'}, "
            f"genbank_submission_dir={'yes' if result.genbank_submission_dir_path else 'no'}."
        )
        results.append(result)

    archival_result = aggregate_archival_data(
        args.output_dir,
        sra_metadata_tsv = args.sra_metadata_tsv, biosample_metadata_tsv = args.biosample_metadata_tsv,
        genbank_metadata_tsv = args.genbank_metadata_tsv, archival_validation_tsv = args.archival_validation_tsv,
        submission_checklist_md = args.submission_checklist_md, checksums_tsv = args.checksums_tsv,
        logger = logger,
    )
    logger.info(
        f"Archival materials aggregated: sra={'yes' if archival_result.sra_metadata_path else 'no'}, "
        f"biosample={'yes' if archival_result.biosample_metadata_path else 'no'}, "
        f"genbank_metadata={'yes' if archival_result.genbank_metadata_path else 'no'}, "
        f"validation={'yes' if archival_result.archival_validation_path else 'no'}, "
        f"checklist={'yes' if archival_result.submission_checklist_path else 'no'}, "
        f"checksums={'yes' if archival_result.checksums_path else 'no'}."
    )

    args.output_dir.mkdir(parents = True, exist_ok = True)
    output_path = args.output_dir / f"{args.job_id}.data_aggregation_result.json"
    with output_path.open("w", encoding = "utf-8") as handle:
        json.dump({
            "job_id": args.job_id,
            "samples": [
                {
                    "sample_id": r.sample_id,
                    "assembly_path": str(r.assembly_path) if r.assembly_path else None,
                    "annotation_image_path": str(r.annotation_image_path) if r.annotation_image_path else None,
                    "annotation_gff_path": str(r.annotation_gff_path) if r.annotation_gff_path else None,
                    "phylogeny_newick_path": str(r.phylogeny_newick_path) if r.phylogeny_newick_path else None,
                    "phylogeny_image_path": str(r.phylogeny_image_path) if r.phylogeny_image_path else None,
                    "genbank_submission_dir_path": str(r.genbank_submission_dir_path) if r.genbank_submission_dir_path else None,
                }
                for r in results
            ],
            "archival": {
                "sra_metadata_path": str(archival_result.sra_metadata_path) if archival_result.sra_metadata_path else None,
                "biosample_metadata_path": str(archival_result.biosample_metadata_path) if archival_result.biosample_metadata_path else None,
                "genbank_metadata_path": str(archival_result.genbank_metadata_path) if archival_result.genbank_metadata_path else None,
                "archival_validation_path": str(archival_result.archival_validation_path) if archival_result.archival_validation_path else None,
                "submission_checklist_path": str(archival_result.submission_checklist_path) if archival_result.submission_checklist_path else None,
                "checksums_path": str(archival_result.checksums_path) if archival_result.checksums_path else None,
            },
        }, handle, indent = 4)

    logger.info(f"Data aggregation complete for {len(results)} sample(s); written to {args.output_dir}.")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
