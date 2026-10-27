"""
run_run_report.py

Execution layer for the run-level report (#360/#441/#443): collects
every sample's own report data plus the batch-level phylogeny
outcome, renders the run-level Markdown report, and exports it to
PDF.
"""

# Imports
import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from mitopipeline.logging.logger_factory import make_logger
from mitopipeline.utils.reporting.collect import collect_phylogeny_summary, collect_run_report_data, collect_sample_report_data
from mitopipeline.utils.reporting.pdf_export import export_markdown_to_pdf
from mitopipeline.utils.reporting.run_report import write_run_report

def parse_args() -> argparse.Namespace:
    """
    Parses the command line arguments for this script.

    Returns:
        argparse.Namespace: The parsed command line arguments.
    """
    parser = argparse.ArgumentParser(description = "Renders the run-level Markdown + PDF report.")

    parser.add_argument("--job-id", type = str, required = True, help = "The pipeline job id.")
    parser.add_argument("--sample-id", type = str, required = True, action = "append", help = "A sample id to include. Repeat for every sample in the run.")
    parser.add_argument("--runtime-manifest", type = Path, required = True, help = "The runtime manifest TSV.")
    parser.add_argument("--job-dir", type = Path, required = True, help = "The pipeline job's own output directory, used to resolve every per-sample result JSON by its known stage-relative path.")
    parser.add_argument("--iqtree-result-json", type = Path, default = None, help = "The batch-level IQ-TREE result JSON, if phylogeny ran.")
    parser.add_argument("--phylogeny-model", type = str, default = None, help = "The substitution model IQ-TREE was run with, if known.")
    parser.add_argument("--pandoc-bin", type = str, default = "pandoc", help = "The pandoc executable name. Defaults to 'pandoc'.")
    parser.add_argument("--pdf-engine", type = str, default = "typst", help = "The PDF engine pandoc should use. Defaults to 'typst'.")
    parser.add_argument("--output-dir", type = Path, required = True, help = "The directory to write the run-level report to.")
    parser.add_argument("--log-file", type = Path, required = True, help = "The path to the log file.")
    parser.add_argument("--global-log-file", type = Path, required = True, help = "The path to the global log file.")

    return parser.parse_args()

def _sample_paths(job_dir: Path, sample_id: str) -> dict:
    """The same real, fixed stage-relative paths every other stage's own rule already uses."""
    return {
        "statistics_json": job_dir / "12-stats_aggregation" / sample_id / f"{sample_id}.statistics.json",
        "biological_validation_json": job_dir / "06-contig_selection" / sample_id / f"{sample_id}.biological_validation.json",
        "visualization_result_json": job_dir / "09-visualization" / sample_id / f"{sample_id}.circular_map_result.json",
        "fastp_json": job_dir / "trimming" / sample_id / f"{sample_id}.fastp.json",
        "getorganelle_config_json": job_dir / "04a-getorganelle_config" / sample_id / f"{sample_id}.getorganelle_config.json",
        "blast_config_json": job_dir / "05a-blast_config" / sample_id / f"{sample_id}.blast_config.json",
        "mitos2_config_json": job_dir / "08a-mitos2_config" / sample_id / f"{sample_id}.mitos2_config.json",
        "sra_metadata_tsv": job_dir / "14a-sra_archival" / "sra_metadata.tsv",
        "biosample_metadata_tsv": job_dir / "14b-biosample_archival" / "biosample_metadata.tsv",
        "genbank_record_json": job_dir / "14c-mitogenome_archival" / f"{sample_id}.genbank_record.json",
    }

def main() -> int:
    """
    Renders the run-level Markdown + PDF report.

    Returns:
        int: The return code for the process.
    """
    args = parse_args()

    logger = make_logger(
        name = "run_report",
        log_file_path = args.log_file,
        global_log_file_path = args.global_log_file,
    )

    samples = [
        collect_sample_report_data(sample_id, args.runtime_manifest, logger = logger,
                                   **_sample_paths(args.job_dir, sample_id))
        for sample_id in args.sample_id
    ]
    phylogeny = collect_phylogeny_summary(args.iqtree_result_json, model = args.phylogeny_model)
    generated_at = datetime.now(timezone.utc).isoformat(timespec = "seconds")
    run_data = collect_run_report_data(args.job_id, generated_at, samples, phylogeny = phylogeny)

    markdown_path = args.output_dir / f"{args.job_id}.run_report.md"
    write_run_report(run_data, markdown_path)

    pdf_success = False
    try:
        result = export_markdown_to_pdf(markdown_path, pandoc_bin = args.pandoc_bin, pdf_engine = args.pdf_engine, logger = logger)
        pdf_success = result.success
    except FileNotFoundError as e:
        logger.warning(f"PDF export skipped: {e}")

    logger.info(f"Run report written to {markdown_path} (PDF: {'yes' if pdf_success else 'no'}), {len(samples)} sample(s).")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
