"""
run_sample_report.py

Execution layer for the per-sample report (#358/#359/#361): collects
one sample's real result JSON into a SampleReportData, renders the
Markdown report (including the live methods section), and exports it
to PDF. A missing PDF export (pandoc/typst unavailable) never blocks
the Markdown report itself from being recorded.
"""

# Imports
import argparse
import json
from pathlib import Path
from mitopipeline.logging.logger_factory import make_logger
from mitopipeline.utils.reporting.collect import collect_sample_report_data
from mitopipeline.utils.reporting.markdown_report import write_sample_report
from mitopipeline.utils.reporting.pdf_export import export_markdown_to_pdf

def parse_args() -> argparse.Namespace:
    """
    Parses the command line arguments for this script.

    Returns:
        argparse.Namespace: The parsed command line arguments.
    """
    parser = argparse.ArgumentParser(description = "Renders one sample's Markdown + PDF report.")

    parser.add_argument("--sample-id", type = str, required = True, help = "The sample id.")
    parser.add_argument("--runtime-manifest", type = Path, required = True, help = "The runtime manifest TSV.")
    parser.add_argument("--statistics-json", type = Path, default = None, help = "The merged statistics JSON written by run_stats_aggregation.py.")
    parser.add_argument("--biological-validation-json", type = Path, default = None, help = "The validation result JSON written by run_biological_validation.py.")
    parser.add_argument("--visualization-result-json", type = Path, default = None, help = "The result JSON written by run_circular_map.py.")
    parser.add_argument("--fastp-json", type = Path, default = None, help = "This sample's own fastp JSON report.")
    parser.add_argument("--getorganelle-config-json", type = Path, default = None, help = "The resolved GetOrganelle config JSON.")
    parser.add_argument("--blast-config-json", type = Path, default = None, help = "The resolved BLAST config JSON.")
    parser.add_argument("--mitos2-config-json", type = Path, default = None, help = "The resolved MITOS2 config JSON.")
    parser.add_argument("--sra-metadata-tsv", type = Path, default = None, help = "The SRA metadata table written by run_sra_metadata.py.")
    parser.add_argument("--biosample-metadata-tsv", type = Path, default = None, help = "The BioSample metadata table written by run_biosample_metadata.py.")
    parser.add_argument("--genbank-record-json", type = Path, default = None, help = "This sample's own GenBank record JSON written by run_mitogenome_metadata.py.")
    parser.add_argument("--pandoc-bin", type = str, default = "pandoc", help = "The pandoc executable name. Defaults to 'pandoc'.")
    parser.add_argument("--pdf-engine", type = str, default = "typst", help = "The PDF engine pandoc should use. Defaults to 'typst'.")
    parser.add_argument("--output-dir", type = Path, required = True, help = "The directory to write this sample's report to.")
    parser.add_argument("--log-file", type = Path, required = True, help = "The path to the log file.")
    parser.add_argument("--global-log-file", type = Path, required = True, help = "The path to the global log file.")

    return parser.parse_args()

def _write_result(output_dir: Path, sample_id: str, success: bool, markdown_path: Path | None,
                   pdf_path: Path | None, pdf_success: bool, reason: str) -> None:
    """Records this sample's report outcome for the run-level report to read next."""
    output_dir.mkdir(parents = True, exist_ok = True)
    path = output_dir / f"{sample_id}.report_result.json"
    with path.open("w", encoding = "utf-8") as handle:
        json.dump({
            "sample_id": sample_id, "success": success,
            "markdown_path": str(markdown_path) if markdown_path else None,
            "pdf_path": str(pdf_path) if (pdf_path and pdf_success) else None,
            "reason": reason,
        }, handle, indent = 4)

def main() -> int:
    """
    Renders one sample's Markdown + PDF report.

    Returns:
        int: The return code for the process.
    """
    args = parse_args()

    logger = make_logger(
        name = "sample_report",
        log_file_path = args.log_file,
        global_log_file_path = args.global_log_file,
    )
    context = f"{{{args.sample_id}}}"

    report_data = collect_sample_report_data(
        args.sample_id, args.runtime_manifest, statistics_json = args.statistics_json,
        biological_validation_json = args.biological_validation_json,
        visualization_result_json = args.visualization_result_json, fastp_json = args.fastp_json,
        getorganelle_config_json = args.getorganelle_config_json, blast_config_json = args.blast_config_json,
        mitos2_config_json = args.mitos2_config_json, sra_metadata_tsv = args.sra_metadata_tsv,
        biosample_metadata_tsv = args.biosample_metadata_tsv, genbank_record_json = args.genbank_record_json,
        logger = logger,
    )

    markdown_path = args.output_dir / f"{args.sample_id}.report.md"
    try:
        write_sample_report(report_data, markdown_path)
    except OSError as e:
        logger.error(f"{context} Could not write Markdown report: {e}")
        _write_result(args.output_dir, args.sample_id, success = False, markdown_path = None,
                      pdf_path = None, pdf_success = False, reason = f"markdown write failed: {e}")
        return 1

    pdf_path, pdf_success, reason = None, False, "completed successfully"
    try:
        result = export_markdown_to_pdf(markdown_path, pandoc_bin = args.pandoc_bin, pdf_engine = args.pdf_engine, logger = logger)
        pdf_path, pdf_success = result.pdf_path, result.success
        if not pdf_success:
            reason = f"pdf export failed: {result.stderr.strip()[:200]}"
    except FileNotFoundError as e:
        logger.warning(f"{context} PDF export skipped: {e}")
        reason = f"pdf export skipped: {e}"

    logger.info(f"{context} Sample report written to {markdown_path} (PDF: {'yes' if pdf_success else 'no'}).")
    _write_result(args.output_dir, args.sample_id, success = True, markdown_path = markdown_path,
                 pdf_path = pdf_path, pdf_success = pdf_success, reason = reason)
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
