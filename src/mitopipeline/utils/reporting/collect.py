"""
collect.py

Assembles report data (#358-361) from the real result JSON every
earlier stage already wrote -- #356's own merged statistics, #414's
biological validation record, #349-352's visualization result, and
each stage's own resolved configuration (for the live methods
section) -- never re-deriving a statistic a typed artifact already has.
"""

# Imports
import csv
import json
from logging import Logger
from pathlib import Path
from mitopipeline.utils.reporting.models import (
    AnnotationConfigSummary,
    ArchivalSummary,
    AssemblyConfigSummary,
    BlastConfigSummary,
    PhylogenySummary,
    RunReportData,
    SampleReportData,
)
from mitopipeline.utils.statistics.models import AnnotationStatistics, AssemblyStatistics, BlastStatistics

def _read_json(path: Path | None) -> dict | None:
    """Reads a JSON file, returning None when it's absent, unset, or unreadable."""
    if path is None or not Path(path).is_file():
        return None
    try:
        return json.loads(Path(path).read_text(encoding = "utf-8"))
    except (OSError, json.JSONDecodeError):
        return None

def read_sample_metadata(runtime_manifest_path: str | Path, sample_id: str,
                          logger: Logger | None = None) -> dict[str, str]:
    """
    Reads one sample's own row from the runtime manifest.

    Args:
        runtime_manifest_path (str | Path): The runtime manifest TSV.
        sample_id (str): The sample id to look up.
        logger (Logger | None): Logger to use. Defaults to None.

    Returns:
        dict[str, str]: This sample's metadata, or {} if the sample or manifest is missing.
    """
    runtime_manifest_path = Path(runtime_manifest_path)
    if not runtime_manifest_path.is_file():
        if logger is not None: logger.warning(f"Runtime manifest not found: {runtime_manifest_path}.")
        return {}

    with runtime_manifest_path.open("r", encoding = "utf-8", newline = "") as handle:
        for row in csv.DictReader(handle, delimiter = "\t"):
            if row.get("sample_id", "").strip() == sample_id:
                return {k: v.strip() for k, v in row.items() if k and v and v.strip()}

    if logger is not None: logger.warning(f"Sample {sample_id} not found in runtime manifest {runtime_manifest_path}.")
    return {}

def parse_fastp_stats(fastp_json_path: str | Path | None) -> dict | None:
    """
    Parses fastp's own real JSON report into the summary fields the
    report displays.

    Args:
        fastp_json_path (str | Path | None): fastp's own JSON report path.

    Returns:
        dict | None: The parsed summary, or None when unavailable/malformed.
    """
    data = _read_json(Path(fastp_json_path) if fastp_json_path else None)
    if data is None:
        return None
    try:
        before, after = data["summary"]["before_filtering"], data["summary"]["after_filtering"]
        filtering = data.get("filtering_result", {})
        reads_in, reads_out = int(before["total_reads"]), int(after["total_reads"])
        return {
            "reads_in": reads_in, "reads_out": reads_out,
            "read_retention_rate": (reads_out / reads_in) if reads_in > 0 else 0.0,
            "q20_rate_out": after["q20_rate"], "q30_rate_out": after["q30_rate"],
            "gc_content_out": after["gc_content"],
            "low_quality_reads": int(filtering.get("low_quality_reads", 0)),
        }
    except (KeyError, TypeError, ValueError):
        return None

def _read_tsv_row(path: Path | None, key_column: str, key_value: str) -> dict[str, str] | None:
    """Reads one row from a TSV, matched by key_column == key_value, or None when absent."""
    if path is None or not Path(path).is_file():
        return None
    with Path(path).open("r", encoding = "utf-8", newline = "") as handle:
        for row in csv.DictReader(handle, delimiter = "\t"):
            if row.get(key_column) == key_value:
                return row
    return None

def collect_archival_summary(sample_id: str, sra_metadata_tsv: str | Path | None = None,
                              biosample_metadata_tsv: str | Path | None = None,
                              genbank_record_json: str | Path | None = None) -> ArchivalSummary | None:
    """
    Reads one sample's real NCBI archival preparation status from the
    SRA/BioSample tables and its own GenBank submission record.

    Args:
        sample_id (str): The sample id.
        sra_metadata_tsv (str | Path | None): The SRA metadata table written by run_sra_metadata.py.
        biosample_metadata_tsv (str | Path | None): The BioSample metadata table written by run_biosample_metadata.py.
        genbank_record_json (str | Path | None): This sample's own GenBank record JSON written by run_mitogenome_metadata.py.

    Returns:
        ArchivalSummary | None: The summary, or None when no archival source was supplied at all.
    """
    if sra_metadata_tsv is None and biosample_metadata_tsv is None and genbank_record_json is None:
        return None

    sra_row = _read_tsv_row(Path(sra_metadata_tsv) if sra_metadata_tsv else None, "sample_name", sample_id)
    biosample_row = _read_tsv_row(Path(biosample_metadata_tsv) if biosample_metadata_tsv else None, "sample_name", sample_id)
    genbank_record = _read_json(Path(genbank_record_json) if genbank_record_json else None)

    return ArchivalSummary(
        sra_review_status = sra_row.get("review_status") if sra_row else None,
        sra_missing_fields = sra_row.get("missing_fields", "") if sra_row else "",
        biosample_review_status = biosample_row.get("review_status") if biosample_row else None,
        biosample_missing_fields = biosample_row.get("missing_fields", "") if biosample_row else "",
        genbank_review_status = genbank_record.get("review_status") if genbank_record else None,
        genbank_missing = ";".join(f for f in [
            (genbank_record or {}).get("missing_files", ""), (genbank_record or {}).get("missing_metadata", ""),
        ] if f),
    )

def collect_sample_report_data(sample_id: str, runtime_manifest_path: str | Path,
                                statistics_json: str | Path | None = None,
                                biological_validation_json: str | Path | None = None,
                                visualization_result_json: str | Path | None = None,
                                fastp_json: str | Path | None = None,
                                getorganelle_config_json: str | Path | None = None,
                                blast_config_json: str | Path | None = None,
                                mitos2_config_json: str | Path | None = None,
                                sra_metadata_tsv: str | Path | None = None,
                                biosample_metadata_tsv: str | Path | None = None,
                                genbank_record_json: str | Path | None = None,
                                logger: Logger | None = None) -> SampleReportData:
    """
    Assembles one sample's report data from its real result JSON.

    Args:
        sample_id (str): The sample id.
        runtime_manifest_path (str | Path): The runtime manifest TSV.
        statistics_json (str | Path | None): The merged statistics JSON written by run_stats_aggregation.py.
        biological_validation_json (str | Path | None): The validation result JSON written by run_biological_validation.py.
        visualization_result_json (str | Path | None): The result JSON written by run_circular_map.py.
        fastp_json (str | Path | None): fastp's own JSON report for this sample.
        getorganelle_config_json (str | Path | None): The resolved GetOrganelle config JSON.
        blast_config_json (str | Path | None): The resolved BLAST config JSON.
        mitos2_config_json (str | Path | None): The resolved MITOS2 config JSON.
        sra_metadata_tsv (str | Path | None): The SRA metadata table written by run_sra_metadata.py.
        biosample_metadata_tsv (str | Path | None): The BioSample metadata table written by run_biosample_metadata.py.
        genbank_record_json (str | Path | None): This sample's own GenBank record JSON written by run_mitogenome_metadata.py.
        logger (Logger | None): Logger to use. Defaults to None.

    Returns:
        SampleReportData: The assembled report data.
    """
    metadata = read_sample_metadata(runtime_manifest_path, sample_id, logger = logger)

    stats = _read_json(Path(statistics_json) if statistics_json else None) or {}
    assembly = AssemblyStatistics(**stats["assembly"]) if stats.get("assembly") else None
    blast = BlastStatistics(**stats["blast"]) if stats.get("blast") else None
    annotation = None
    if stats.get("annotation"):
        annotation_fields = dict(stats["annotation"])
        annotation_fields["missing_core_genes"] = tuple(annotation_fields.get("missing_core_genes", ()))
        annotation = AnnotationStatistics(**annotation_fields)

    validation = _read_json(Path(biological_validation_json) if biological_validation_json else None)
    overall_verdict, validation_checks = None, ()
    if validation is not None:
        selected_attempt = validation.get("selected_attempt")
        for result in validation.get("results", []):
            if result.get("attempt") == selected_attempt:
                overall_verdict = result.get("overall_verdict")
                validation_checks = tuple(result.get("checks", []))
                break

    visualization = _read_json(Path(visualization_result_json) if visualization_result_json else None)
    visualization_png = None
    if visualization is not None and visualization.get("success") and visualization.get("output_png"):
        visualization_png = Path(visualization["output_png"])

    assembly_config = None
    go_config = _read_json(Path(getorganelle_config_json) if getorganelle_config_json else None)
    if go_config is not None:
        assembly_config = AssemblyConfigSummary(
            organelle_type = go_config.get("organelle_type", "animal_mt"),
            seeded = bool(go_config.get("seeded", False)), seed_file = go_config.get("seed_file"),
        )

    blast_config = None
    b_config = _read_json(Path(blast_config_json) if blast_config_json else None)
    if b_config is not None:
        blast_config = BlastConfigSummary(database = b_config.get("database", ""), mode = b_config.get("mode", ""))

    annotation_config = None
    m_config = _read_json(Path(mitos2_config_json) if mitos2_config_json else None)
    if m_config is not None:
        annotation_config = AnnotationConfigSummary(
            genetic_code = m_config.get("genetic_code", 2), refseqver = m_config.get("refseqver", ""),
        )

    archival = collect_archival_summary(sample_id, sra_metadata_tsv, biosample_metadata_tsv, genbank_record_json)

    return SampleReportData(
        sample_id = sample_id, metadata = metadata, fastp_stats = parse_fastp_stats(fastp_json),
        assembly = assembly, blast = blast, annotation = annotation,
        overall_verdict = overall_verdict, validation_checks = validation_checks,
        visualization_png = visualization_png, assembly_config = assembly_config,
        blast_config = blast_config, annotation_config = annotation_config, archival = archival,
    )

def collect_phylogeny_summary(iqtree_result_json: str | Path | None, model: str | None = None,
                               tree_figure_png: str | Path | None = None) -> PhylogenySummary | None:
    """
    Reads the batch-level phylogeny outcome.

    Args:
        iqtree_result_json (str | Path | None): The result JSON written by run_iqtree.py.
        model (str | None): The substitution model IQ-TREE was run with, when known.
        tree_figure_png (str | Path | None): The rendered tree figure PNG, when run_tree_figure.py produced one.

    Returns:
        PhylogenySummary | None: The outcome, or None if phylogeny didn't run.
    """
    result = _read_json(Path(iqtree_result_json) if iqtree_result_json else None)
    if result is None:
        return None
    figure_path = Path(tree_figure_png) if tree_figure_png else None
    return PhylogenySummary(
        success = result.get("success", False), tip_count = result.get("tip_count"),
        newick_path = Path(result["newick_path"]) if result.get("newick_path") else None,
        reason = result.get("reason", ""), model = model,
        tree_figure_png = figure_path if (figure_path and figure_path.is_file()) else None,
    )

def collect_run_report_data(job_id: str, generated_at: str, samples: list[SampleReportData],
                             phylogeny: PhylogenySummary | None = None) -> RunReportData:
    """
    Assembles the run-level report data from every sample's own
    already-collected report data.

    Args:
        job_id (str): The pipeline job id.
        generated_at (str): An ISO-8601 UTC timestamp.
        samples (list[SampleReportData]): Every sample's report data, in manifest order.
        phylogeny (PhylogenySummary | None): The batch-level phylogeny outcome, if available.

    Returns:
        RunReportData: The assembled run-level report data.
    """
    return RunReportData(job_id = job_id, generated_at = generated_at, samples = tuple(samples), phylogeny = phylogeny)
