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
from datetime import datetime
from logging import Logger
from pathlib import Path
from mitopipeline.utils.fastqc.parse import extract_fastqc_image, parse_fastqc_zip
from mitopipeline.models.pipeline_job.pipeline_job import aggregate_job_history
from mitopipeline.utils.reporting.compute_metrics import collect_compute_summary, collect_machine_info
from mitopipeline.utils.reporting.models import (
    AnnotationConfigSummary,
    ArchivalSummary,
    AssemblyConfigSummary,
    BlastConfigSummary,
    JobMetadataSummary,
    PhylogenySummary,
    ReadQCSummary,
    RescueAttemptSummary,
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

def _parse_fastp_settings(command: str) -> dict:
    """
    Extracts the real trimming thresholds fastp was actually invoked
    with, from its own command string (fastp's JSON report records
    its own exact real command line, confirmed directly) -- never a
    second, independently-maintained copy of the configured values
    that could drift from what actually ran.

    Args:
        command (str): fastp's own real command string.

    Returns:
        dict: The settings found (missing/unparseable flags are simply absent, never guessed).
    """
    tokens = command.split()
    flag_map = {
        "--qualified_quality_phred": "qualified_quality_phred",
        "--length_required": "length_required",
        "--cut_window_size": "cut_window_size",
        "--cut_mean_quality": "cut_mean_quality",
        "--n_base_limit": "n_base_limit",
        "--unqualified_percent_limit": "unqualified_percent_limit",
    }
    settings: dict = {}
    for index, token in enumerate(tokens):
        key = flag_map.get(token)
        if key is not None and index + 1 < len(tokens):
            try:
                settings[key] = int(tokens[index + 1])
            except ValueError:
                pass
    settings["adapter_detection"] = "--detect_adapter_for_pe" in tokens
    return settings

def parse_fastp_stats(fastp_json_path: str | Path | None) -> dict | None:
    """
    Parses fastp's own real JSON report into the summary fields the
    report displays, including the real trimming settings it was
    actually invoked with.

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
            **_parse_fastp_settings(str(data.get("command", ""))),
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

    # missing_files and missing_metadata are each already de-duplicated
    # within themselves (genbank_metadata.py), but the two lists can
    # legitimately share a field NAME for two different real reasons
    # (e.g. "assembly_fasta" missing as a FILE and "assembly_fasta"
    # missing as a required METADATA field) -- de-duplicating the
    # UNION here is what keeps the report's own displayed list from
    # showing the same name twice, without touching either source
    # list's own real, independently-correct meaning.
    genbank_missing_parts: dict[str, None] = {}
    for raw in ((genbank_record or {}).get("missing_files", ""), (genbank_record or {}).get("missing_metadata", "")):
        for field in raw.split(";"):
            field = field.strip()
            if field:
                genbank_missing_parts.setdefault(field, None)

    return ArchivalSummary(
        sra_review_status = sra_row.get("review_status") if sra_row else None,
        sra_missing_fields = sra_row.get("missing_fields", "") if sra_row else "",
        biosample_review_status = biosample_row.get("review_status") if biosample_row else None,
        biosample_missing_fields = biosample_row.get("missing_fields", "") if biosample_row else "",
        genbank_review_status = genbank_record.get("review_status") if genbank_record else None,
        genbank_missing = ";".join(genbank_missing_parts),
    )

def collect_sample_report_data(sample_id: str, runtime_manifest_path: str | Path,
                                statistics_json: str | Path | None = None,
                                biological_validation_json: str | Path | None = None,
                                visualization_result_json: str | Path | None = None,
                                fastp_json: str | Path | None = None,
                                getorganelle_config_json: str | Path | None = None,
                                seed_selection_json: str | Path | None = None,
                                blast_config_json: str | Path | None = None,
                                mitos2_config_json: str | Path | None = None,
                                sra_metadata_tsv: str | Path | None = None,
                                biosample_metadata_tsv: str | Path | None = None,
                                genbank_record_json: str | Path | None = None,
                                iqtree_result_json: str | Path | None = None,
                                phylogeny_model: str | None = None,
                                tree_figure_png: str | Path | None = None,
                                raw_r1_fastqc_zip: str | Path | None = None,
                                raw_r2_fastqc_zip: str | Path | None = None,
                                trimmed_r1_fastqc_zip: str | Path | None = None,
                                trimmed_r2_fastqc_zip: str | Path | None = None,
                                recruited_r1_fastqc_zip: str | Path | None = None,
                                recruited_r2_fastqc_zip: str | Path | None = None,
                                read_qc_image_dir: str | Path | None = None,
                                aggregated_data_dir: str | Path | None = None,
                                rescue_dir: str | Path | None = None,
                                metrics_dir: str | Path | None = None,
                                log_dir: str | Path | None = None,
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
        seed_selection_json (str | Path | None): This sample's own real seed-selection result JSON, written by run_seed_selection.py.
        blast_config_json (str | Path | None): The resolved BLAST config JSON.
        mitos2_config_json (str | Path | None): The resolved MITOS2 config JSON.
        sra_metadata_tsv (str | Path | None): The SRA metadata table written by run_sra_metadata.py.
        biosample_metadata_tsv (str | Path | None): The BioSample metadata table written by run_biosample_metadata.py.
        genbank_record_json (str | Path | None): This sample's own GenBank record JSON written by run_mitogenome_metadata.py.
        iqtree_result_json (str | Path | None): This sample's own result JSON written by run_iqtree.py.
        phylogeny_model (str | None): The substitution model IQ-TREE was run with, when known.
        tree_figure_png (str | Path | None): This sample's own rendered tree figure PNG, if run_tree_figure.py produced one.
        raw_r1_fastqc_zip (str | Path | None): This sample's own raw R1 FastQC zip, written by run_fastqc_raw.py.
        raw_r2_fastqc_zip (str | Path | None): This sample's own raw R2 FastQC zip, written by run_fastqc_raw.py.
        trimmed_r1_fastqc_zip (str | Path | None): This sample's own trimmed R1 FastQC zip, written by run_fastqc_trimmed.py.
        trimmed_r2_fastqc_zip (str | Path | None): This sample's own trimmed R2 FastQC zip, written by run_fastqc_trimmed.py.
        recruited_r1_fastqc_zip (str | Path | None): The selected attempt's own recruited-reads R1 FastQC zip, written by GetOrganelleRunner.
        recruited_r2_fastqc_zip (str | Path | None): The selected attempt's own recruited-reads R2 FastQC zip, written by GetOrganelleRunner.
        read_qc_image_dir (str | Path | None): Where to extract each read-set's own representative FastQC chart. None skips image extraction.
        aggregated_data_dir (str | Path | None): The run's own consolidated aggregated_data/ directory, written by run_data_aggregation.py, or None if unavailable.
        rescue_dir (str | Path | None): This sample's own real 07-assembly_rescue/{sample} directory, or None if unavailable.
        metrics_dir (str | Path | None): The job's own metrics/ directory, filtered to this sample's own real compute-metrics records, or None if unavailable.
        log_dir (str | Path | None): The job's own logs/ directory, for linking to this sample's own real per-stage log files, or None if unavailable.
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
        # seed_selection.json carries the real seed gene/accession/
        # source -- getorganelle_config.json only ever records the
        # seed FILE path, never which real reference it came from.
        seed_selection = _read_json(Path(seed_selection_json) if seed_selection_json else None)
        seeded = bool(go_config.get("seeded", False))
        assembly_config = AssemblyConfigSummary(
            organelle_type = go_config.get("organelle_type", "animal_mt"),
            seeded = seeded, seed_file = go_config.get("seed_file"),
            seed_gene = seed_selection.get("gene") if seed_selection else None,
            seed_accession = seed_selection.get("accession") if (seed_selection and seeded) else None,
            seed_source = seed_selection.get("source") if seed_selection else None,
            seed_fallback_reason = seed_selection.get("fallback_reason") if (seed_selection and not seeded) else None,
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
    phylogeny = collect_phylogeny_summary(iqtree_result_json, model = phylogeny_model, tree_figure_png = tree_figure_png)

    image_dir = Path(read_qc_image_dir) if read_qc_image_dir else None
    raw_read_qc = collect_read_qc(
        raw_r1_fastqc_zip, raw_r2_fastqc_zip,
        image_output_path = (image_dir / "raw_per_base_quality.png") if image_dir else None, logger = logger,
    )
    trimmed_read_qc = collect_read_qc(
        trimmed_r1_fastqc_zip, trimmed_r2_fastqc_zip,
        image_output_path = (image_dir / "trimmed_per_base_quality.png") if image_dir else None, logger = logger,
    )
    recruited_read_qc = collect_read_qc(
        recruited_r1_fastqc_zip, recruited_r2_fastqc_zip,
        image_output_path = (image_dir / "recruited_per_base_quality.png") if image_dir else None, logger = logger,
    )

    return SampleReportData(
        sample_id = sample_id, metadata = metadata, fastp_stats = parse_fastp_stats(fastp_json),
        assembly = assembly, blast = blast, annotation = annotation,
        overall_verdict = overall_verdict, validation_checks = validation_checks,
        visualization_png = visualization_png, assembly_config = assembly_config,
        blast_config = blast_config, annotation_config = annotation_config, archival = archival,
        phylogeny = phylogeny,
        raw_read_qc = raw_read_qc, trimmed_read_qc = trimmed_read_qc, recruited_read_qc = recruited_read_qc,
        aggregated_data_dir = Path(aggregated_data_dir) if aggregated_data_dir else None,
        rescue_attempts = collect_rescue_summary(rescue_dir),
        compute_summary = collect_compute_summary(Path(metrics_dir), sample_id = sample_id) if metrics_dir else (),
        log_dir = Path(log_dir) if log_dir else None,
        metrics_dir = Path(metrics_dir) if metrics_dir else None,
        machine_info = collect_machine_info(),
    )

def collect_phylogeny_summary(iqtree_result_json: str | Path | None, model: str | None = None,
                               tree_figure_png: str | Path | None = None) -> PhylogenySummary | None:
    """
    Reads one sample's own real phylogeny outcome.

    Args:
        iqtree_result_json (str | Path | None): This sample's own result JSON written by run_iqtree.py.
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

def collect_read_qc(r1_zip: str | Path | None, r2_zip: str | Path | None,
                     image_output_path: str | Path | None = None,
                     logger: Logger | None = None) -> ReadQCSummary | None:
    """
    Combines one read-set's R1/R2 FastQC zips into one display-ready
    summary, extracting R1's own real per-base-quality chart as a
    standalone PNG for the report to embed as a figure.

    Args:
        r1_zip (str | Path | None): This read-set's own R1 FastQC zip.
        r2_zip (str | Path | None): This read-set's own R2 FastQC zip.
        image_output_path (str | Path | None): Where to write R1's extracted per-base-quality chart. None skips image extraction.
        logger (Logger | None): Logger to use. Defaults to None.

    Returns:
        ReadQCSummary | None: The combined summary, or None when neither zip could be read.
    """
    r1_summary = parse_fastqc_zip(r1_zip, logger = logger) if r1_zip else None
    r2_summary = parse_fastqc_zip(r2_zip, logger = logger) if r2_zip else None
    if r1_summary is None and r2_summary is None:
        return None

    counts = [s.total_sequences for s in (r1_summary, r2_summary) if s and s.total_sequences is not None]
    total_sequences = sum(counts) if counts else None

    gc_values = [s.percent_gc for s in (r1_summary, r2_summary) if s and s.percent_gc is not None]
    percent_gc = (sum(gc_values) / len(gc_values)) if gc_values else None

    sequence_length = next((s.sequence_length for s in (r1_summary, r2_summary) if s and s.sequence_length), None)

    verdict_rank = {"PASS": 0, "WARN": 1, "FAIL": 2}
    verdicts = [s.overall_verdict for s in (r1_summary, r2_summary) if s is not None]
    overall_verdict = max(verdicts, key = lambda v: verdict_rank.get(v, 0)) if verdicts else "PASS"

    # Combines R1's and R2's own per-module calls into the single
    # worst verdict per module name, so a module that FAILs on R1 but
    # only WARNs on R2 is correctly reported as a FAIL reason, not
    # silently dropped by whichever read happened second.
    worst_module_verdict: dict[str, str] = {}
    for summary in (r1_summary, r2_summary):
        if summary is None:
            continue
        for module, verdict in summary.module_verdicts.items():
            if verdict_rank.get(verdict, 0) > verdict_rank.get(worst_module_verdict.get(module, "PASS"), 0):
                worst_module_verdict[module] = verdict
    failed_modules = tuple(sorted(m for m, v in worst_module_verdict.items() if v == "FAIL"))
    warned_modules = tuple(sorted(m for m, v in worst_module_verdict.items() if v == "WARN"))

    image_path = None
    if r1_zip is not None and image_output_path is not None:
        image_path = extract_fastqc_image(r1_zip, "per_base_quality.png", image_output_path, logger = logger)

    return ReadQCSummary(
        total_sequences = total_sequences, percent_gc = percent_gc, sequence_length = sequence_length,
        overall_verdict = overall_verdict, failed_modules = failed_modules, warned_modules = warned_modules,
        representative_image_png = image_path,
    )

def _attempt_number(attempt_dir_name: str) -> int:
    """Extracts the integer attempt number from a real "attempt_N" directory name."""
    try:
        return int(attempt_dir_name.rsplit("_", 1)[-1])
    except ValueError:
        return 0

def collect_rescue_summary(rescue_dir: str | Path | None) -> tuple[RescueAttemptSummary, ...]:
    """
    Reads every real rescue attempt made for one sample, in attempt
    order -- whether each attempt actually ran (its own real
    rescue_mapping.json) and, if it did, what GetOrganelle actually
    produced (the reassembly's own real assembly_status.json).
    Covers every attempt the pipeline made, not only one that was
    ultimately selected, so the methods section can report rescue
    honestly even when every attempt failed or was skipped.

    Args:
        rescue_dir (str | Path | None): This sample's own real 07-assembly_rescue/{sample} directory.

    Returns:
        tuple[RescueAttemptSummary, ...]: Every real attempt found, in attempt order. Empty when rescue_dir is unavailable or has no attempts recorded.
    """
    if rescue_dir is None:
        return ()
    rescue_dir = Path(rescue_dir)
    if not rescue_dir.is_dir():
        return ()

    attempts: list[RescueAttemptSummary] = []
    for attempt_dir in sorted(rescue_dir.glob("attempt_*"), key = lambda p: _attempt_number(p.name)):
        mapping_json_path = next(attempt_dir.glob("*.rescue_mapping.json"), None)
        mapping = _read_json(mapping_json_path) if mapping_json_path else None
        if mapping is None:
            continue

        reassembly_status = None
        sample_id = mapping.get("sample_id")
        if sample_id:
            reassembly_status = _read_json(attempt_dir / "reassembly" / f"{sample_id}.assembly_status.json")

        attempts.append(RescueAttemptSummary(
            attempt = _attempt_number(attempt_dir.name),
            attempted = bool(mapping.get("attempted", False)),
            reason = mapping.get("reason", ""),
            profile_name = reassembly_status.get("profile_name") if reassembly_status else None,
            status = reassembly_status.get("status") if reassembly_status else None,
        ))
    return tuple(attempts)

def collect_job_metadata(current_metadata_json: str | Path | None, historical_metadata_json: str | Path | None,
                          job_id: str) -> "JobMetadataSummary | None":
    """
    Reads this job's own real lifecycle metadata -- the current
    attempt's own status/runtime from current_metadata.json, plus a
    cumulative summary across every real attempt this job_id has ever
    had from historical_metadata.json. Never derived from a live
    PipelineJob instance, which has real side effects (creates
    directories, writes a new "created" entry) a read-only report must
    never trigger.

    The report-generation rule itself runs as part of the same
    Snakemake DAG it's reporting on, so current_metadata.json's own
    terminal "completed"/"failed" entry -- written only after the
    whole run (including this very report) finishes -- can never be
    visible yet; its own "status" will structurally always still read
    "created" with the original runtime_seconds=0.0 placeholder. Real
    elapsed time (now minus the recorded start_time) is substituted
    in that case, so the report shows an honest "how long has this
    taken so far" instead of a stale, confusing zero.

    Args:
        current_metadata_json (str | Path | None): The job's own current_metadata.json, or None if unavailable.
        historical_metadata_json (str | Path | None): The job's own historical_metadata.json, or None if unavailable.
        job_id (str): The job id.

    Returns:
        JobMetadataSummary | None: The assembled summary, or None when neither real source is available.
    """
    current = _read_json(Path(current_metadata_json)) if current_metadata_json else None
    history = (
        aggregate_job_history(Path(historical_metadata_json), job_id)
        if historical_metadata_json and Path(historical_metadata_json).is_file() else None
    )
    if current is None and history is None:
        return None
    history = history or {}

    runtime_seconds = current.get("runtime_seconds") if current else None
    if current is not None and current.get("status") == "created" and current.get("start_time"):
        try:
            elapsed = (datetime.now() - datetime.fromisoformat(current["start_time"])).total_seconds()
            runtime_seconds = max(0.0, elapsed)
        except ValueError:
            pass

    return JobMetadataSummary(
        status = current.get("status") if current else None,
        start_time = current.get("start_time") if current else None,
        runtime_seconds = runtime_seconds,
        sample_count = current.get("sample_count") if current else None,
        attempt_count = history.get("attempt_count", 0),
        cumulative_runtime_seconds = history.get("cumulative_runtime_seconds", 0.0),
        incomplete_attempts = history.get("incomplete_attempts", 0),
    )

def collect_run_report_data(job_id: str, generated_at: str, samples: list[SampleReportData],
                             job_dir: str | Path | None = None,
                             current_metadata_json: str | Path | None = None,
                             historical_metadata_json: str | Path | None = None,
                             global_log_file: str | Path | None = None) -> RunReportData:
    """
    Assembles the run-level report data from every sample's own
    already-collected report data (each one carrying its own real
    phylogeny outcome).

    Args:
        job_id (str): The pipeline job id.
        generated_at (str): An ISO-8601 UTC timestamp.
        samples (list[SampleReportData]): Every sample's report data, in manifest order.
        job_dir (str | Path | None): The job's own output root, used to locate metrics/ for a batch-wide compute summary and to link to logs/metrics directly.
        current_metadata_json (str | Path | None): The job's own current_metadata.json, or None if unavailable.
        historical_metadata_json (str | Path | None): The job's own historical_metadata.json, or None if unavailable.
        global_log_file (str | Path | None): The job's own shared global log file, or None if unavailable.

    Returns:
        RunReportData: The assembled run-level report data.
    """
    return RunReportData(
        job_id = job_id, generated_at = generated_at, samples = tuple(samples),
        job_metadata = collect_job_metadata(current_metadata_json, historical_metadata_json, job_id),
        compute_summary = collect_compute_summary(Path(job_dir) / "metrics") if job_dir else (),
        job_dir = Path(job_dir) if job_dir else None,
        global_log_file = Path(global_log_file) if global_log_file else None,
        machine_info = collect_machine_info(),
    )
