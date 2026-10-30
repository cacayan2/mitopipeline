"""
13-report_pdf.smk

Renders each sample's own Markdown + PDF report,
including that sample's own real phylogenetic tree, restored to
v0.1.0's per-sample design), then the run-level Markdown + PDF report
aggregating every sample.
"""

def _recruited_reads_fastqc_dir(selection_path: Path) -> str:
    """
    The selected attempt's own recruited-reads FastQC directory --
    sibling of the same contig_selection.json fasta_path every other
    rule in this file already derives its own selected-attempt path
    from (see 12a-stats_aggregation.smk's identical real reasoning for
    why this can't be a fixed path: the selected attempt may be the
    primary run or any rescue reassembly).
    """
    import json
    try:
        selection = json.loads(Path(selection_path).read_text(encoding = "utf-8"))
    except (OSError, json.JSONDecodeError):
        return ""
    if not selection.get("selected") or not selection.get("fasta_path"):
        return ""
    return str(Path(selection["fasta_path"]).parent / "recruited_reads_fastqc")

def _sample_report_inputs(wildcards):
    inputs = {
        "stats_done": str(JOB_DIR / "12a-stats_aggregation" / wildcards.sample / f"{wildcards.sample}.stats_aggregation.done"),
        "validation_done": str(JOB_DIR / "06-contig_selection" / wildcards.sample / f"{wildcards.sample}.biological_validation.done"),
        "visualization_done": str(JOB_DIR / "09-visualization" / wildcards.sample / f"{wildcards.sample}.visualization.done"),
        "selection": str(JOB_DIR / "06-contig_selection" / wildcards.sample / f"{wildcards.sample}.contig_selection.json"),
        # A real dependency edge (not just a params path) so editing
        # the manifest to add a missing archival field and re-running
        # regenerates this sample's report on its own.
        "runtime_manifest": str(RUNTIME_MANIFEST_PATH),
    }
    if stage_enabled("12b-data_aggregation"):
        inputs["data_aggregation_done"] = str(JOB_DIR / "aggregated_data" / "data_aggregation.done")
    # Raw/trimmed FastQC are only real dependency edges when those
    # stages are actually enabled -- same reasoning as phylogeny/
    # archival below.
    if stage_enabled("01-qc_raw"):
        inputs["qc_raw_done"] = str(JOB_DIR / "qc" / "raw" / wildcards.sample / f"{wildcards.sample}.qc_raw.done")
    if stage_enabled("03-qc_trimmed"):
        inputs["qc_trimmed_done"] = str(JOB_DIR / "qc" / "trimmed" / wildcards.sample / f"{wildcards.sample}.qc_trimmed.done")
    # Phylogeny/archival dependencies are only real dependency edges
    # when those stages are actually enabled -- otherwise their own
    # .done markers are never produced by anything else in the DAG.
    if stage_enabled("11-phylogeny"):
        inputs["tree_figure_done"] = str(JOB_DIR / "11-phylogeny" / wildcards.sample / f"{wildcards.sample}.tree_figure.done")
    if stage_enabled("14a-sra_archival"):
        inputs["sra_archival_done"] = str(JOB_DIR / "14a-sra_archival" / "sra_archival.done")
    if stage_enabled("14b-biosample_archival"):
        inputs["biosample_archival_done"] = str(JOB_DIR / "14b-biosample_archival" / "biosample_archival.done")
    if stage_enabled("14c-mitogenome_archival"):
        inputs["mitogenome_archival_done"] = str(
            JOB_DIR / "14c-mitogenome_archival" / wildcards.sample / f"{wildcards.sample}.mitogenome_archival.done"
        )
    return inputs

rule sample_report:
    input:
        unpack(_sample_report_inputs)
    output:
        # {sample}.report_result.json is NOT declared -- the exec
        # script always writes it (missing stats simply render as
        # "not available" sections), so the done marker alone is a
        # sufficient, always-true dependency edge.
        done = str(JOB_DIR / "13-report_pdf" / "{sample}" / "{sample}.report_pdf.done")
    params:
        runtime_manifest = str(RUNTIME_MANIFEST_PATH),
        statistics_json = str(JOB_DIR / "12a-stats_aggregation" / "{sample}" / "{sample}.statistics.json"),
        biological_validation_json = str(JOB_DIR / "06-contig_selection" / "{sample}" / "{sample}.biological_validation.json"),
        visualization_result_json = str(JOB_DIR / "09-visualization" / "{sample}" / "{sample}.circular_map_result.json"),
        fastp_json = str(JOB_DIR / "trimming" / "{sample}" / "{sample}.fastp.json"),
        getorganelle_config_json = str(JOB_DIR / "04a-getorganelle_config" / "{sample}" / "{sample}.getorganelle_config.json"),
        seed_selection_json = str(JOB_DIR / "04c-seed_sequences" / "{sample}" / "{sample}.seed_selection.json"),
        blast_config_json = str(JOB_DIR / "05a-blast_config" / "{sample}" / "{sample}.blast_config.json"),
        mitos2_config_json = str(JOB_DIR / "08a-mitos2_config" / "{sample}" / "{sample}.mitos2_config.json"),
        sra_metadata_tsv = str(JOB_DIR / "14a-sra_archival" / "sra_metadata.tsv"),
        biosample_metadata_tsv = str(JOB_DIR / "14b-biosample_archival" / "biosample_metadata.tsv"),
        genbank_record_json = str(JOB_DIR / "14c-mitogenome_archival" / "{sample}.genbank_record.json"),
        iqtree_result_json = str(JOB_DIR / "11-phylogeny" / "{sample}" / "iqtree_result.json"),
        phylogeny_model = config["stage_options"]["phylogeny"]["iqtree"].get("model", "") if stage_enabled("11-phylogeny") else "",
        tree_figure_png = str(JOB_DIR / "11-phylogeny" / "{sample}" / "tree_figure.png"),
        raw_r1_fastqc_zip = str(JOB_DIR / "qc" / "raw" / "{sample}" / "{sample}.R1_fastqc.zip"),
        raw_r2_fastqc_zip = str(JOB_DIR / "qc" / "raw" / "{sample}" / "{sample}.R2_fastqc.zip"),
        trimmed_r1_fastqc_zip = str(JOB_DIR / "qc" / "trimmed" / "{sample}" / "{sample}.R1_fastqc.zip"),
        trimmed_r2_fastqc_zip = str(JOB_DIR / "qc" / "trimmed" / "{sample}" / "{sample}.R2_fastqc.zip"),
        recruited_reads_fastqc_dir = lambda wc, input: _recruited_reads_fastqc_dir(input.selection),
        aggregated_data_dir = str(JOB_DIR / "aggregated_data"),
        rescue_dir = str(JOB_DIR / "07-assembly_rescue" / "{sample}"),
        metrics_dir = str(JOB_DIR / "metrics"),
        log_dir = str(JOB_DIR / "logs"),
        pandoc_bin = config["stage_options"]["report_export"]["pandoc"]["bin"],
        pdf_engine = config["stage_options"]["report_export"]["pandoc"]["engine"],
        output_dir = str(JOB_DIR / "13-report_pdf" / "{sample}"),
        log_file = str(JOB_DIR / "logs" / "13-report_pdf" / "{sample}.log"),
        global_log_file = str(GLOBAL_LOG_FILE),
        exec_script = str(SRC_DIR / "mitopipeline" / "exec" / "13-report_pdf" / "run_sample_report.py")
    conda:
        "../../envs/13-reporting.yaml"
    shell:
        """
        MODEL_ARGS=""
        if [ -n "{params.phylogeny_model}" ]; then
            MODEL_ARGS="--phylogeny-model {params.phylogeny_model:q}"
        fi
        RECRUITED_FASTQC_DIR="{params.recruited_reads_fastqc_dir}"
        {PYTHON_BIN:q} {params.exec_script:q} \
            --sample-id {wildcards.sample} \
            --runtime-manifest {params.runtime_manifest:q} \
            --statistics-json {params.statistics_json:q} \
            --biological-validation-json {params.biological_validation_json:q} \
            --visualization-result-json {params.visualization_result_json:q} \
            --fastp-json {params.fastp_json:q} \
            --getorganelle-config-json {params.getorganelle_config_json:q} \
            --seed-selection-json {params.seed_selection_json:q} \
            --blast-config-json {params.blast_config_json:q} \
            --mitos2-config-json {params.mitos2_config_json:q} \
            --sra-metadata-tsv {params.sra_metadata_tsv:q} \
            --biosample-metadata-tsv {params.biosample_metadata_tsv:q} \
            --genbank-record-json {params.genbank_record_json:q} \
            --iqtree-result-json {params.iqtree_result_json:q} \
            $MODEL_ARGS \
            --tree-figure-png {params.tree_figure_png:q} \
            --raw-r1-fastqc-zip {params.raw_r1_fastqc_zip:q} \
            --raw-r2-fastqc-zip {params.raw_r2_fastqc_zip:q} \
            --trimmed-r1-fastqc-zip {params.trimmed_r1_fastqc_zip:q} \
            --trimmed-r2-fastqc-zip {params.trimmed_r2_fastqc_zip:q} \
            --recruited-r1-fastqc-zip "$RECRUITED_FASTQC_DIR/recruited_R1_fastqc.zip" \
            --recruited-r2-fastqc-zip "$RECRUITED_FASTQC_DIR/recruited_R2_fastqc.zip" \
            --aggregated-data-dir {params.aggregated_data_dir:q} \
            --rescue-dir {params.rescue_dir:q} \
            --metrics-dir {params.metrics_dir:q} \
            --log-dir {params.log_dir:q} \
            --pandoc-bin {params.pandoc_bin:q} \
            --pdf-engine {params.pdf_engine:q} \
            --output-dir {params.output_dir:q} \
            --log-file {params.log_file:q} \
            --global-log-file {params.global_log_file:q} \
            || true
        mkdir -p $(dirname {output.done:q})
        touch {output.done:q}
        """

def _run_report_inputs(wildcards):
    inputs = {
        # Every sample's own report already depends on its own real
        # phylogeny outcome (above) -- phylogeny is a transitive
        # dependency here now, not a separate batch-level edge.
        "sample_reports": expand(
            str(JOB_DIR / "13-report_pdf" / "{sample}" / "{sample}.report_pdf.done"), sample = SAMPLES
        )
    }
    if stage_enabled("14a-sra_archival"):
        inputs["sra_archival_done"] = str(JOB_DIR / "14a-sra_archival" / "sra_archival.done")
    if stage_enabled("14b-biosample_archival"):
        inputs["biosample_archival_done"] = str(JOB_DIR / "14b-biosample_archival" / "biosample_archival.done")
    if stage_enabled("14c-mitogenome_archival"):
        inputs["mitogenome_archival_done"] = sample_target("14c-mitogenome_archival", "{sample}.mitogenome_archival.done")
    return inputs

rule run_report:
    input:
        unpack(_run_report_inputs)
    output:
        done = str(JOB_DIR / "13-report_pdf" / "run_report.done")
    params:
        job_id = JOB_ID,
        sample_id_flags = " ".join(f"--sample-id {s}" for s in SAMPLES),
        runtime_manifest = str(RUNTIME_MANIFEST_PATH),
        job_dir = str(JOB_DIR),
        phylogeny_model = config["stage_options"]["phylogeny"]["iqtree"].get("model", "") if stage_enabled("11-phylogeny") else "",
        pandoc_bin = config["stage_options"]["report_export"]["pandoc"]["bin"],
        pdf_engine = config["stage_options"]["report_export"]["pandoc"]["engine"],
        output_dir = str(JOB_DIR / "13-report_pdf"),
        log_file = str(JOB_DIR / "logs" / "13-report_pdf" / "run_report.log"),
        global_log_file = str(GLOBAL_LOG_FILE),
        exec_script = str(SRC_DIR / "mitopipeline" / "exec" / "13-report_pdf" / "run_run_report.py")
    conda:
        "../../envs/13-reporting.yaml"
    shell:
        """
        MODEL_ARGS=""
        if [ -n "{params.phylogeny_model}" ]; then
            MODEL_ARGS="--phylogeny-model {params.phylogeny_model:q}"
        fi
        {PYTHON_BIN:q} {params.exec_script:q} \
            --job-id {params.job_id:q} \
            {params.sample_id_flags} \
            --runtime-manifest {params.runtime_manifest:q} \
            --job-dir {params.job_dir:q} \
            $MODEL_ARGS \
            --pandoc-bin {params.pandoc_bin:q} \
            --pdf-engine {params.pdf_engine:q} \
            --output-dir {params.output_dir:q} \
            --log-file {params.log_file:q} \
            --global-log-file {params.global_log_file:q} \
            || true
        mkdir -p $(dirname {output.done:q})
        touch {output.done:q}
        """
