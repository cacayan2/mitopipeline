"""
13-report_pdf.smk

Renders each sample's own Markdown + PDF report (#358/#359/#361),
then the run-level Markdown + PDF report (#360/#441/#443) aggregating
every sample plus the batch-level phylogeny outcome.
"""

def _sample_report_inputs(wildcards):
    inputs = {
        "stats_done": str(JOB_DIR / "12-stats_aggregation" / wildcards.sample / f"{wildcards.sample}.stats_aggregation.done"),
        "validation_done": str(JOB_DIR / "06-contig_selection" / wildcards.sample / f"{wildcards.sample}.biological_validation.done"),
        "visualization_done": str(JOB_DIR / "09-visualization" / wildcards.sample / f"{wildcards.sample}.visualization.done"),
    }
    # Archival dependencies are only real dependency edges when those
    # stages are actually enabled -- otherwise their own .done markers
    # are never produced by anything else in the DAG (see run_report's
    # own identical reasoning for phylogeny, below).
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
        statistics_json = str(JOB_DIR / "12-stats_aggregation" / "{sample}" / "{sample}.statistics.json"),
        biological_validation_json = str(JOB_DIR / "06-contig_selection" / "{sample}" / "{sample}.biological_validation.json"),
        visualization_result_json = str(JOB_DIR / "09-visualization" / "{sample}" / "{sample}.circular_map_result.json"),
        fastp_json = str(JOB_DIR / "trimming" / "{sample}" / "{sample}.fastp.json"),
        getorganelle_config_json = str(JOB_DIR / "04a-getorganelle_config" / "{sample}" / "{sample}.getorganelle_config.json"),
        blast_config_json = str(JOB_DIR / "05a-blast_config" / "{sample}" / "{sample}.blast_config.json"),
        mitos2_config_json = str(JOB_DIR / "08a-mitos2_config" / "{sample}" / "{sample}.mitos2_config.json"),
        sra_metadata_tsv = str(JOB_DIR / "14a-sra_archival" / "sra_metadata.tsv"),
        biosample_metadata_tsv = str(JOB_DIR / "14b-biosample_archival" / "biosample_metadata.tsv"),
        genbank_record_json = str(JOB_DIR / "14c-mitogenome_archival" / "{sample}.genbank_record.json"),
        output_dir = str(JOB_DIR / "13-report_pdf" / "{sample}"),
        log_file = str(JOB_DIR / "logs" / "13-report_pdf" / "{sample}.log"),
        global_log_file = str(GLOBAL_LOG_FILE),
        exec_script = str(SRC_DIR / "mitopipeline" / "exec" / "13-report_pdf" / "run_sample_report.py")
    conda:
        "../../envs/13-reporting.yaml"
    shell:
        """
        {PYTHON_BIN:q} {params.exec_script:q} \
            --sample-id {wildcards.sample} \
            --runtime-manifest {params.runtime_manifest:q} \
            --statistics-json {params.statistics_json:q} \
            --biological-validation-json {params.biological_validation_json:q} \
            --visualization-result-json {params.visualization_result_json:q} \
            --fastp-json {params.fastp_json:q} \
            --getorganelle-config-json {params.getorganelle_config_json:q} \
            --blast-config-json {params.blast_config_json:q} \
            --mitos2-config-json {params.mitos2_config_json:q} \
            --sra-metadata-tsv {params.sra_metadata_tsv:q} \
            --biosample-metadata-tsv {params.biosample_metadata_tsv:q} \
            --genbank-record-json {params.genbank_record_json:q} \
            --output-dir {params.output_dir:q} \
            --log-file {params.log_file:q} \
            --global-log-file {params.global_log_file:q} \
            || true
        mkdir -p $(dirname {output.done:q})
        touch {output.done:q}
        """

def _run_report_inputs(wildcards):
    inputs = {
        "sample_reports": expand(
            str(JOB_DIR / "13-report_pdf" / "{sample}" / "{sample}.report_pdf.done"), sample = SAMPLES
        )
    }
    # Phylogeny is only a real dependency edge when the stage is
    # actually enabled -- otherwise phylogeny.done is never produced
    # by anything else in the DAG, and requiring it here would force
    # phylogeny to run even when a user explicitly disabled it.
    if stage_enabled("11-phylogeny"):
        inputs["phylogeny_done"] = str(JOB_DIR / "11-phylogeny" / "phylogeny.done")
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
        iqtree_result_json = str(JOB_DIR / "11-phylogeny" / "iqtree_result.json"),
        phylogeny_model = config["stage_options"]["phylogeny"]["iqtree"].get("model", "") if stage_enabled("11-phylogeny") else "",
        output_dir = str(JOB_DIR / "13-report_pdf"),
        log_file = str(JOB_DIR / "logs" / "13-report_pdf" / "run_report.log"),
        global_log_file = str(GLOBAL_LOG_FILE),
        exec_script = str(SRC_DIR / "mitopipeline" / "exec" / "13-report_pdf" / "run_run_report.py")
    conda:
        "../../envs/13-reporting.yaml"
    shell:
        """
        MODEL_ARGS=""
        if [ -n {params.phylogeny_model:q} ]; then
            MODEL_ARGS="--phylogeny-model {params.phylogeny_model:q}"
        fi
        {PYTHON_BIN:q} {params.exec_script:q} \
            --job-id {params.job_id:q} \
            {params.sample_id_flags} \
            --runtime-manifest {params.runtime_manifest:q} \
            --job-dir {params.job_dir:q} \
            --iqtree-result-json {params.iqtree_result_json:q} \
            $MODEL_ARGS \
            --output-dir {params.output_dir:q} \
            --log-file {params.log_file:q} \
            --global-log-file {params.global_log_file:q} \
            || true
        mkdir -p $(dirname {output.done:q})
        touch {output.done:q}
        """
