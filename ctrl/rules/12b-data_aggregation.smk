"""
12b-data_aggregation.smk

Consolidates every sample's own final assembly/annotation/phylogeny
outputs into one flat aggregated_data/ tree (Emil's own direct ask,
2026-10-28) -- assemblies/, annotations/ (+ annotations/gff/), and
phylogenies/ -- a job-level step covering every sample in one call,
the same pattern run_report.py's own aggregation already uses.

Also consolidates the run's real NCBI archival materials into
archival/ (Emil's own direct follow-up ask), which means this stage
now genuinely depends on the archival stages (#14a-#15) completing
first -- despite its "12b" numeric prefix placing it earlier in the
stage listing than "14a"/"15". That numbering reflects when this
step was first built (assembly/annotation/phylogeny consolidation,
before archival existed), not its real position in the DAG today;
left as-is rather than re-renumbered, since Snakemake schedules off
real declared dependencies, not directory-name order.
"""

def _data_aggregation_inputs(wildcards):
    inputs = {
        "stats_done": expand(
            str(JOB_DIR / "12a-stats_aggregation" / "{sample}" / "{sample}.stats_aggregation.done"), sample = SAMPLES
        ),
    }
    # Visualization/phylogeny/archival are only real dependency edges
    # when those stages are actually enabled -- stats_aggregation's
    # own inputs don't cover any of these, so they're declared
    # explicitly here (the same reasoning 13-report_pdf.smk's own
    # _sample_report_inputs() already uses for the identical gap).
    if stage_enabled("09-visualization"):
        inputs["visualization_done"] = expand(
            str(JOB_DIR / "09-visualization" / "{sample}" / "{sample}.visualization.done"), sample = SAMPLES
        )
    if stage_enabled("11-phylogeny"):
        inputs["tree_figure_done"] = expand(
            str(JOB_DIR / "11-phylogeny" / "{sample}" / "{sample}.tree_figure.done"), sample = SAMPLES
        )
    if stage_enabled("14a-sra_archival"):
        inputs["sra_archival_done"] = str(JOB_DIR / "14a-sra_archival" / "sra_archival.done")
    if stage_enabled("14b-biosample_archival"):
        inputs["biosample_archival_done"] = str(JOB_DIR / "14b-biosample_archival" / "biosample_archival.done")
    if stage_enabled("14c-mitogenome_archival"):
        inputs["mitogenome_archival_done"] = expand(
            str(JOB_DIR / "14c-mitogenome_archival" / "{sample}" / "{sample}.mitogenome_archival.done"), sample = SAMPLES
        )
    # archival_validation is the real source of the submission
    # checklist/readiness report/checksums this step copies in, and
    # itself already depends on 14a/14b/14c -- declaring it alone
    # would transitively cover those three, but each is still
    # declared explicitly above too, since this step also pulls real
    # files directly from 14a/14b/14c (not just from 15's own output).
    if stage_enabled("15-archival_validation"):
        inputs["archival_validation_done"] = str(JOB_DIR / "15-archival_validation" / "archival_validation.done")
    return inputs

rule data_aggregation:
    input:
        unpack(_data_aggregation_inputs)
    output:
        # {job_id}.data_aggregation_result.json is NOT declared -- the
        # exec script always writes it (even when every sample had
        # nothing to aggregate, it records that honestly), so the done
        # marker alone is a sufficient, always-true dependency edge.
        done = str(JOB_DIR / "aggregated_data" / "data_aggregation.done")
    params:
        job_id = JOB_ID,
        sample_id_flags = " ".join(f"--sample-id {s}" for s in SAMPLES),
        job_dir = str(JOB_DIR),
        output_dir = str(JOB_DIR / "aggregated_data"),
        sra_metadata_tsv = str(JOB_DIR / "14a-sra_archival" / "sra_metadata.tsv") if stage_enabled("14a-sra_archival") else "",
        biosample_metadata_tsv = str(JOB_DIR / "14b-biosample_archival" / "biosample_metadata.tsv") if stage_enabled("14b-biosample_archival") else "",
        genbank_metadata_tsv = str(JOB_DIR / "15-archival_validation" / "genbank_metadata.tsv") if stage_enabled("15-archival_validation") else "",
        archival_validation_tsv = str(JOB_DIR / "15-archival_validation" / "archival_validation.tsv") if stage_enabled("15-archival_validation") else "",
        submission_checklist_md = str(JOB_DIR / "15-archival_validation" / "submission_checklist.md") if stage_enabled("15-archival_validation") else "",
        checksums_tsv = str(JOB_DIR / "15-archival_validation" / "md5_checksums.tsv") if stage_enabled("15-archival_validation") else "",
        log_file = str(JOB_DIR / "logs" / "12b-data_aggregation" / "data_aggregation.log"),
        global_log_file = str(GLOBAL_LOG_FILE),
        exec_script = str(SRC_DIR / "mitopipeline" / "exec" / "12b-data_aggregation" / "run_data_aggregation.py")
    conda:
        "../../envs/mitopipeline.yaml"
    shell:
        """
        ARCHIVAL_ARGS=""
        if [ -n "{params.sra_metadata_tsv}" ]; then
            ARCHIVAL_ARGS="$ARCHIVAL_ARGS --sra-metadata-tsv {params.sra_metadata_tsv:q}"
        fi
        if [ -n "{params.biosample_metadata_tsv}" ]; then
            ARCHIVAL_ARGS="$ARCHIVAL_ARGS --biosample-metadata-tsv {params.biosample_metadata_tsv:q}"
        fi
        if [ -n "{params.genbank_metadata_tsv}" ]; then
            ARCHIVAL_ARGS="$ARCHIVAL_ARGS --genbank-metadata-tsv {params.genbank_metadata_tsv:q}"
        fi
        if [ -n "{params.archival_validation_tsv}" ]; then
            ARCHIVAL_ARGS="$ARCHIVAL_ARGS --archival-validation-tsv {params.archival_validation_tsv:q}"
        fi
        if [ -n "{params.submission_checklist_md}" ]; then
            ARCHIVAL_ARGS="$ARCHIVAL_ARGS --submission-checklist-md {params.submission_checklist_md:q}"
        fi
        if [ -n "{params.checksums_tsv}" ]; then
            ARCHIVAL_ARGS="$ARCHIVAL_ARGS --checksums-tsv {params.checksums_tsv:q}"
        fi
        {PYTHON_BIN:q} {params.exec_script:q} \
            --job-id {params.job_id:q} \
            {params.sample_id_flags} \
            --job-dir {params.job_dir:q} \
            --output-dir {params.output_dir:q} \
            $ARCHIVAL_ARGS \
            --log-file {params.log_file:q} \
            --global-log-file {params.global_log_file:q} \
            || true
        mkdir -p $(dirname {output.done:q})
        touch {output.done:q}
        """
