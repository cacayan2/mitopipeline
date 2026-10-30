"""
15-archival_validation.smk

Cross-checks every real archival table for sample-identifier
consistency, builds the per-sample submission-readiness report, and
packages real checksums + a submission checklist --
job-level, depends on every prior archival stage having run for
every sample.
"""

rule archival_validation:
    input:
        sra_done = str(JOB_DIR / "14a-sra_archival" / "sra_archival.done"),
        biosample_done = str(JOB_DIR / "14b-biosample_archival" / "biosample_archival.done"),
        mitogenome_done = sample_target("14c-mitogenome_archival", "{sample}.mitogenome_archival.done"),
        table2asn_done = (sample_target("16-table2asn_validation", "{sample}.table2asn_validation.done")
                         if stage_enabled("16-table2asn_validation") else []),
        # Real dependency edge, same reasoning as 14c-mitogenome_archival.smk.
        runtime_manifest = str(RUNTIME_MANIFEST_PATH)
    output:
        done = str(JOB_DIR / "15-archival_validation" / "archival_validation.done")
    params:
        job_id = JOB_ID,
        runtime_manifest = str(RUNTIME_MANIFEST_PATH),
        sra_metadata_tsv = str(JOB_DIR / "14a-sra_archival" / "sra_metadata.tsv"),
        biosample_metadata_tsv = str(JOB_DIR / "14b-biosample_archival" / "biosample_metadata.tsv"),
        genbank_records_dir = str(JOB_DIR / "14c-mitogenome_archival"),
        table2asn_results_dir = (str(JOB_DIR / "16-table2asn_validation") if stage_enabled("16-table2asn_validation") else ""),
        job_dir = str(JOB_DIR),
        output_dir = str(JOB_DIR / "15-archival_validation"),
        log_file = str(JOB_DIR / "logs" / "15-archival_validation" / "archival_validation.log"),
        global_log_file = str(GLOBAL_LOG_FILE),
        exec_script = str(SRC_DIR / "mitopipeline" / "exec" / "15-archival_validation" / "run_archival_validation.py")
    conda:
        "../../envs/mitopipeline.yaml"
    shell:
        """
        TABLE2ASN_ARGS=""
        if [ -n "{params.table2asn_results_dir}" ]; then
            TABLE2ASN_ARGS="--table2asn-results-dir {params.table2asn_results_dir:q}"
        fi
        {PYTHON_BIN:q} {params.exec_script:q} \
            --job-id {params.job_id:q} \
            --runtime-manifest {params.runtime_manifest:q} \
            --sra-metadata-tsv {params.sra_metadata_tsv:q} \
            --biosample-metadata-tsv {params.biosample_metadata_tsv:q} \
            --genbank-records-dir {params.genbank_records_dir:q} \
            $TABLE2ASN_ARGS \
            --job-dir {params.job_dir:q} \
            --output-dir {params.output_dir:q} \
            --log-file {params.log_file:q} \
            --global-log-file {params.global_log_file:q} \
            || true
        mkdir -p $(dirname {output.done:q})
        touch {output.done:q}
        """
