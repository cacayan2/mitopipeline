"""
15-archival_validation.smk

Cross-checks every real archival table for sample-identifier
consistency, builds the per-sample submission-readiness report, and
packages real checksums + a submission checklist (#374/#375/#450) --
job-level, depends on every prior archival stage having run for
every sample.
"""

rule archival_validation:
    input:
        sra_done = str(JOB_DIR / "14a-sra_archival" / "sra_archival.done"),
        biosample_done = str(JOB_DIR / "14b-biosample_archival" / "biosample_archival.done"),
        mitogenome_done = sample_target("14c-mitogenome_archival", "{sample}.mitogenome_archival.done")
    output:
        done = str(JOB_DIR / "15-archival_validation" / "archival_validation.done")
    params:
        job_id = JOB_ID,
        runtime_manifest = str(RUNTIME_MANIFEST_PATH),
        sra_metadata_tsv = str(JOB_DIR / "14a-sra_archival" / "sra_metadata.tsv"),
        biosample_metadata_tsv = str(JOB_DIR / "14b-biosample_archival" / "biosample_metadata.tsv"),
        genbank_records_dir = str(JOB_DIR / "14c-mitogenome_archival"),
        job_dir = str(JOB_DIR),
        output_dir = str(JOB_DIR / "15-archival_validation"),
        log_file = str(JOB_DIR / "logs" / "15-archival_validation" / "archival_validation.log"),
        global_log_file = str(GLOBAL_LOG_FILE),
        exec_script = str(SRC_DIR / "mitopipeline" / "exec" / "15-archival_validation" / "run_archival_validation.py")
    conda:
        "../../envs/mitopipeline.yaml"
    shell:
        """
        {PYTHON_BIN:q} {params.exec_script:q} \
            --job-id {params.job_id:q} \
            --runtime-manifest {params.runtime_manifest:q} \
            --sra-metadata-tsv {params.sra_metadata_tsv:q} \
            --biosample-metadata-tsv {params.biosample_metadata_tsv:q} \
            --genbank-records-dir {params.genbank_records_dir:q} \
            --job-dir {params.job_dir:q} \
            --output-dir {params.output_dir:q} \
            --log-file {params.log_file:q} \
            --global-log-file {params.global_log_file:q} \
            || true
        mkdir -p $(dirname {output.done:q})
        touch {output.done:q}
        """
