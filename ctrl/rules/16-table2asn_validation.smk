"""
16-table2asn_validation.smk

Locally validates one sample's real GenBank submission folder (the
FASTA + 5-column feature table 14c already prepared) with NCBI's own
table2asn tool -- per-sample, since each one depends only on that
sample's own real GenBank record. Entirely local; never contacts or
submits anything to NCBI.
"""

TABLE2ASN_CONFIG = config["stage_options"]["archival"]["table2asn"]

rule table2asn_validation:
    input:
        mitogenome_archival_done = str(JOB_DIR / "14c-mitogenome_archival" / "{sample}" / "{sample}.mitogenome_archival.done")
    output:
        done = str(JOB_DIR / "16-table2asn_validation" / "{sample}" / "{sample}.table2asn_validation.done")
    params:
        genbank_record_json = str(JOB_DIR / "14c-mitogenome_archival" / "{sample}.genbank_record.json"),
        binary_path = str((PROJECT_ROOT / TABLE2ASN_CONFIG["binary_path"]).resolve()),
        submission_template = str((PROJECT_ROOT / TABLE2ASN_CONFIG["submission_template"]).resolve()),
        topology = TABLE2ASN_CONFIG.get("topology", "circular"),
        genome_location = TABLE2ASN_CONFIG.get("genome_location", "mitochondrion"),
        working_dir = str(JOB_DIR / "16-table2asn_validation"),
        output_dir = str(JOB_DIR / "16-table2asn_validation"),
        metrics_dir = str(JOB_DIR / "metrics" / "16-table2asn_validation"),
        log_file = str(JOB_DIR / "logs" / "16-table2asn_validation" / "{sample}.log"),
        global_log_file = str(GLOBAL_LOG_FILE),
        exec_script = str(SRC_DIR / "mitopipeline" / "exec" / "16-table2asn_validation" / "run_table2asn_validation.py")
    conda:
        "../../envs/mitopipeline.yaml"
    shell:
        """
        {PYTHON_BIN:q} {params.exec_script:q} \
            --sample-id {wildcards.sample} \
            --genbank-record-json {params.genbank_record_json:q} \
            --binary-path {params.binary_path:q} \
            --submission-template {params.submission_template:q} \
            --topology {params.topology:q} \
            --genome-location {params.genome_location:q} \
            --working-dir {params.working_dir:q} \
            --output-dir {params.output_dir:q} \
            --metrics-dir {params.metrics_dir:q} \
            --log-file {params.log_file:q} \
            --global-log-file {params.global_log_file:q} \
            || true
        mkdir -p $(dirname {output.done:q})
        touch {output.done:q}
        """
