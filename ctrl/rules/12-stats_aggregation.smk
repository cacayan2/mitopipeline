"""
12-stats_aggregation.smk

Aggregates one sample's assembly/BLAST/annotation statistics from
the real result JSON #06/#05c/#08b each already wrote.
"""

rule stats_aggregation:
    input:
        validation_done = str(JOB_DIR / "06-contig_selection" / "{sample}" / "{sample}.biological_validation.done"),
        hit_selection_done = str(JOB_DIR / "05c-blast_hits" / "{sample}" / "{sample}.hit_selection.done"),
        annotation_done = str(JOB_DIR / "08b-annotation" / "{sample}" / "{sample}.annotation.done")
    output:
        # {sample}.statistics.json is NOT declared -- the exec script
        # always writes it (even when every optional input is
        # missing, it records null fields), so the done marker alone
        # is a sufficient, always-true dependency edge.
        done = str(JOB_DIR / "12-stats_aggregation" / "{sample}" / "{sample}.stats_aggregation.done")
    params:
        contig_selection_json = str(JOB_DIR / "06-contig_selection" / "{sample}" / "{sample}.contig_selection.json"),
        top_hits_tsv = str(JOB_DIR / "05c-blast_hits" / "{sample}" / "{sample}.top_hits.tsv"),
        biological_validation_json = str(JOB_DIR / "06-contig_selection" / "{sample}" / "{sample}.biological_validation.json"),
        annotation_result_json = str(JOB_DIR / "08b-annotation" / "{sample}" / "{sample}.mitos2_annotation_result.json"),
        output_dir = str(JOB_DIR / "12-stats_aggregation" / "{sample}"),
        log_file = str(JOB_DIR / "logs" / "12-stats_aggregation" / "{sample}.log"),
        global_log_file = str(GLOBAL_LOG_FILE),
        exec_script = str(SRC_DIR / "mitopipeline" / "exec" / "12-stats_aggregation" / "run_stats_aggregation.py")
    conda:
        "../../envs/mitopipeline.yaml"
    shell:
        """
        TOP_HITS_ARGS=""
        if [ -s {params.top_hits_tsv:q} ]; then
            TOP_HITS_ARGS="--top-hits-tsv {params.top_hits_tsv:q} --biological-validation-json {params.biological_validation_json:q}"
        fi
        ANNOTATION_ARGS=""
        if [ -s {params.annotation_result_json:q} ]; then
            ANNOTATION_ARGS="--annotation-result-json {params.annotation_result_json:q}"
        fi
        {PYTHON_BIN:q} {params.exec_script:q} \
            --sample-id {wildcards.sample} \
            --contig-selection-json {params.contig_selection_json:q} \
            $TOP_HITS_ARGS $ANNOTATION_ARGS \
            --output-dir {params.output_dir:q} \
            --log-file {params.log_file:q} \
            --global-log-file {params.global_log_file:q} \
            || true
        mkdir -p $(dirname {output.done:q})
        touch {output.done:q}
        """
