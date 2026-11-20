"""
12a-stats_aggregation.smk

Aggregates one sample's assembly/BLAST/annotation statistics from
the real result JSON #06/#05c/#08b each already wrote.
"""

BLAST_CONFIG = config["stage_options"]["phylogeny"]["blast"]

def _recruited_reads_qc_json(selection_path: Path) -> str:
    """
    Reads a sample's own selected assembly gfa_path from its real
    contig_selection.json and returns the sibling recruited-reads QC
    record's own path -- GetOrganelleRunner always writes
    {sample}.recruited_reads_qc.json into the exact same output_dir
    as the selected attempt's own data.fasta/graph.gfa, whether that
    attempt was the primary run or a rescue reassembly, so deriving
    it this way (rather than guessing a fixed path) is correct for
    both cases. Empty string when nothing was selected, mirroring
    the project's own _selected_fasta() helpers elsewhere.

    Confirmed directly, 2026-11-19, against real run output: must use
    gfa_path, not fasta_path -- contig_selection.json's own fasta_path
    points to contig_selection's own copy (06-contig_selection/
    {sample}/{sample}.selected.fasta), never the real GetOrganelle
    attempt directory; gfa_path is the one field that genuinely stays
    inside that attempt's own output_dir for both a primary run
    (04d-assembly/{sample}/primary/) and a rescue reassembly
    (07-assembly_rescue/{sample}/attempt_N/reassembly/). Using
    fasta_path here silently produced a path that never existed --
    recruited-reads QC was missing for every sample, not a few.
    """
    import json
    try:
        selection = json.loads(Path(selection_path).read_text(encoding = "utf-8"))
    except (OSError, json.JSONDecodeError):
        return ""
    if not selection.get("selected") or not selection.get("gfa_path"):
        return ""
    sample_id = selection.get("sample_id", "")
    return str(Path(selection["gfa_path"]).parent / f"{sample_id}.recruited_reads_qc.json")

rule stats_aggregation:
    input:
        selection = str(JOB_DIR / "06-contig_selection" / "{sample}" / "{sample}.contig_selection.json"),
        validation_done = str(JOB_DIR / "06-contig_selection" / "{sample}" / "{sample}.biological_validation.done"),
        hit_selection_done = str(JOB_DIR / "05c-blast_hits" / "{sample}" / "{sample}.hit_selection.done"),
        circular_coverage_done = str(JOB_DIR / "05c-blast_hits" / "{sample}" / "{sample}.circular_coverage.done"),
        annotation_done = str(JOB_DIR / "08b-annotation" / "{sample}" / "{sample}.annotation.done")
    output:
        # {sample}.statistics.json is NOT declared -- the exec script
        # always writes it (even when every optional input is
        # missing, it records null fields), so the done marker alone
        # is a sufficient, always-true dependency edge.
        done = str(JOB_DIR / "12a-stats_aggregation" / "{sample}" / "{sample}.stats_aggregation.done")
    params:
        contig_selection_json = str(JOB_DIR / "06-contig_selection" / "{sample}" / "{sample}.contig_selection.json"),
        recruited_reads_qc_json = lambda wc, input: _recruited_reads_qc_json(input.selection),
        top_hits_tsv = str(JOB_DIR / "05c-blast_hits" / "{sample}" / "{sample}.top_hits.tsv"),
        biological_validation_json = str(JOB_DIR / "06-contig_selection" / "{sample}" / "{sample}.biological_validation.json"),
        circular_coverage_json = str(JOB_DIR / "05c-blast_hits" / "{sample}" / "{sample}.circular_coverage_result.json"),
        annotation_result_json = str(JOB_DIR / "08b-annotation" / "{sample}" / "{sample}.mitos2_annotation_result.json"),
        reference_fasta = str((PROJECT_ROOT / BLAST_CONFIG["reference_fasta"]).resolve()),
        output_dir = str(JOB_DIR / "12a-stats_aggregation" / "{sample}"),
        log_file = str(JOB_DIR / "logs" / "12a-stats_aggregation" / "{sample}.log"),
        global_log_file = str(GLOBAL_LOG_FILE),
        exec_script = str(SRC_DIR / "mitopipeline" / "exec" / "12a-stats_aggregation" / "run_stats_aggregation.py")
    conda:
        "../../envs/mitopipeline.yaml"
    shell:
        """
        TOP_HITS_ARGS=""
        if [ -s {params.top_hits_tsv:q} ]; then
            TOP_HITS_ARGS="--top-hits-tsv {params.top_hits_tsv:q} --biological-validation-json {params.biological_validation_json:q}"
            if [ -s {params.circular_coverage_json:q} ]; then
                TOP_HITS_ARGS="$TOP_HITS_ARGS --circular-coverage-json {params.circular_coverage_json:q}"
            fi
            TOP_HITS_ARGS="$TOP_HITS_ARGS --reference-fasta {params.reference_fasta:q}"
        fi
        ANNOTATION_ARGS=""
        if [ -s {params.annotation_result_json:q} ]; then
            ANNOTATION_ARGS="--annotation-result-json {params.annotation_result_json:q}"
        fi
        RECRUITED_QC_ARGS=""
        if [ -n "{params.recruited_reads_qc_json}" ] && [ -s "{params.recruited_reads_qc_json}" ]; then
            RECRUITED_QC_ARGS="--recruited-reads-qc-json {params.recruited_reads_qc_json:q}"
        fi
        {PYTHON_BIN:q} {params.exec_script:q} \
            --sample-id {wildcards.sample} \
            --contig-selection-json {params.contig_selection_json:q} \
            $RECRUITED_QC_ARGS \
            $TOP_HITS_ARGS $ANNOTATION_ARGS \
            --output-dir {params.output_dir:q} \
            --log-file {params.log_file:q} \
            --global-log-file {params.global_log_file:q} \
            || true
        mkdir -p $(dirname {output.done:q})
        touch {output.done:q}
        """
