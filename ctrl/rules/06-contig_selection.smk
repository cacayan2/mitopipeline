"""
06-contig_selection.smk

Selects one sample's final mitochondrial contig across the primary
assembly and every rescue attempt, then runs biological validation
(#414) against every attempt considered -- identically, primary and
rescue alike -- wiring BLAST/MITOS2 evidence onto whichever attempt
was actually selected.
"""

# Imports
from pathlib import Path

ASSESSMENT_CONFIG = config["stage_options"]["assembly_post"]["assembly_assessment"]

def _all_attempt_dirs(sample: str) -> list[str]:
    dirs = [str(JOB_DIR / "04d-assembly" / sample / "primary")]
    # Each rescue attempt's real data.fasta/status.json live in its
    # nested reassembly/ directory (see 07-assembly_rescue.smk's
    # _reassembly_dir()), not attempt_dir itself.
    dirs += [str(_reassembly_dir(sample, int(n))) for n in RESCUE_ATTEMPT_NUMS]
    return dirs

def _attempt_dir_flags(sample: str) -> str:
    """Builds the repeated "--attempt-dir X --attempt-dir Y ..." flag
    string both exec scripts expect, in primary-then-rescue order."""
    return " ".join(f"--attempt-dir {d}" for d in _all_attempt_dirs(sample))

rule contig_selection:
    input:
        # The primary's own status.json/data.fasta are NOT declared --
        # neither is a formal output of the assembly rule (see
        # 04d-assembly.smk), so the real dependency edge is its done
        # marker; run_contig_selection.py's own _load_candidate()
        # already handles a missing status file gracefully either way.
        primary_done = str(JOB_DIR / "04d-assembly" / "{sample}" / "{sample}.assembly.done"),
        rescue_done = str(JOB_DIR / "07-assembly_rescue" / "{sample}" / "{sample}.assembly_rescue.done")
    output:
        selection = str(JOB_DIR / "06-contig_selection" / "{sample}" / "{sample}.contig_selection.json"),
        done = str(JOB_DIR / "06-contig_selection" / "{sample}" / "{sample}.contig_selection.done")
    params:
        attempt_dir_flags = lambda wc: _attempt_dir_flags(wc.sample),
        output_dir = str(JOB_DIR / "06-contig_selection" / "{sample}"),
        log_file = str(JOB_DIR / "logs" / "06-contig_selection" / "{sample}.log"),
        global_log_file = str(GLOBAL_LOG_FILE),
        exec_script = str(SRC_DIR / "mitopipeline" / "exec" / "06-contig_selection" / "run_contig_selection.py")
    conda:
        "../../envs/mitopipeline.yaml"
    shell:
        """
        {PYTHON_BIN:q} {params.exec_script:q} \
            --sample-id {wildcards.sample} \
            {params.attempt_dir_flags} \
            --output-dir {params.output_dir:q} \
            --log-file {params.log_file:q} \
            --global-log-file {params.global_log_file:q} \
            || true
        touch {output.done:q}
        """

rule biological_validation:
    input:
        selection = str(JOB_DIR / "06-contig_selection" / "{sample}" / "{sample}.contig_selection.json"),
        # blast_result.json/mitos2_annotation_result.json themselves
        # are NOT declared -- neither is a formal output of any rule
        # (both are conditionally produced), so their own ".done"
        # markers are the real dependency edges; their paths are
        # referenced in params and existence-checked in the shell.
        blast_done = str(JOB_DIR / "05b-blast" / "{sample}" / "{sample}.blast.done"),
        circular_coverage_done = str(JOB_DIR / "05c-blast_hits" / "{sample}" / "{sample}.circular_coverage.done"),
        annotation_done = str(JOB_DIR / "08b-annotation" / "{sample}" / "{sample}.annotation.done")
    output:
        validation = str(JOB_DIR / "06-contig_selection" / "{sample}" / "{sample}.biological_validation.json"),
        done = str(JOB_DIR / "06-contig_selection" / "{sample}" / "{sample}.biological_validation.done")
    params:
        attempt_dir_flags = lambda wc: _attempt_dir_flags(wc.sample),
        min_length = ASSESSMENT_CONFIG.get("complete_min_length", 14000),
        max_length = ASSESSMENT_CONFIG.get("complete_max_length", 22000),
        blast_result_json = str(JOB_DIR / "05b-blast" / "{sample}" / "{sample}.blast_result.json"),
        top_hits_tsv = str(JOB_DIR / "05c-blast_hits" / "{sample}" / "{sample}.top_hits.tsv"),
        circular_coverage_json = str(JOB_DIR / "05c-blast_hits" / "{sample}" / "{sample}.circular_coverage_result.json"),
        annotation_result_json = str(JOB_DIR / "08b-annotation" / "{sample}" / "{sample}.mitos2_annotation_result.json"),
        output_dir = str(JOB_DIR / "06-contig_selection" / "{sample}"),
        log_file = str(JOB_DIR / "logs" / "06-contig_selection" / "{sample}.validation.log"),
        global_log_file = str(GLOBAL_LOG_FILE),
        exec_script = str(SRC_DIR / "mitopipeline" / "exec" / "06-contig_selection" / "run_biological_validation.py")
    conda:
        "../../envs/mitopipeline.yaml"
    shell:
        """
        BLAST_ARGS=""
        if [ -s {params.blast_result_json:q} ]; then
            BLAST_ARGS="--blast-result-json {params.blast_result_json:q} --blast-top-hits-tsv {params.top_hits_tsv:q}"
        fi
        ANNOTATION_ARGS=""
        if [ -s {params.annotation_result_json:q} ]; then
            ANNOTATION_ARGS="--annotation-result-json {params.annotation_result_json:q}"
        fi
        CIRCULAR_COVERAGE_ARGS=""
        if [ -s {params.circular_coverage_json:q} ]; then
            CIRCULAR_COVERAGE_ARGS="--circular-coverage-json {params.circular_coverage_json:q}"
        fi
        {PYTHON_BIN:q} {params.exec_script:q} \
            --sample-id {wildcards.sample} \
            {params.attempt_dir_flags} \
            --contig-selection-json {input.selection:q} \
            --min-length {params.min_length} --max-length {params.max_length} \
            $BLAST_ARGS $ANNOTATION_ARGS $CIRCULAR_COVERAGE_ARGS \
            --output-dir {params.output_dir:q} \
            --log-file {params.log_file:q} \
            --global-log-file {params.global_log_file:q} \
            || true
        touch {output.done:q}
        """
