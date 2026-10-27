"""
05c-blast_hits.smk

Ranks/selects #05b's raw BLAST output into a top-hit TSV, then
retrieves the top hit's reference FASTA from the same local database.
"""

# Imports
from pathlib import Path

BLAST_CONFIG = config["stage_options"]["phylogeny"]["blast"]

rule blast_hit_selection:
    input:
        blast_done = str(JOB_DIR / "05b-blast" / "{sample}" / "{sample}.blast.done")
    output:
        done = str(JOB_DIR / "05c-blast_hits" / "{sample}" / "{sample}.hit_selection.done")
    params:
        blast_tsv = str(JOB_DIR / "05b-blast" / "{sample}" / "{sample}.blast.tsv"),
        maximum_hits = BLAST_CONFIG.get("maximum_matches", 6),
        output_dir = str(JOB_DIR / "05c-blast_hits" / "{sample}"),
        log_file = str(JOB_DIR / "logs" / "05c-blast_hits" / "{sample}.selection.log"),
        global_log_file = str(GLOBAL_LOG_FILE),
        exec_script = str(SRC_DIR / "mitopipeline" / "exec" / "05c-blast_hits" / "run_blast_hit_selection.py")
    conda:
        "../../envs/05-blast.yaml"
    shell:
        """
        if [ -s {params.blast_tsv:q} ]; then
            {PYTHON_BIN:q} {params.exec_script:q} \
                --sample-id {wildcards.sample} \
                --blast-tsv {params.blast_tsv:q} \
                --maximum-hits {params.maximum_hits} \
                --output-dir {params.output_dir:q} \
                --log-file {params.log_file:q} \
                --global-log-file {params.global_log_file:q} \
                || true
        fi
        mkdir -p $(dirname {output.done:q})
        touch {output.done:q}
        """

rule circular_coverage:
    input:
        selection = str(JOB_DIR / "06-contig_selection" / "{sample}" / "{sample}.contig_selection.json"),
        # Reads the BLAST database directly, not through blast_config's
        # own chain -- needs this edge declared separately.
        database_done = str(JOB_DIR / "05a-blast_config" / "blast_database.done")
    output:
        # {sample}.circular_coverage_result.json is NOT declared --
        # the exec script always writes it, including a clean
        # "not circular" record for a non-circular/unselected sample,
        # so the done marker alone is a sufficient dependency edge.
        done = str(JOB_DIR / "05c-blast_hits" / "{sample}" / "{sample}.circular_coverage.done")
    params:
        database = str((PROJECT_ROOT / BLAST_CONFIG["database"]).resolve()),
        output_dir = str(JOB_DIR / "05c-blast_hits" / "{sample}"),
        working_dir = str(Path.cwd()),
        log_file = str(JOB_DIR / "logs" / "05c-blast_hits" / "{sample}.circular_coverage.log"),
        global_log_file = str(GLOBAL_LOG_FILE),
        exec_script = str(SRC_DIR / "mitopipeline" / "exec" / "05c-blast_hits" / "run_circular_coverage.py")
    conda:
        "../../envs/05-blast.yaml"
    shell:
        """
        {PYTHON_BIN:q} {params.exec_script:q} \
            --sample-id {wildcards.sample} \
            --contig-selection-json {input.selection:q} \
            --database {params.database:q} \
            --working-dir {params.working_dir:q} \
            --output-dir {params.output_dir:q} \
            --log-file {params.log_file:q} \
            --global-log-file {params.global_log_file:q} \
            || true
        mkdir -p $(dirname {output.done:q})
        touch {output.done:q}
        """

rule blast_reference_retrieval:
    input:
        hit_selection_done = str(JOB_DIR / "05c-blast_hits" / "{sample}" / "{sample}.hit_selection.done")
    output:
        done = str(JOB_DIR / "05c-blast_hits" / "{sample}" / "{sample}.blast_hits.done")
    params:
        top_hits_tsv = str(JOB_DIR / "05c-blast_hits" / "{sample}" / "{sample}.top_hits.tsv"),
        database = str((PROJECT_ROOT / BLAST_CONFIG["database"]).resolve()),
        output_dir = str(JOB_DIR / "05c-blast_hits" / "{sample}"),
        working_dir = str(Path.cwd()),
        log_file = str(JOB_DIR / "logs" / "05c-blast_hits" / "{sample}.retrieval.log"),
        global_log_file = str(GLOBAL_LOG_FILE),
        exec_script = str(SRC_DIR / "mitopipeline" / "exec" / "05c-blast_hits" / "run_blast_reference_retrieval.py")
    conda:
        "../../envs/05-blast.yaml"
    shell:
        """
        if [ -s {params.top_hits_tsv:q} ]; then
            {PYTHON_BIN:q} {params.exec_script:q} \
                --sample-id {wildcards.sample} \
                --top-hits-tsv {params.top_hits_tsv:q} \
                --database {params.database:q} \
                --output-dir {params.output_dir:q} \
                --working-dir {params.working_dir:q} \
                --log-file {params.log_file:q} \
                --global-log-file {params.global_log_file:q} \
                || true
        fi
        touch {output.done:q}
        """
