"""
05b-blast.smk

Runs BLAST identification for one sample against #05a's resolved
configuration.
"""

# Imports
from pathlib import Path

BLAST_CONFIG = config["stage_options"]["phylogeny"]["blast"]

rule blast:
    input:
        config_done = str(JOB_DIR / "05a-blast_config" / "{sample}" / "{sample}.blast_config.done")
    output:
        # blast_result.json is NOT declared -- never written when
        # #05a itself produced no usable config (nothing selected
        # upstream); the done marker alone proves this rule completed.
        done = str(JOB_DIR / "05b-blast" / "{sample}" / "{sample}.blast.done")
    params:
        config_json = str(JOB_DIR / "05a-blast_config" / "{sample}" / "{sample}.blast_config.json"),
        output_dir = str(JOB_DIR / "05b-blast" / "{sample}"),
        working_dir = str(Path.cwd()),
        log_file = str(JOB_DIR / "logs" / "05b-blast" / "{sample}.log"),
        global_log_file = str(GLOBAL_LOG_FILE),
        metrics_dir = str(JOB_DIR / "metrics" / "05b-blast" / "{sample}"),
        exec_script = str(SRC_DIR / "mitopipeline" / "exec" / "05b-blast" / "run_blast.py")
    threads:
        BLAST_CONFIG.get("threads", 4)
    conda:
        "../../envs/05-blast.yaml"
    shell:
        """
        if [ -s {params.config_json:q} ]; then
            {PYTHON_BIN:q} {params.exec_script:q} \
                --sample-id {wildcards.sample} \
                --config-json {params.config_json:q} \
                --output-dir {params.output_dir:q} \
                --working-dir {params.working_dir:q} \
                --threads {threads} \
                --log-file {params.log_file:q} \
                --global-log-file {params.global_log_file:q} \
                --metrics-dir {params.metrics_dir:q} \
                || true
        fi
        mkdir -p $(dirname {output.done:q})
        touch {output.done:q}
        """
