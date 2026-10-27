"""
08b-annotation.smk

Runs MITOS2 annotation for one sample against #08a's resolved
configuration.
"""

# Imports
from pathlib import Path

rule annotation:
    input:
        config_done = str(JOB_DIR / "08a-mitos2_config" / "{sample}" / "{sample}.mitos2_config.done")
    output:
        # mitos2_annotation_result.json is NOT declared -- never
        # written when #08a itself produced no usable config.
        done = str(JOB_DIR / "08b-annotation" / "{sample}" / "{sample}.annotation.done")
    params:
        config_json = str(JOB_DIR / "08a-mitos2_config" / "{sample}" / "{sample}.mitos2_config.json"),
        working_dir = str(Path.cwd()),
        log_file = str(JOB_DIR / "logs" / "08b-annotation" / "{sample}.log"),
        global_log_file = str(GLOBAL_LOG_FILE),
        metrics_dir = str(JOB_DIR / "metrics" / "08b-annotation" / "{sample}"),
        exec_script = str(SRC_DIR / "mitopipeline" / "exec" / "08b-annotation" / "run_mitos2_annotation.py")
    conda:
        "../../envs/08b-annotation.yaml"
    shell:
        """
        if [ -s {params.config_json:q} ]; then
            {PYTHON_BIN:q} {params.exec_script:q} \
                --sample-id {wildcards.sample} \
                --config-json {params.config_json:q} \
                --working-dir {params.working_dir:q} \
                --log-file {params.log_file:q} \
                --global-log-file {params.global_log_file:q} \
                --metrics-dir {params.metrics_dir:q} \
                || true
        fi
        mkdir -p $(dirname {output.done:q})
        touch {output.done:q}
        """
