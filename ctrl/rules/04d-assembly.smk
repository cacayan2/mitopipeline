"""
04d-assembly.smk

Runs the primary GetOrganelle assembly attempt and classifies its
result.
"""

# Imports
from pathlib import Path

ASSEMBLY_CONFIG = config["stage_options"]["assembly"]["getorganelle"]

rule assembly:
    input:
        r1 = str(JOB_DIR / "trimming" / "{sample}" / "{sample}.R1.trimmed.fastq.gz"),
        r2 = str(JOB_DIR / "trimming" / "{sample}" / "{sample}.R2.trimmed.fastq.gz"),
        config_json = str(JOB_DIR / "04a-getorganelle_config" / "{sample}" / "{sample}.getorganelle_config.json")
    output:
        # Neither data.fasta NOR status.json is declared here, and
        # deliberately so -- confirmed by a real failure: Snakemake
        # auto-creates the PARENT directory of every declared output
        # before the rule runs, and status.json's parent is exactly
        # "primary/", the same leaf GetOrganelleRunner passes to
        # GetOrganelle as -o. GetOrganelle itself refuses to run if
        # that directory already exists (even empty) -- so declaring
        # status.json here made Snakemake itself sabotage every real
        # assembly attempt before it could start. done lives one
        # level up, in a directory Snakemake pre-creating is harmless.
        done = str(JOB_DIR / "04d-assembly" / "{sample}" / "{sample}.assembly.done")
    params:
        working_dir = str(Path.cwd()),
        log_file = str(JOB_DIR / "logs" / "04d-assembly" / "{sample}.log"),
        global_log_file = str(GLOBAL_LOG_FILE),
        metrics_dir = str(JOB_DIR / "metrics" / "04d-assembly" / "{sample}"),
        exec_script = str(SRC_DIR / "mitopipeline" / "exec" / "04d-assembly" / "run_getorganelle_assembly.py")
    threads:
        ASSEMBLY_CONFIG.get("threads", 4)
    conda:
        "../../envs/04c-assembly.yaml"
    shell:
        """
        {PYTHON_BIN:q} {params.exec_script:q} \
            --sample-id {wildcards.sample} \
            --r1 {input.r1:q} \
            --r2 {input.r2:q} \
            --config-json {input.config_json:q} \
            --working-dir {params.working_dir:q} \
            --log-file {params.log_file:q} \
            --global-log-file {params.global_log_file:q} \
            --threads {threads} \
            --metrics-dir {params.metrics_dir:q} \
            || true
        touch {output.done:q}
        """
