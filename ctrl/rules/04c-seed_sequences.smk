"""
04c-seed_sequences.smk

Selects one sample's mitochondrial seed sequence before assembly.
Needs only the sample's scientific name -- no dependency on reads,
so it can run in parallel with trimming.
"""

# Imports
from pathlib import Path

SEED_CONFIG = config["stage_options"]["assembly"]["seed"]

rule seed_selection:
    output:
        record = str(JOB_DIR / "04c-seed_sequences" / "{sample}" / "{sample}.seed_selection.json"),
        done = str(JOB_DIR / "04c-seed_sequences" / "{sample}" / "{sample}.seed_sequences.done")
    params:
        scientific_name = lambda wc: SAMPLE_TABLE.loc[wc.sample, "organism"],
        database_dir = str((PROJECT_ROOT / SEED_CONFIG["database_dir"]).resolve()),
        output_dir = str(JOB_DIR / "04c-seed_sequences" / "{sample}"),
        log_file = str(JOB_DIR / "logs" / "04c-seed_sequences" / "{sample}.log"),
        global_log_file = str(GLOBAL_LOG_FILE),
        exec_script = str(SRC_DIR / "mitopipeline" / "exec" / "04c-seed_sequences" / "run_seed_selection.py")
    conda:
        "../../envs/04c-assembly.yaml"
    shell:
        """
        {PYTHON_BIN:q} {params.exec_script:q} \
            --sample-id {wildcards.sample} \
            --scientific-name {params.scientific_name:q} \
            --database-dir {params.database_dir:q} \
            --output-dir {params.output_dir:q} \
            --log-file {params.log_file:q} \
            --global-log-file {params.global_log_file:q}
        touch {output.done:q}
        """
