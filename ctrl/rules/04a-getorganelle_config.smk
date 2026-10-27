"""
04a-getorganelle_config.smk

Resolves one sample's GetOrganelle assembly configuration (organelle
type, seed/no-seed input, output path) before the real assembly runs.
"""

# Imports
from pathlib import Path

ASSEMBLY_CONFIG = config["stage_options"]["assembly"]["getorganelle"]

def getorganelle_config_input(wc):
    """
    Depends on #04c's seed-selection record when seed selection is
    enabled, otherwise has no real file dependency.
    """
    if stage_enabled("04c-seed_sequences"):
        return {"seed_record": str(JOB_DIR / "04c-seed_sequences" / "{sample}" / f"{wc.sample}.seed_selection.json")}
    return {}

rule getorganelle_config:
    input:
        unpack(getorganelle_config_input)
    output:
        config_json = str(JOB_DIR / "04a-getorganelle_config" / "{sample}" / "{sample}.getorganelle_config.json"),
        done = str(JOB_DIR / "04a-getorganelle_config" / "{sample}" / "{sample}.getorganelle_config.done")
    params:
        organelle_type = ASSEMBLY_CONFIG.get("database_type", "animal_mt"),
        assembly_output_dir = str(JOB_DIR / "04d-assembly" / "{sample}" / "primary"),
        output_dir = str(JOB_DIR / "04a-getorganelle_config" / "{sample}"),
        log_file = str(JOB_DIR / "logs" / "04a-getorganelle_config" / "{sample}.log"),
        global_log_file = str(GLOBAL_LOG_FILE),
        exec_script = str(SRC_DIR / "mitopipeline" / "exec" / "04a-getorganelle_config" / "run_getorganelle_config.py"),
        # The seed FASTA path itself (not the JSON record) -- empty when
        # seed selection didn't produce one (a fallback result, or the
        # stage is disabled), in which case GetOrganelle falls back to
        # its own default database.
        seed_fasta = lambda wc: (
            _seed_fasta_from_record(JOB_DIR / "04c-seed_sequences" / wc.sample / f"{wc.sample}.seed_selection.json")
            if stage_enabled("04c-seed_sequences") else ""
        )
    conda:
        "../../envs/04c-assembly.yaml"
    shell:
        """
        {PYTHON_BIN:q} {params.exec_script:q} \
            --sample-id {wildcards.sample} \
            --organelle-type {params.organelle_type:q} \
            {params.seed_fasta} \
            --assembly-output-dir {params.assembly_output_dir:q} \
            --output-dir {params.output_dir:q} \
            --log-file {params.log_file:q} \
            --global-log-file {params.global_log_file:q}
        touch {output.done:q}
        """

def _seed_fasta_from_record(record_path: Path) -> str:
    """Reads a seed_selection.json record and returns a --seed-file
    argument fragment, or an empty string for a fallback result."""
    import json
    record_path = Path(record_path)
    if not record_path.is_file():
        return ""
    record = json.loads(record_path.read_text())
    fasta_path = record.get("fasta_path")
    return f"--seed-file {fasta_path}" if fasta_path else ""
