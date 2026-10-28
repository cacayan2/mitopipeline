"""
04c-seed_sequences.smk

Ensures the local MitoFish reference database is downloaded once
(job-level, shared across every sample -- the exact pattern used
elsewhere for a shared resource no per-sample rule should race to
build in parallel), then selects each sample's own mitochondrial
seed sequence before assembly. Seed selection needs only the sample's
scientific name -- no dependency on reads, so it can run in parallel
with trimming.
"""

# Imports
import json
from pathlib import Path

SEED_CONFIG = config["stage_options"]["assembly"]["seed"]

def _mitofish_database_dir(result_path: Path, fallback: str) -> str:
    """Reads the resolved MitoFish database directory mitofish_database_setup
    actually used -- never re-derived, since the real path is version-
    subdirectoried and only mitofish_database_setup resolves the version.
    Falls back to the configured root (never populated, so seed
    selection's own existing missing-file handling naturally falls back
    to unseeded) if setup never produced a usable result."""
    result_path = Path(result_path)
    if result_path.is_file():
        result = json.loads(result_path.read_text())
        if result.get("success"):
            return result.get("database_dir", fallback)
    return fallback

rule mitofish_database_setup:
    output:
        result = str(JOB_DIR / "04c-seed_sequences" / "mitofish_database.json"),
        done = str(JOB_DIR / "04c-seed_sequences" / "mitofish_database.done")
    params:
        database_root = str((PROJECT_ROOT / SEED_CONFIG["database_dir"]).resolve()),
        version = str(SEED_CONFIG.get("version") or ""),
        output_dir = str(JOB_DIR / "04c-seed_sequences"),
        log_file = str(JOB_DIR / "logs" / "04c-seed_sequences" / "mitofish_database.log"),
        global_log_file = str(GLOBAL_LOG_FILE),
        exec_script = str(SRC_DIR / "mitopipeline" / "exec" / "04c-seed_sequences" / "run_mitofish_setup.py")
    conda:
        "../../envs/04c-assembly.yaml"
    shell:
        """
        VERSION_ARGS=""
        if [ -n "{params.version}" ]; then
            VERSION_ARGS="--version {params.version:q}"
        fi
        {PYTHON_BIN:q} {params.exec_script:q} \
            --database-root {params.database_root:q} \
            $VERSION_ARGS \
            --output-dir {params.output_dir:q} \
            --log-file {params.log_file:q} \
            --global-log-file {params.global_log_file:q} \
            || true
        touch {output.done:q}
        """

rule seed_selection:
    input:
        mitofish_database_done = str(JOB_DIR / "04c-seed_sequences" / "mitofish_database.done")
    output:
        record = str(JOB_DIR / "04c-seed_sequences" / "{sample}" / "{sample}.seed_selection.json"),
        done = str(JOB_DIR / "04c-seed_sequences" / "{sample}" / "{sample}.seed_sequences.done")
    params:
        scientific_name = lambda wc: SAMPLE_TABLE.loc[wc.sample, "organism"],
        database_dir = lambda wc, input: _mitofish_database_dir(
            Path(input.mitofish_database_done).parent / "mitofish_database.json",
            fallback = str((PROJECT_ROOT / SEED_CONFIG["database_dir"]).resolve()),
        ),
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
