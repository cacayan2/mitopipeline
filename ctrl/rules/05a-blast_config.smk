"""
05a-blast_config.smk

Ensures the local BLAST reference database is downloaded and built
once (job-level, shared across every sample -- the same race-
avoidance pattern used for the MitoFish database), then resolves one
sample's BLAST configuration against the contig contig_selection
(#06) actually chose.
"""

# Imports
import json
from pathlib import Path

BLAST_CONFIG = config["stage_options"]["phylogeny"]["blast"]

rule blast_database_setup:
    output:
        # blast_database.json is NOT declared -- never written when
        # initialize_database is False (the exec script is skipped
        # entirely); the done marker alone is a sufficient, always-
        # true dependency edge, same pattern as circular_coverage's
        # own conditionally-produced result.
        done = str(JOB_DIR / "05a-blast_config" / "blast_database.done")
    params:
        reference_fasta = str((PROJECT_ROOT / BLAST_CONFIG["reference_fasta"]).resolve()),
        database_prefix = str((PROJECT_ROOT / BLAST_CONFIG["database"]).resolve()),
        database_title = BLAST_CONFIG["database_title"],
        taxon_scope = BLAST_CONFIG["taxon_scope"],
        entrez_email = BLAST_CONFIG["entrez_email"],
        entrez_api_key = str(BLAST_CONFIG.get("entrez_api_key") or ""),
        max_records = BLAST_CONFIG.get("ncbi_max_records", 10000),
        batch_size = BLAST_CONFIG.get("ncbi_batch_size", 50),
        parse_seqids_flag = "--parse-seqids" if BLAST_CONFIG.get("parse_seqids", True) else "",
        overwrite_flag = "--overwrite" if BLAST_CONFIG.get("overwrite_database", False) else "",
        run_setup = "true" if BLAST_CONFIG.get("initialize_database", True) else "false",
        output_dir = str(JOB_DIR / "05a-blast_config"),
        log_file = str(JOB_DIR / "logs" / "05a-blast_config" / "blast_database_setup.log"),
        global_log_file = str(GLOBAL_LOG_FILE),
        exec_script = str(SRC_DIR / "mitopipeline" / "exec" / "05a-blast_config" / "run_blast_database_setup.py")
    conda:
        "../../envs/05-blast.yaml"
    shell:
        """
        if [ {params.run_setup:q} = "true" ]; then
            API_KEY_ARGS=""
            if [ -n {params.entrez_api_key:q} ]; then
                API_KEY_ARGS="--entrez-api-key {params.entrez_api_key:q}"
            fi
            {PYTHON_BIN:q} {params.exec_script:q} \
                --reference-fasta {params.reference_fasta:q} \
                --database-prefix {params.database_prefix:q} \
                --database-title {params.database_title:q} \
                --taxon-scope {params.taxon_scope:q} \
                --entrez-email {params.entrez_email:q} \
                $API_KEY_ARGS \
                --max-records {params.max_records} \
                --batch-size {params.batch_size} \
                {params.parse_seqids_flag} \
                {params.overwrite_flag} \
                --output-dir {params.output_dir:q} \
                --log-file {params.log_file:q} \
                --global-log-file {params.global_log_file:q} \
                || true
        fi
        touch {output.done:q}
        """

def _selected_fasta(selection_path: Path) -> str:
    selection_path = Path(selection_path)
    if not selection_path.is_file():
        return ""
    selection = json.loads(selection_path.read_text())
    return selection.get("fasta_path", "") if selection.get("selected") else ""

rule blast_config:
    input:
        selection = str(JOB_DIR / "06-contig_selection" / "{sample}" / "{sample}.contig_selection.json"),
        database_done = str(JOB_DIR / "05a-blast_config" / "blast_database.done")
    output:
        # config_json is NOT declared -- never written when the
        # selected FASTA is missing/empty (contig_selection found
        # nothing usable); that's a valid, non-error outcome this
        # rule must still complete for.
        done = str(JOB_DIR / "05a-blast_config" / "{sample}" / "{sample}.blast_config.done")
    params:
        query_fasta = lambda wc, input: _selected_fasta(input.selection),
        database = str((PROJECT_ROOT / BLAST_CONFIG["database"]).resolve()),
        mode = BLAST_CONFIG.get("mode", "local"),
        blast_output_dir = str(JOB_DIR / "05b-blast" / "{sample}"),
        output_dir = str(JOB_DIR / "05a-blast_config" / "{sample}"),
        log_file = str(JOB_DIR / "logs" / "05a-blast_config" / "{sample}.log"),
        global_log_file = str(GLOBAL_LOG_FILE),
        exec_script = str(SRC_DIR / "mitopipeline" / "exec" / "05a-blast_config" / "run_blast_config.py")
    conda:
        "../../envs/05-blast.yaml"
    shell:
        """
        {PYTHON_BIN:q} {params.exec_script:q} \
            --sample-id {wildcards.sample} \
            --query-fasta {params.query_fasta:q} \
            --database {params.database:q} \
            --mode {params.mode:q} \
            --blast-output-dir {params.blast_output_dir:q} \
            --output-dir {params.output_dir:q} \
            --log-file {params.log_file:q} \
            --global-log-file {params.global_log_file:q} \
            || true
        touch {output.done:q}
        """
