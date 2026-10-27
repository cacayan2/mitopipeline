"""
05a-blast_config.smk

Resolves one sample's BLAST configuration against the contig
contig_selection (#06) actually chose.
"""

# Imports
import json
from pathlib import Path

BLAST_CONFIG = config["stage_options"]["phylogeny"]["blast"]

def _selected_fasta(selection_path: Path) -> str:
    selection_path = Path(selection_path)
    if not selection_path.is_file():
        return ""
    selection = json.loads(selection_path.read_text())
    return selection.get("fasta_path", "") if selection.get("selected") else ""

rule blast_config:
    input:
        selection = str(JOB_DIR / "06-contig_selection" / "{sample}" / "{sample}.contig_selection.json")
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
