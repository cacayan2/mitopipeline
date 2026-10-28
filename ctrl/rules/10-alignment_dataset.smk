"""
10-alignment_dataset.smk

Builds each sample's own real phylogenetic dataset (restoring
v0.1.0's per-sample design, not the earlier batch-level combined
FASTA): the sample's own assembled genome plus its own top-ranked
BLAST reference hits, fetched live from NCBI Entrez.
"""

# Imports
from pathlib import Path

BLAST_CONFIG = config["stage_options"]["phylogeny"]["blast"]

def _selected_fasta(selection_path: Path) -> str:
    """Reads a sample's own selected assembly FASTA path from its real contig_selection.json -- empty string if nothing was selected (the same pattern 05a-blast_config.smk/08a-mitos2_config.smk already use)."""
    import json
    try:
        selection = json.loads(Path(selection_path).read_text(encoding = "utf-8"))
    except (OSError, json.JSONDecodeError):
        return ""
    return selection.get("fasta_path", "") if selection.get("selected") else ""

rule phylogeny_dataset:
    input:
        selection = str(JOB_DIR / "06-contig_selection" / "{sample}" / "{sample}.contig_selection.json"),
        selection_done = str(JOB_DIR / "06-contig_selection" / "{sample}" / "{sample}.contig_selection.done"),
        hit_selection_done = str(JOB_DIR / "05c-blast_hits" / "{sample}" / "{sample}.hit_selection.done")
    output:
        # {sample}.phylogeny_dataset.fasta is NOT declared -- never
        # written when this sample has no selected assembly or no
        # BLAST hits, a valid, non-error outcome the done marker
        # alone proves this rule completed for.
        done = str(JOB_DIR / "10-alignment_dataset" / "{sample}" / "{sample}.phylogeny_dataset.done")
    params:
        # Empty string (no assembly selected) is handled by the exec
        # script itself, the same "missing input is a normal outcome"
        # pattern used throughout this stage -- never an empty-string
        # path error here.
        assembly_fasta = lambda wc, input: _selected_fasta(input.selection) or "/dev/null",
        top_hits_tsv = str(JOB_DIR / "05c-blast_hits" / "{sample}" / "{sample}.top_hits.tsv"),
        output_fasta = str(JOB_DIR / "10-alignment_dataset" / "{sample}" / "{sample}.phylogeny_dataset.fasta"),
        entrez_email = BLAST_CONFIG["entrez_email"],
        entrez_api_key = BLAST_CONFIG.get("entrez_api_key") or "",
        output_dir = str(JOB_DIR / "10-alignment_dataset" / "{sample}"),
        log_file = str(JOB_DIR / "logs" / "10-alignment_dataset" / "{sample}.log"),
        global_log_file = str(GLOBAL_LOG_FILE),
        exec_script = str(SRC_DIR / "mitopipeline" / "exec" / "10-alignment_dataset" / "run_build_sample_dataset.py")
    resources:
        entrez = 1
    conda:
        "../../envs/08-phylogeny.yaml"
    shell:
        """
        API_KEY_ARGS=""
        if [ -n "{params.entrez_api_key}" ]; then
            API_KEY_ARGS="--entrez-api-key {params.entrez_api_key:q}"
        fi
        {PYTHON_BIN:q} {params.exec_script:q} \
            --sample-id {wildcards.sample} \
            --assembly-fasta {params.assembly_fasta:q} \
            --top-hits-tsv {params.top_hits_tsv:q} \
            --output-fasta {params.output_fasta:q} \
            --entrez-email {params.entrez_email:q} \
            $API_KEY_ARGS \
            --output-dir {params.output_dir:q} \
            --log-file {params.log_file:q} \
            --global-log-file {params.global_log_file:q} \
            || true
        mkdir -p $(dirname {output.done:q})
        touch {output.done:q}
        """
