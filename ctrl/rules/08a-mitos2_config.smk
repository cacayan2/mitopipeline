"""
08a-mitos2_config.smk

Ensures the local MITOS2 reference data is downloaded once (job-
level, shared across every sample -- the same race-avoidance pattern
used for the MitoFish/BLAST databases), then resolves one sample's
MITOS2 configuration against the contig contig_selection (#06)
actually chose.
"""

# Imports
import json
from pathlib import Path

MITOS2_CONFIG = config["stage_options"]["annotation"]["mitos2"]

rule mitos_database_setup:
    output:
        # mitos_database.json is NOT declared -- the exec script
        # always writes it (success or failure); the done marker
        # alone is a sufficient dependency edge.
        done = str(JOB_DIR / "08a-mitos2_config" / "mitos_database.done")
    params:
        refdir = str((PROJECT_ROOT / MITOS2_CONFIG["refdir"]).resolve()),
        refseqver = MITOS2_CONFIG.get("refseqver", "refseq89m"),
        zenodo_record = MITOS2_CONFIG["zenodo_record"],
        output_dir = str(JOB_DIR / "08a-mitos2_config"),
        log_file = str(JOB_DIR / "logs" / "08a-mitos2_config" / "mitos_database_setup.log"),
        global_log_file = str(GLOBAL_LOG_FILE),
        exec_script = str(SRC_DIR / "mitopipeline" / "exec" / "08a-mitos2_config" / "run_mitos_database_setup.py")
    conda:
        "../../envs/08b-annotation.yaml"
    shell:
        """
        {PYTHON_BIN:q} {params.exec_script:q} \
            --refdir {params.refdir:q} \
            --refseqver {params.refseqver:q} \
            --zenodo-record {params.zenodo_record:q} \
            --output-dir {params.output_dir:q} \
            --log-file {params.log_file:q} \
            --global-log-file {params.global_log_file:q} \
            || true
        touch {output.done:q}
        """

def _selected_fasta(selection_path: Path) -> str:
    selection_path = Path(selection_path)
    if not selection_path.is_file():
        return ""
    selection = json.loads(selection_path.read_text())
    return selection.get("fasta_path", "") if selection.get("selected") else ""

rule mitos2_config:
    input:
        selection = str(JOB_DIR / "06-contig_selection" / "{sample}" / "{sample}.contig_selection.json"),
        database_done = str(JOB_DIR / "08a-mitos2_config" / "mitos_database.done")
    output:
        # config_json is NOT declared -- never written when the
        # selected FASTA is missing (contig_selection found nothing
        # usable); that's a valid, non-error outcome.
        done = str(JOB_DIR / "08a-mitos2_config" / "{sample}" / "{sample}.mitos2_config.done")
    params:
        input_fasta = lambda wc, input: _selected_fasta(input.selection),
        genetic_code = MITOS2_CONFIG.get("genetic_code", 2),
        refdir = str((PROJECT_ROOT / MITOS2_CONFIG["refdir"]).resolve()),
        refseqver = MITOS2_CONFIG.get("refseqver", "refseq89m"),
        circular_flag = lambda wc: "" if MITOS2_CONFIG.get("circular", True) else "--linear",
        noplots_flag = lambda wc: "--noplots" if MITOS2_CONFIG.get("noplots", False) else "",
        best_flag = lambda wc: "--best" if MITOS2_CONFIG.get("best", False) else "",
        ncbicode_flag = lambda wc: "--ncbicode" if MITOS2_CONFIG.get("ncbicode", False) else "",
        zip_flag = lambda wc: "--zip-output" if MITOS2_CONFIG.get("zip_output", False) else "",
        annotation_output_dir = str(JOB_DIR / "08b-annotation" / "{sample}"),
        output_dir = str(JOB_DIR / "08a-mitos2_config" / "{sample}"),
        log_file = str(JOB_DIR / "logs" / "08a-mitos2_config" / "{sample}.log"),
        global_log_file = str(GLOBAL_LOG_FILE),
        exec_script = str(SRC_DIR / "mitopipeline" / "exec" / "08a-mitos2_config" / "run_mitos2_config.py")
    conda:
        "../../envs/08b-annotation.yaml"
    shell:
        """
        if [ -s {params.input_fasta:q} ]; then
            {PYTHON_BIN:q} {params.exec_script:q} \
                --sample-id {wildcards.sample} \
                --input-fasta {params.input_fasta:q} \
                --genetic-code {params.genetic_code} \
                --refdir {params.refdir:q} \
                --refseqver {params.refseqver:q} \
                {params.circular_flag} {params.noplots_flag} {params.best_flag} {params.ncbicode_flag} {params.zip_flag} \
                --annotation-output-dir {params.annotation_output_dir:q} \
                --output-dir {params.output_dir:q} \
                --log-file {params.log_file:q} \
                --global-log-file {params.global_log_file:q} \
                || true
        fi
        mkdir -p $(dirname {output.done:q})
        touch {output.done:q}
        """
