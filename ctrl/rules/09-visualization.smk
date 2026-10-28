"""
09-visualization.smk

Renders one sample's circular mitochondrial genome map against the
same selected assembly attempt and MITOS2 annotation #08b already
produced.
"""

# Imports
import json
from pathlib import Path

def _selected_fasta(selection_path: Path) -> str:
    selection_path = Path(selection_path)
    if not selection_path.is_file():
        return ""
    selection = json.loads(selection_path.read_text())
    return selection.get("fasta_path", "") if selection.get("selected") else ""

def _annotation_gff(sample: str) -> str:
    """Reads a sample's own real MITOS2 annotation GFF path from its real mitos2_annotation_result.json -- never a hardcoded filename, since Mitos2Runner.postprocess_outputs() renames runmitos.py's own fixed "result.gff" to a sample-prefixed name."""
    result_path = JOB_DIR / "08b-annotation" / sample / f"{sample}.mitos2_annotation_result.json"
    if not result_path.is_file():
        return ""
    result = json.loads(result_path.read_text())
    return result.get("gff_path", "") if result.get("success") else ""

rule visualization:
    input:
        selection = str(JOB_DIR / "06-contig_selection" / "{sample}" / "{sample}.contig_selection.json"),
        annotation_done = str(JOB_DIR / "08b-annotation" / "{sample}" / "{sample}.annotation.done")
    output:
        # circular_map_result.json is NOT declared -- never written
        # when the selected FASTA or MITOS2's own GFF is missing.
        done = str(JOB_DIR / "09-visualization" / "{sample}" / "{sample}.visualization.done")
    params:
        fasta_path = lambda wc, input: _selected_fasta(input.selection),
        gff_path = lambda wc: _annotation_gff(wc.sample),
        dpi = config["stage_options"]["annotation_viz"]["dpi"],
        output_dir = str(JOB_DIR / "09-visualization" / "{sample}"),
        log_file = str(JOB_DIR / "logs" / "09-visualization" / "{sample}.log"),
        global_log_file = str(GLOBAL_LOG_FILE),
        exec_script = str(SRC_DIR / "mitopipeline" / "exec" / "09-visualization" / "run_circular_map.py")
    conda:
        "../../envs/09-visualization.yaml"
    shell:
        """
        if [ -s {params.fasta_path:q} ] && [ -s {params.gff_path:q} ]; then
            {PYTHON_BIN:q} {params.exec_script:q} \
                --sample-id {wildcards.sample} \
                --gff-path {params.gff_path:q} \
                --fasta-path {params.fasta_path:q} \
                --dpi {params.dpi} \
                --output-dir {params.output_dir:q} \
                --log-file {params.log_file:q} \
                --global-log-file {params.global_log_file:q} \
                || true
        fi
        mkdir -p $(dirname {output.done:q})
        touch {output.done:q}
        """
