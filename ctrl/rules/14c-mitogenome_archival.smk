"""
14c-mitogenome_archival.smk

Prepares one sample's real GenBank submission folder (#370/#373) --
per-sample, since each one depends only on that sample's own
selected assembly and MITOS2 annotation.
"""

# Imports
import json
from pathlib import Path

GENBANK_DEFAULTS = config["stage_options"]["archival"]["biosample"]

_GENBANK_DEFAULTS_PATH = JOB_DIR / "14c-mitogenome_archival" / "genbank_defaults.json"
if stage_enabled("14c-mitogenome_archival"):
    _GENBANK_DEFAULTS_PATH.parent.mkdir(parents = True, exist_ok = True)
    _GENBANK_DEFAULTS_PATH.write_text(json.dumps(GENBANK_DEFAULTS))

rule mitogenome_archival:
    input:
        contig_selection_done = str(JOB_DIR / "06-contig_selection" / "{sample}" / "{sample}.contig_selection.done"),
        annotation_done = str(JOB_DIR / "08b-annotation" / "{sample}" / "{sample}.annotation.done")
    output:
        done = str(JOB_DIR / "14c-mitogenome_archival" / "{sample}" / "{sample}.mitogenome_archival.done")
    params:
        runtime_manifest = str(RUNTIME_MANIFEST_PATH),
        contig_selection_json = str(JOB_DIR / "06-contig_selection" / "{sample}" / "{sample}.contig_selection.json"),
        annotation_result_json = str(JOB_DIR / "08b-annotation" / "{sample}" / "{sample}.mitos2_annotation_result.json"),
        defaults_json = str(_GENBANK_DEFAULTS_PATH),
        output_dir = str(JOB_DIR / "14c-mitogenome_archival"),
        log_file = str(JOB_DIR / "logs" / "14c-mitogenome_archival" / "{sample}.log"),
        global_log_file = str(GLOBAL_LOG_FILE),
        exec_script = str(SRC_DIR / "mitopipeline" / "exec" / "14c-mitogenome_archival" / "run_mitogenome_metadata.py")
    conda:
        "../../envs/mitopipeline.yaml"
    shell:
        """
        {PYTHON_BIN:q} {params.exec_script:q} \
            --sample-id {wildcards.sample} \
            --runtime-manifest {params.runtime_manifest:q} \
            --contig-selection-json {params.contig_selection_json:q} \
            --annotation-result-json {params.annotation_result_json:q} \
            --defaults-json {params.defaults_json:q} \
            --output-dir {params.output_dir:q} \
            --log-file {params.log_file:q} \
            --global-log-file {params.global_log_file:q} \
            || true
        mkdir -p $(dirname {output.done:q})
        touch {output.done:q}
        """
