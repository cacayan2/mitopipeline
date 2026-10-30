"""
14b-biosample_archival.smk

Builds the run's real BioSample submission metadata table --
job-level, one table covering every sample.
"""

# Imports
import json
from pathlib import Path

BIOSAMPLE_DEFAULTS = config["stage_options"]["archival"]["biosample"]

_BIOSAMPLE_DEFAULTS_PATH = JOB_DIR / "14b-biosample_archival" / "biosample_defaults.json"
if stage_enabled("14b-biosample_archival"):
    _BIOSAMPLE_DEFAULTS_PATH.parent.mkdir(parents = True, exist_ok = True)
    _BIOSAMPLE_DEFAULTS_PATH.write_text(json.dumps(BIOSAMPLE_DEFAULTS))

rule biosample_archival:
    input:
        manifest = str(RUNTIME_MANIFEST_PATH)
    output:
        done = str(JOB_DIR / "14b-biosample_archival" / "biosample_archival.done")
    params:
        defaults_json = str(_BIOSAMPLE_DEFAULTS_PATH),
        output_dir = str(JOB_DIR / "14b-biosample_archival"),
        log_file = str(JOB_DIR / "logs" / "14b-biosample_archival" / "biosample_archival.log"),
        global_log_file = str(GLOBAL_LOG_FILE),
        exec_script = str(SRC_DIR / "mitopipeline" / "exec" / "14b-biosample_archival" / "run_biosample_metadata.py")
    conda:
        "../../envs/mitopipeline.yaml"
    shell:
        """
        {PYTHON_BIN:q} {params.exec_script:q} \
            --runtime-manifest {input.manifest:q} \
            --defaults-json {params.defaults_json:q} \
            --output-dir {params.output_dir:q} \
            --log-file {params.log_file:q} \
            --global-log-file {params.global_log_file:q} \
            || true
        mkdir -p $(dirname {output.done:q})
        touch {output.done:q}
        """
