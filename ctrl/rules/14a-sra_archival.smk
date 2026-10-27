"""
14a-sra_archival.smk

Builds the run's real SRA submission metadata table (#368/#371) --
job-level, one table covering every sample.
"""

# Imports
import json
from pathlib import Path

SRA_DEFAULTS = config["stage_options"]["archival"]["sra"]

# Written once at parse time, only when actually enabled -- the exec
# script reads real JSON, not a config dict it has no access to once
# invoked as a subprocess.
_SRA_DEFAULTS_PATH = JOB_DIR / "14a-sra_archival" / "sra_defaults.json"
if stage_enabled("14a-sra_archival"):
    _SRA_DEFAULTS_PATH.parent.mkdir(parents = True, exist_ok = True)
    _SRA_DEFAULTS_PATH.write_text(json.dumps(SRA_DEFAULTS))

rule sra_archival:
    input:
        manifest = str(RUNTIME_MANIFEST_PATH)
    output:
        done = str(JOB_DIR / "14a-sra_archival" / "sra_archival.done")
    params:
        defaults_json = str(_SRA_DEFAULTS_PATH),
        output_dir = str(JOB_DIR / "14a-sra_archival"),
        log_file = str(JOB_DIR / "logs" / "14a-sra_archival" / "sra_archival.log"),
        global_log_file = str(GLOBAL_LOG_FILE),
        exec_script = str(SRC_DIR / "mitopipeline" / "exec" / "14a-sra_archival" / "run_sra_metadata.py")
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
