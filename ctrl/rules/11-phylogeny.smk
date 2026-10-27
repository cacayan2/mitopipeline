"""
11-phylogeny.smk

Aligns the combined, batch-level FASTA with MAFFT (#348) and infers
a maximum-likelihood phylogeny with IQ-TREE (#436) -- both single
job-level artifacts, not per-sample.
"""

# Imports
from pathlib import Path

MAFFT_CONFIG = config["stage_options"]["phylogeny"]["mafft"]
IQTREE_CONFIG = config["stage_options"]["phylogeny"]["iqtree"]

rule mafft_alignment:
    input:
        dataset_done = str(JOB_DIR / "10-alignment_dataset" / "alignment_dataset.done")
    output:
        # aligned.fasta is NOT declared -- never written when the
        # combined dataset itself was never produced (zero accepted
        # samples); the done marker alone proves this rule completed.
        done = str(JOB_DIR / "11-phylogeny" / "mafft_alignment.done")
    params:
        input_fasta = str(JOB_DIR / "10-alignment_dataset" / "combined.fasta"),
        output_fasta = str(JOB_DIR / "11-phylogeny" / "aligned.fasta"),
        working_dir = str(Path.cwd()),
        output_dir = str(JOB_DIR / "11-phylogeny"),
        log_file = str(JOB_DIR / "logs" / "11-phylogeny" / "mafft_alignment.log"),
        global_log_file = str(GLOBAL_LOG_FILE),
        metrics_dir = str(JOB_DIR / "metrics" / "11-phylogeny"),
        exec_script = str(SRC_DIR / "mitopipeline" / "exec" / "11-phylogeny" / "run_mafft_alignment.py")
    threads:
        MAFFT_CONFIG.get("threads", 4)
    conda:
        "../../envs/08-phylogeny.yaml"
    shell:
        """
        if [ -s {params.input_fasta:q} ]; then
            {PYTHON_BIN:q} {params.exec_script:q} \
                --run-id job \
                --input-fasta {params.input_fasta:q} \
                --output-fasta {params.output_fasta:q} \
                --strategy {MAFFT_CONFIG[strategy]:q} \
                --threads {threads} \
                --working-dir {params.working_dir:q} \
                --output-dir {params.output_dir:q} \
                --log-file {params.log_file:q} \
                --global-log-file {params.global_log_file:q} \
                --metrics-dir {params.metrics_dir:q} \
                || true
        fi
        mkdir -p $(dirname {output.done:q})
        touch {output.done:q}
        """

rule iqtree_inference:
    input:
        mafft_done = str(JOB_DIR / "11-phylogeny" / "mafft_alignment.done")
    output:
        done = str(JOB_DIR / "11-phylogeny" / "phylogeny.done")
    params:
        alignment = str(JOB_DIR / "11-phylogeny" / "aligned.fasta"),
        # IQ-TREE's own -o target (output_prefix's parent) must stay
        # fresh the same way GetOrganelle's did -- nested one level
        # under the rule's own output_dir, never equal to the
        # directory a declared output's parent would pre-create.
        output_prefix = str(JOB_DIR / "11-phylogeny" / "tree" / "job"),
        output_newick = str(JOB_DIR / "11-phylogeny" / "job.treefile.nwk"),
        working_dir = str(Path.cwd()),
        output_dir = str(JOB_DIR / "11-phylogeny"),
        log_file = str(JOB_DIR / "logs" / "11-phylogeny" / "iqtree.log"),
        global_log_file = str(GLOBAL_LOG_FILE),
        metrics_dir = str(JOB_DIR / "metrics" / "11-phylogeny"),
        exec_script = str(SRC_DIR / "mitopipeline" / "exec" / "11-phylogeny" / "run_iqtree.py")
    threads:
        IQTREE_CONFIG.get("threads", 4)
    conda:
        "../../envs/08-phylogeny.yaml"
    shell:
        """
        if [ -s {params.alignment:q} ]; then
            {PYTHON_BIN:q} {params.exec_script:q} \
                --run-id job \
                --alignment {params.alignment:q} \
                --output-prefix {params.output_prefix:q} \
                --output-newick {params.output_newick:q} \
                --sequence-type {IQTREE_CONFIG[sequence_type]:q} \
                --model {IQTREE_CONFIG[model]:q} \
                --ultrafast-bootstrap {IQTREE_CONFIG[ultrafast_bootstrap]} \
                --sh-alrt {IQTREE_CONFIG[sh_alrt]} \
                --random-seed {IQTREE_CONFIG[random_seed]} \
                --threads {threads} \
                --working-dir {params.working_dir:q} \
                --output-dir {params.output_dir:q} \
                --log-file {params.log_file:q} \
                --global-log-file {params.global_log_file:q} \
                --metrics-dir {params.metrics_dir:q} \
                || true
        fi
        mkdir -p $(dirname {output.done:q})
        touch {output.done:q}
        """
