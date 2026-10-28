"""
11-phylogeny.smk

Aligns each sample's own real phylogenetic dataset with MAFFT (#348)
and infers a maximum-likelihood phylogeny with IQ-TREE (#436) --
restored to v0.1.0's per-sample design: a sample's own tree (its
assembled genome plus its own live-fetched BLAST references) never
depends on how many other samples are in the same run.
"""

# Imports
from pathlib import Path

MAFFT_CONFIG = config["stage_options"]["phylogeny"]["mafft"]
IQTREE_CONFIG = config["stage_options"]["phylogeny"]["iqtree"]

rule mafft_alignment:
    input:
        dataset_done = str(JOB_DIR / "10-alignment_dataset" / "{sample}" / "{sample}.phylogeny_dataset.done")
    output:
        # {sample}.aligned.fasta is NOT declared -- never written when
        # this sample's own phylogeny dataset was never produced (no
        # selected assembly, or no BLAST hits); the done marker alone
        # proves this rule completed.
        done = str(JOB_DIR / "11-phylogeny" / "{sample}" / "{sample}.mafft_alignment.done")
    params:
        input_fasta = str(JOB_DIR / "10-alignment_dataset" / "{sample}" / "{sample}.phylogeny_dataset.fasta"),
        output_fasta = str(JOB_DIR / "11-phylogeny" / "{sample}" / "{sample}.aligned.fasta"),
        working_dir = str(Path.cwd()),
        output_dir = str(JOB_DIR / "11-phylogeny" / "{sample}"),
        log_file = str(JOB_DIR / "logs" / "11-phylogeny" / "{sample}.mafft_alignment.log"),
        global_log_file = str(GLOBAL_LOG_FILE),
        metrics_dir = str(JOB_DIR / "metrics" / "11-phylogeny" / "{sample}"),
        exec_script = str(SRC_DIR / "mitopipeline" / "exec" / "11-phylogeny" / "run_mafft_alignment.py")
    threads:
        MAFFT_CONFIG.get("threads", 4)
    conda:
        "../../envs/08-phylogeny.yaml"
    shell:
        """
        if [ -s {params.input_fasta:q} ]; then
            {PYTHON_BIN:q} {params.exec_script:q} \
                --run-id {wildcards.sample} \
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
        mafft_done = str(JOB_DIR / "11-phylogeny" / "{sample}" / "{sample}.mafft_alignment.done")
    output:
        done = str(JOB_DIR / "11-phylogeny" / "{sample}" / "{sample}.phylogeny.done")
    params:
        alignment = str(JOB_DIR / "11-phylogeny" / "{sample}" / "{sample}.aligned.fasta"),
        # IQ-TREE's own -o target (output_prefix's parent) must stay
        # fresh the same way GetOrganelle's did -- nested one level
        # under the rule's own output_dir, never equal to the
        # directory a declared output's parent would pre-create.
        output_prefix = str(JOB_DIR / "11-phylogeny" / "{sample}" / "tree" / "{sample}"),
        output_newick = str(JOB_DIR / "11-phylogeny" / "{sample}" / "{sample}.treefile.nwk"),
        working_dir = str(Path.cwd()),
        output_dir = str(JOB_DIR / "11-phylogeny" / "{sample}"),
        log_file = str(JOB_DIR / "logs" / "11-phylogeny" / "{sample}.iqtree.log"),
        global_log_file = str(GLOBAL_LOG_FILE),
        metrics_dir = str(JOB_DIR / "metrics" / "11-phylogeny" / "{sample}"),
        exec_script = str(SRC_DIR / "mitopipeline" / "exec" / "11-phylogeny" / "run_iqtree.py")
    threads:
        IQTREE_CONFIG.get("threads", 4)
    conda:
        "../../envs/08-phylogeny.yaml"
    shell:
        """
        if [ -s {params.alignment:q} ]; then
            {PYTHON_BIN:q} {params.exec_script:q} \
                --run-id {wildcards.sample} \
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

rule tree_figure:
    input:
        phylogeny_done = str(JOB_DIR / "11-phylogeny" / "{sample}" / "{sample}.phylogeny.done")
    output:
        # {sample}.tree_figure.png/tree_figure_result.json are NOT
        # declared -- never written when IQ-TREE itself produced no
        # usable tree for this sample; the done marker alone proves
        # this rule completed.
        done = str(JOB_DIR / "11-phylogeny" / "{sample}" / "{sample}.tree_figure.done")
    params:
        iqtree_result_json = str(JOB_DIR / "11-phylogeny" / "{sample}" / "iqtree_result.json"),
        dpi = IQTREE_CONFIG.get("render_dpi", 600),
        midpoint_flag = "--midpoint-root" if IQTREE_CONFIG.get("midpoint_root_figure", True) else "",
        output_dir = str(JOB_DIR / "11-phylogeny" / "{sample}"),
        log_file = str(JOB_DIR / "logs" / "11-phylogeny" / "{sample}.tree_figure.log"),
        global_log_file = str(GLOBAL_LOG_FILE),
        exec_script = str(SRC_DIR / "mitopipeline" / "exec" / "11-phylogeny" / "run_tree_figure.py")
    conda:
        "../../envs/08-phylogeny.yaml"
    shell:
        """
        {PYTHON_BIN:q} {params.exec_script:q} \
            --iqtree-result-json {params.iqtree_result_json:q} \
            --dpi {params.dpi} \
            {params.midpoint_flag} \
            --output-dir {params.output_dir:q} \
            --log-file {params.log_file:q} \
            --global-log-file {params.global_log_file:q} \
            || true
        mkdir -p $(dirname {output.done:q})
        touch {output.done:q}
        """
