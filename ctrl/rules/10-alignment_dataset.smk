"""
10-alignment_dataset.smk

Builds the combined, batch-level phylogeny FASTA (#344/#345) across
every sample in this run -- a single job-level artifact, not
per-sample.
"""

# Imports
from pathlib import Path

rule build_combined_fasta:
    input:
        # Every sample's contig_selection/biological_validation done
        # markers -- the real dependency edge; the actual JSON records
        # this reads are referenced by path inside the exec script
        # itself (same "conditionally produced" pattern used
        # throughout 06-contig_selection.smk).
        selection_done = expand(str(JOB_DIR / "06-contig_selection" / "{sample}" / "{sample}.contig_selection.done"), sample = SAMPLES),
        validation_done = expand(str(JOB_DIR / "06-contig_selection" / "{sample}" / "{sample}.biological_validation.done"), sample = SAMPLES)
    output:
        # combined.fasta is NOT declared -- never written when zero
        # samples end up accepted; that's a valid, non-error outcome
        # this rule must still complete for.
        done = str(JOB_DIR / "10-alignment_dataset" / "alignment_dataset.done")
    params:
        sample_flags = " ".join(f"--sample-id {s}" for s in SAMPLES),
        contig_selection_dir = str(JOB_DIR / "06-contig_selection"),
        output_dir = str(JOB_DIR / "10-alignment_dataset"),
        log_file = str(JOB_DIR / "logs" / "10-alignment_dataset" / "build_combined_fasta.log"),
        global_log_file = str(GLOBAL_LOG_FILE),
        exec_script = str(SRC_DIR / "mitopipeline" / "exec" / "10-alignment_dataset" / "run_build_combined_fasta.py")
    conda:
        "../../envs/09-visualization.yaml"
    shell:
        """
        {PYTHON_BIN:q} {params.exec_script:q} \
            {params.sample_flags} \
            --contig-selection-dir {params.contig_selection_dir:q} \
            --output-dir {params.output_dir:q} \
            --log-file {params.log_file:q} \
            --global-log-file {params.global_log_file:q} \
            || true
        mkdir -p $(dirname {output.done:q})
        touch {output.done:q}
        """
