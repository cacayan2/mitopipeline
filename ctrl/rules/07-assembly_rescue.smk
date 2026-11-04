"""
07-assembly_rescue.smk

Runs assembly rescue for one sample: a bounded sequence of mapping +
reassembly attempts, one per configured getorganelle_profiles entry,
each depending on the PREVIOUS attempt's own classification (the
primary assembly for attempt 1). Each attempt's own exec script
already decides for itself whether rescue is actually warranted
(decide_rescue()) and cheaply skips the expensive tools when it
isn't -- this file's job is only to wire that bounded chain into the
DAG, never to re-decide it.
"""

# Imports
import json
from pathlib import Path

RESCUE_CONFIG = config["stage_options"]["assembly_post"]["assembly_rescue"]
RESCUE_PROFILES = RESCUE_CONFIG.get("getorganelle_profiles") or []
RESCUE_ATTEMPT_NUMS = [str(i) for i in range(1, len(RESCUE_PROFILES) + 1)]
TRIGGER_STATUSES = ",".join(RESCUE_CONFIG.get("trigger_statuses") or [])

wildcard_constraints:
    attempt_num = r"\d+"

def _attempt_dir(sample: str, attempt_num: int) -> Path:
    return JOB_DIR / "07-assembly_rescue" / sample / f"attempt_{attempt_num}"

def _reassembly_dir(sample: str, attempt_num: int) -> Path:
    """The actual GetOrganelle -o target for a rescue reassembly --
    deliberately nested one level under _attempt_dir(), never equal
    to it. Confirmed by a real failure (the same one that hit the
    primary assembly, see 04d-assembly.smk): Snakemake auto-creates
    the parent of every declared output before a rule runs, and
    GetOrganelle refuses to use a -o directory that already exists.
    attempt_dir/ itself is safe to pre-exist (mapping's own artifacts
    already live there); this nested leaf is not."""
    return _attempt_dir(sample, attempt_num) / "reassembly"

def _previous_status_and_fasta(sample: str, attempt_num: int) -> tuple[Path, Path]:
    """The previous attempt's own status.json and data.fasta -- the
    primary assembly's for attempt 1, the prior rescue attempt's for
    every attempt after that."""
    if attempt_num == 1:
        primary_dir = JOB_DIR / "04d-assembly" / sample / "primary"
        return primary_dir / f"{sample}.assembly_status.json", primary_dir / "data.fasta"
    previous_dir = _attempt_dir(sample, attempt_num - 1)
    return previous_dir / f"{sample}.assembly_status.json", previous_dir / "data.fasta"

def _read_status(status_path: Path) -> str:
    try:
        return json.loads(Path(status_path).read_text())["status"]
    except (OSError, json.JSONDecodeError, KeyError):
        return "missing"

def _previous_done_marker(sample: str, attempt_num: int) -> Path:
    """The real, always-produced dependency edge for the PREVIOUS
    attempt -- its own status.json/data.fasta are conditionally
    produced (never declared as formal outputs), so the chain must
    depend on a done marker that's guaranteed to exist either way."""
    if attempt_num == 1:
        return JOB_DIR / "04d-assembly" / sample / f"{sample}.assembly.done"
    return _attempt_dir(sample, attempt_num - 1) / f"{sample}.rescue_reassembly.done"

def _primary_seed_info(sample: str) -> tuple[str, str]:
    """The primary attempt's own real seed_selection.json -- its
    fasta_path (empty string when none/fallback) and gene (empty
    string when none) -- so a rescue profile can either reuse it
    as-is (reseed: false) or walk past it to a genuinely different
    locus (reseed: true), never re-deriving what select_seed() for
    the primary attempt already resolved."""
    path = JOB_DIR / "04c-seed_sequences" / sample / f"{sample}.seed_selection.json"
    try:
        record = json.loads(Path(path).read_text())
    except (OSError, json.JSONDecodeError):
        return "", ""
    if record.get("fallback_used") or not record.get("fasta_path"):
        return "", ""
    return record["fasta_path"], record.get("gene") or ""

rule rescue_mapping:
    input:
        r1 = str(JOB_DIR / "02-trimming" / "{sample}" / "{sample}.R1.trimmed.fastq.gz"),
        r2 = str(JOB_DIR / "02-trimming" / "{sample}" / "{sample}.R2.trimmed.fastq.gz"),
        previous_done = lambda wc: str(_previous_done_marker(wc.sample, int(wc.attempt_num)))
    output:
        # rescued_r1/r2 are NOT declared -- only ever produced when
        # decide_rescue() actually warrants rescue; a skipped attempt
        # (already circular, or trigger status not matched) legitimately
        # produces neither, which is a valid outcome, not a failure.
        decision = str(_attempt_dir("{sample}", "{attempt_num}") / "{sample}.rescue_mapping.json"),
        done = str(_attempt_dir("{sample}", "{attempt_num}") / "{sample}.rescue_mapping.done")
    params:
        attempt_idx = lambda wc: int(wc.attempt_num) - 1,   # decide_rescue() is 0-indexed
        status = lambda wc: _read_status(_previous_status_and_fasta(wc.sample, int(wc.attempt_num))[0]),
        reference_fasta = lambda wc: str(_previous_status_and_fasta(wc.sample, int(wc.attempt_num))[1]),
        max_attempts = len(RESCUE_PROFILES),
        trigger_statuses = TRIGGER_STATUSES,
        recruitment = lambda wc: RESCUE_PROFILES[int(wc.attempt_num) - 1].get("recruitment", "candidate_mapped"),
        output_dir = lambda wc: str(_attempt_dir(wc.sample, int(wc.attempt_num))),
        working_dir = str(Path.cwd()),
        log_file = str(JOB_DIR / "logs" / "07-assembly_rescue" / "{sample}.mapping.attempt{attempt_num}.log"),
        global_log_file = str(GLOBAL_LOG_FILE),
        metrics_dir = str(JOB_DIR / "metrics" / "07-assembly_rescue" / "{sample}"),
        exec_script = str(SRC_DIR / "mitopipeline" / "exec" / "07-assembly_rescue" / "run_rescue_mapping.py")
    threads:
        RESCUE_CONFIG.get("threads", 4)
    conda:
        "../../envs/07-assembly_rescue.yaml"
    shell:
        """
        {PYTHON_BIN:q} {params.exec_script:q} \
            --sample-id {wildcards.sample} \
            --r1 {input.r1:q} --r2 {input.r2:q} \
            --reference-fasta {params.reference_fasta:q} \
            --status {params.status:q} \
            --attempt {params.attempt_idx} \
            --trigger-statuses {params.trigger_statuses:q} \
            --max-attempts {params.max_attempts} \
            --recruitment {params.recruitment:q} \
            --output-dir {params.output_dir:q} \
            --working-dir {params.working_dir:q} \
            --threads {threads} \
            --log-file {params.log_file:q} \
            --global-log-file {params.global_log_file:q} \
            --metrics-dir {params.metrics_dir:q} \
            || true
        touch {output.done:q}
        """

rule rescue_reassembly:
    input:
        mapping_done = str(_attempt_dir("{sample}", "{attempt_num}") / "{sample}.rescue_mapping.done")
    output:
        # data.fasta/assembly_status.json are NOT declared -- the
        # reassembly step only runs real GetOrganelle when mapping
        # actually produced rescued reads; a skipped attempt still
        # completes this rule, just with nothing to show for it.
        done = str(_attempt_dir("{sample}", "{attempt_num}") / "{sample}.rescue_reassembly.done")
    params:
        reassembly_dir = lambda wc: str(_reassembly_dir(wc.sample, int(wc.attempt_num))),
        rescued_r1 = lambda wc: str(_attempt_dir(wc.sample, int(wc.attempt_num)) / f"{wc.sample}.rescue.R1.fastq"),
        rescued_r2 = lambda wc: str(_attempt_dir(wc.sample, int(wc.attempt_num)) / f"{wc.sample}.rescue.R2.fastq"),
        profile = lambda wc: RESCUE_PROFILES[int(wc.attempt_num) - 1],
        reseed_flag = lambda wc: "1" if RESCUE_PROFILES[int(wc.attempt_num) - 1].get("reseed") else "",
        # A profile may override the job's seed mode for its reseed (e.g.
        # retry a whole-genome seed, or GetOrganelle's bundled set).
        reseed_mode = lambda wc: str(RESCUE_PROFILES[int(wc.attempt_num) - 1].get("seed_mode") or SEED_MODE),
        primary_seed_fasta = lambda wc: _primary_seed_info(wc.sample)[0],
        primary_seed_gene = lambda wc: _primary_seed_info(wc.sample)[1],
        scientific_name = lambda wc: SAMPLE_TABLE.loc[wc.sample, "organism"],
        database_dir = lambda wc: _mitofish_database_dir(
            JOB_DIR / "04c-seed_sequences" / "mitofish_database.json",
            fallback = str((PROJECT_ROOT / SEED_CONFIG["database_dir"]).resolve()),
        ),
        working_dir = str(Path.cwd()),
        log_file = str(JOB_DIR / "logs" / "07-assembly_rescue" / "{sample}.reassembly.attempt{attempt_num}.log"),
        global_log_file = str(GLOBAL_LOG_FILE),
        metrics_dir = str(JOB_DIR / "metrics" / "07-assembly_rescue" / "{sample}"),
        exec_script = str(SRC_DIR / "mitopipeline" / "exec" / "07-assembly_rescue" / "run_rescue_reassembly.py")
    threads:
        RESCUE_CONFIG.get("threads", 4)
    conda:
        "../../envs/04c-assembly.yaml"
    shell:
        """
        if [ -s {params.rescued_r1:q} ] && [ -s {params.rescued_r2:q} ]; then
            SEED_ARGS=""
            if [ -n "{params.reseed_flag}" ]; then
                SEED_ARGS="--reseed --exclude-gene {params.primary_seed_gene:q} --scientific-name {params.scientific_name:q} --database-dir {params.database_dir:q} --seed-mode {params.reseed_mode:q}"
            elif [ -n "{params.primary_seed_fasta}" ]; then
                SEED_ARGS="--seed-file {params.primary_seed_fasta:q}"
            fi
            {PYTHON_BIN:q} {params.exec_script:q} \
                --sample-id {wildcards.sample} \
                --r1 {params.rescued_r1:q} --r2 {params.rescued_r2:q} \
                --profile-name {params.profile[name]:q} \
                --kmer {params.profile[kmer]:q} \
                --target-coverage {params.profile[target_coverage]} \
                --max-rounds {params.profile[max_rounds]} \
                $SEED_ARGS \
                --output-dir {params.reassembly_dir:q} \
                --working-dir {params.working_dir:q} \
                --threads {threads} \
                --log-file {params.log_file:q} \
                --global-log-file {params.global_log_file:q} \
                --metrics-dir {params.metrics_dir:q} \
                || true
        fi
        touch {output.done:q}
        """

rule assembly_rescue_summary:
    input:
        expand(str(_attempt_dir("{{sample}}", "{attempt_num}") / "{{sample}}.rescue_reassembly.done"),
               attempt_num = RESCUE_ATTEMPT_NUMS) if RESCUE_ATTEMPT_NUMS else []
    output:
        done = str(JOB_DIR / "07-assembly_rescue" / "{sample}" / "{sample}.assembly_rescue.done")
    shell:
        "mkdir -p $(dirname {output.done:q}) && touch {output.done:q}"
