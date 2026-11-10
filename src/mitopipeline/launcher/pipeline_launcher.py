"""
pipeline_launcher.py

Launch a pipeline run with the specified configuration and parameters from the config file.

The launcher creates or reuses job-specific runtime manifest 
(outputs/<job_id>/validated_samples.tsv) and writes a runtime config
that points Snakemake to that manifest for the duration of the run.
"""

# Imports
from __future__ import  annotations
import argparse
import logging
import os
from pathlib import Path
import re
import shutil
import subprocess
import pandas as pd
import yaml
import sys
from mitopipeline.utils.manifest.manifest_parser import parse_sample_manifest
from mitopipeline.models.pipeline_job.pipeline_job import PipelineJob
from mitopipeline.models.sample.sample import Sample
from mitopipeline.utils.system.resources import (
    cgroup_memory_max_mb,
    current_cgroup_path,
    oom_monitor_threshold_mb,
    recommended_concurrent_slots,
)

# Set on the environment of the re-launched (wrapped) child process by
# _reexec_under_memory_cgroup, so that child's own main() knows not to
# wrap itself again -- the outer, unwrapped invocation never has this set.
_CGROUP_WRAPPED_ENV = "MITOPIPELINE_CGROUP_WRAPPED"

# Matches a paired-end FASTQ filename and captures the sample ID (everything
# before the read-pair token) and the mate number. Recognises Illumina-style
# "<sample>_R1_001.fastq.gz" / "<sample>_R1.fastq.gz" and SRA-style
# "<sample>_1.fastq.gz"; see discover_fastq_pairs().
FASTQ_PAIR_PATTERN = re.compile(
    r"^(?P<sample>.+?)_R?(?P<mate>[12])(?:_001)?\.(?:fastq|fq)(?:\.gz)?$"
)

def parse_args() -> argparse.Namespace:
    """
    Adds arguments to the command line parser.

    Args:
        None
    
    Returns:
        argparse.Namespace: The parsed arguments.
    """
    # Creating the command line parser. 
    parser = argparse.ArgumentParser(description = "Launches a pipeline run.")

    # Adding each argument to the parser: 
    # --config: The path to the config file.
    # --snakefile: The path to the snakefile (optional).
    # --cores: The number of cores to use for the job (optional, but 1 is default so recommend changing this).
    # --no-conda: The flag to disable the conda environment (optional, but recommend not messing with this.
    # --dry-run: The flag to run the pipeline in dry-run mode (optional).
    # --printshellcmds: The flag to print shell commands (optional - useful for debugging).
    # --rerun-triggers: Passed straight through to Snakemake (optional). "mtime" alone is how
    #     to resume a run after its output directories were moved/renamed without Snakemake
    #     rerunning every already-finished job just because its recorded input paths changed.
    # The minimum command to run the pipeline is:
    ## python pipeline_launcher.py --config <path/to/config.yaml>
    # Recommended:
    ## python pipeline_launcher.py --config <path/to/config,yaml> --cores <number>
    parser.add_argument("--config", type = Path, required = True, help = "The path to the config file.")
    parser.add_argument("--snakefile", type = Path, default = "ctrl/Snakefile", help = "The path to the snakefile (optional).")
    parser.add_argument("--cores", type = int, default = 1, help = "The number of cores to use for the job (optional, but 1 is default so recommend changing this).")
    parser.add_argument("--no-conda", action = "store_true", help = "The flag to disable the conda environment (optional, but recommend not messing with this).")
    parser.add_argument("--dry-run", action = "store_true", help = "The flag to run the pipeline in dry-run mode (optional).")
    parser.add_argument("--printshellcmds", action = "store_true", help = "The flag to print shell commands (optional - useful for debugging).")
    parser.add_argument("--rerun-triggers", nargs = "+", choices = ["mtime", "params", "input", "software-env", "code"], default = None,
                        help = "Passed through to Snakemake's own --rerun-triggers (optional). Defaults to Snakemake's own default (all triggers).")

    # Returning the parsed arguments.
    return parser.parse_args()

def load_config(config_path: Path) -> dict:
    """
    Loads a YAML configuration file. Runs before any PipelineJob
    exists -- the job itself can't be created until the config's own
    output_root/job_id are known -- so this logs through the module
    logger directly rather than a job logger.

    Args:
        config_path (Path): The path to the YAML configuration file.

    Returns:
        dict: The loaded configuration as a dictionary.
    """
    # Normalizing the config path and validating that it exists.
    config_path = Path(config_path)
    if not config_path.exists() or not config_path.is_file():
        logging.error(f"Config file does not exist: {config_path}.")
        sys.exit(1)

    # Loading the config file and validating that it is not empty.
    with open(config_path, "r", encoding = "utf-8") as handle:
        config = yaml.safe_load(handle)
    if config is None:
        logging.error(f"Config file is empty: {config_path}.")
        sys.exit(1)

    # Returning the loaded config.
    return config

def _job_configuration(config: dict) -> dict:
    """
    Reads the config's own job_configuration block -- output_root,
    manifest, input_directory, job_id, validate_input_files -- the
    same block every config file documents under its own
    "job_configuration:" heading. Falls back to the config itself
    when that block is absent, so a flat config still works.

    Args:
        config (dict): The pipeline configuration.

    Returns:
        dict: The job_configuration block, or config itself as a fallback.
    """
    return config.get("job_configuration") or config

def discover_fastq_pairs(input_dir: Path, job: PipelineJob) -> list[dict[str, str]]:
    """
    Discover paired-end FASTQ files from an input directory if a manifest with a list of these files is not provided.

    Args:
        input_dir (Path): Directory containing paired-end FASTQ files.
        job (PipelineJob): The current pipeline job.

    Returns: 
        list[dict[str, str]]: Discovered sample records. 
    """
    # Normalizing and validating the input directory.
    input_dir = Path(input_dir)
    if not input_dir.exists() or not input_dir.is_dir():
        job.job_logger.error(f"Input directory does not exist: {input_dir}.")
        sys.exit(1)

    # Collecting every FASTQ in the directory (not recursive), sorted and de-duplicated.
    fastq_suffixes = [
        ".fastq.gz",
        ".fq.gz",
        ".fastq",
        ".fq"
    ]
    fastq_files: list[Path] = []
    for suffix in fastq_suffixes:
        fastq_files.extend(path for path in input_dir.glob(f"*{suffix}") if path.is_file())
    fastq_files = sorted(set(fastq_files))

    # Sorting each file into R1 or R2 by its read-pair token. The sample
    # ID is everything before that token, so both mates of one pair map
    # to the same ID. Illumina-style "_R1"/"_R2" (optionally followed by
    # "_001") and SRA-style "_1"/"_2" are recognised; the token must sit
    # right before the extension (or before "_001"), so a stray "_1"
    # inside a sample name (e.g. "FMNH_124727") can't be mistaken for it.
    r1_files: dict[str, Path] = {}
    r2_files: dict[str, Path] = {}
    unpaired: list[str] = []
    for fastq in fastq_files:
        match = FASTQ_PAIR_PATTERN.match(fastq.name)
        if match is None:
            unpaired.append(fastq.name)
            continue
        sample_id, mate = match.group("sample"), match.group("mate")
        target = r1_files if mate == "1" else r2_files
        if sample_id in target:
            job.job_logger.error(f"Two R{mate} files resolve to the same sample ID '{sample_id}': "
                                 f"{target[sample_id].name} and {fastq.name}. Provide a manifest to disambiguate.")
            sys.exit(1)
        target[sample_id] = fastq

    if unpaired:
        job.job_logger.warning(f"Ignoring {len(unpaired)} FASTQ file(s) with no recognisable read-pair token "
                               f"(_R1/_R2 or _1/_2 before the extension): {', '.join(unpaired)}.")

    # The finalized list of sample ids are the union of the r1 and r2 sample ids and validated that it is not empty.
    sample_ids = sorted(set(r1_files) | set(r2_files))
    if len(sample_ids) == 0:
        job.job_logger.error(f"No paired FASTQ files found in {input_dir}.")
        sys.exit(1)

    # Then we create the list of records - a dict containing each of the r1 and r2 files mapped to each other.
    # We validate each along the way, ensuring that each sample has a valid r1 and r2 file.
    records: list[dict[str, str]] = []
    for sample_id in sample_ids:
        if sample_id not in r1_files:
            job.job_logger.error(f"Missing R1 file for sample {sample_id}.")
            sys.exit(1)
        if sample_id not in r2_files:
            job.job_logger.error(f"Missing R2 file for sample {sample_id}.")
            sys.exit(1)
        records.append({
            "sample_id": sample_id,
            "r1": str(r1_files[sample_id]),
            "r2": str(r2_files[sample_id]),
            "source": "discovered_fastq"
        })

    job.job_logger.info(f"Discovered {len(records)} paired-end sample(s) in {input_dir}.")
    job.job_logger.warning("No manifest was provided, so discovered samples have no organism: MitoFish seed "
                           "selection falls back to GetOrganelle's default database for every sample, and "
                           "archival metadata (BioSample/GenBank) will be incomplete.")

    # Returning the list of records.
    return records

def write_discovered_manifest(
        manifest_path: Path,
        records: list[dict[str, str]],
        job: PipelineJob
) -> None:
    """
    Write discovered FASTQ records to a temporary runtime manifest.

    This manfiest is subsequently parsed and validated by the parse_sample_manifest() function 
    in manifest_parser.py so that discovered FASTQ inputs follow the same validation pathway as user-provided manifests.

    Args:
        manifest_path (Path): The path to the runtime manifest file.
        records (list[dict[str, str]]): A list of sample records.
        job (PipelineJob): The current pipeline job.
    """
    # Normalizing and creating the parent directory of the manifest path. 
    manifest_path = Path(manifest_path)
    manifest_path.parent.mkdir(parents = True, exist_ok = True)

    # List of columns to populate the manifest with. 
    columns = [
        "sample_id",
        "r1",
        "r2",
        "genus",
        "species",
        "source"
    ]

    # Creating a table and writing it to tsv format.
    table = pd.DataFrame(records, columns = columns)
    table.to_csv(manifest_path, sep = "\t", index = False)

    # Logging that tha manifest has been written. 
    job.job_logger.info(f"Discovered manifest written to {manifest_path}.")

def write_validated_manifest(
        manifest_path: Path,
        samples: list[Sample],
        job: PipelineJob
) -> None:
    """
    Write validated Sample objects to the runtime manifest.

    Required sample fields are taken directly from each Sample object.
    Additional nonblank metadata retained by the manifest parser is preserved.

    Args:
        manifest_path (Path): Path to the validated runtime manifest.
        samples (list[Sample]): Validated Sample objects.
        job (PipelineJob): The current pipeline job.
    """
    # Normalizing and creating the parent directory of the manifest path. 
    manifest_path = Path(manifest_path)
    manifest_path.parent.mkdir(parents = True, exist_ok = True)

    # Setting required columns for manifest. "organism" is always written
    # (from Sample.scientific_name, which also derives it from genus/species
    # when no explicit organism column was given, and is blank for a
    # discovered sample with no manifest at all) because the Snakemake rule
    # layer reads SAMPLE_TABLE["organism"] directly for seed selection and
    # rescue reseeding -- a runtime manifest without the column would raise
    # a KeyError there regardless of how the samples were provided.
    required_columns = ["sample_id", "r1", "r2", "organism"]

    # Creating an empty list for the metadata columns.
    metadata_columns: list[str] = []

    # Iterating through the list of samples, appending required columns and nonblank metadata columns to the list.
    for sample in samples:
        for column in sample.metadata:
            if column not in required_columns and column not in metadata_columns:
                metadata_columns.append(column)

    # Now creating a complete list of columns.
    columns = (required_columns + metadata_columns)

    # Creating an empty records list.
    records: list[dict[str, str]] = []

    # Iterating through the list of samples, appending the required columns to the list using the Sample data structure from sample.py. 
    for sample in samples:
        record = dict(sample.metadata)
        record["sample_id"] = sample.sample_id
        record["r1"] = str(sample.r1)
        record["r2"] = str(sample.r2)
        record["organism"] = sample.scientific_name or ""
        records.append(record)

    # Writing the records to a pandas dataframe and then exporting to tsv. 
    table = pd.DataFrame(records, columns = columns)
    table.to_csv(manifest_path, sep = "\t", index = False)
    
    # Logging that the manifest has been written.
    job.job_logger.info(f"Validated manifest written to {manifest_path}.")

def prepare_runtime_manifest(config: dict, job: PipelineJob) -> tuple[Path, list[Sample]]:
    """
    Prepares and validates a job-specific runtime manifest.

    If a user-provided manifest is configured, it is parsed directly by manifest parser in manifest_parser.py.
    If no manifest is configured, paired FASTQ files are discovered from the input directory, written to a temporary
    runtime manifest, and then passed through the manifest parser. 

    Both pathways therefore produce validated Sample objects using the same parser and validation logic.

    Args:
        config (dict): The pipeline configuration.
        job (PipelineJob): The current pipeline job.

    Returns:
        tuple[Path, list[Sample]]: A tuple containing the path to the runtime manifest and a list of validated Sample objects.
    
    """
    # Creating path to validated runtime manifest.
    runtime_manifest = job.job_dir / "runtime_manifest.tsv"
    job_configuration = _job_configuration(config)
    manifest_value = job_configuration.get("manifest")

    # Obtaining and validating the input directory if provided.
    input_dir_value = job_configuration.get("input_directory")
    input_directory = Path(input_dir_value).expanduser().resolve() if input_dir_value is not None else None

    # If a manifest is provided, this pathway is followed.
    if manifest_value is not None:
        manifest_path = Path(manifest_value)
        samples = parse_sample_manifest(manifest_path = manifest_path,
                                        logger = job.job_logger,
                                        input_directory = input_directory)
        write_validated_manifest(manifest_path = runtime_manifest,
                                 samples = samples,
                                 job = job)
        job.job_logger.info(f"Prepared validated runtime manifest from configured manifest: {runtime_manifest}.")
        return runtime_manifest, samples

    # If no manifest is provided, this pathway is followed.
    if input_directory is None:
        job.job_logger.error("Neither sample manifest nor input directory is provided. Aborting.")
        sys.exit(1)
    records = discover_fastq_pairs(input_dir = input_directory, job = job)
    write_discovered_manifest(manifest_path = runtime_manifest,
                              records = records,
                              job = job)
    samples = parse_sample_manifest(manifest_path = runtime_manifest, logger = job.job_logger)

    # Rewriting the manifest from validated Sample objects so that the manifest always contains validated Sample values produced by manifest parser.
    write_validated_manifest(manifest_path = runtime_manifest,
                             samples = samples,
                             job = job)
    job.job_logger.info(f"Prepared validated runtime manifest from discovered FASTQ files: {runtime_manifest}.")
    return runtime_manifest, samples

def write_runtime_config(
        config: dict,
        job: PipelineJob,
        runtime_manifest: Path
) -> Path:
    """
    Write a runtime configuration file for Snakemake.

    Args:
        config (dict): Original production configuration.
        job (PipelineJob): The current pipeline job.
        runtime_manifest (Path): The path to the runtime manifest.

    Returns:
        Path: The path to the runtime configuration file.
    """
    # Obtaining job-related fields from the config file. 
    runtime_config = dict(config)
    runtime_config["job_id"] = job.job_id
    runtime_config["manifest"] = str(runtime_manifest)
    runtime_config["runtime_manifest"] = str(runtime_manifest)
    runtime_config["output_root"] = str(job.parent_dir)
    runtime_config_path = job.job_dir / "runtime_config.yaml"

    # Writing the runtime config file.
    with runtime_config_path.open("w", encoding = "utf-8") as handle:
        yaml.safe_dump(
            runtime_config,
            handle,
            sort_keys = False
        )

    # Logging and returning the path to the runtime config file.
    job.job_logger.info(f"Runtime config written to {runtime_config_path}.")
    return runtime_config_path

def _extract_failure_reason(output: str, return_code: int) -> str:
    """
    Picks a short, specific description of why Snakemake failed out of
    its captured output. Snakemake names the failing rule on a
    line starting with "Error in rule <name>:"; that line (plus a
    couple that usually follow with the exit code/log path) is the
    most useful single thing to surface. Falls back to a generic,
    honest message when no such line is found (e.g. a config error
    before any rule ran) rather than guessing.

    Args:
        output (str): The full captured stdout+stderr from the Snakemake run.
        return_code (int): Snakemake's own process exit code.

    Returns:
        str: A short, specific failure reason.
    """
    lines = output.splitlines()
    for i, line in enumerate(lines):
        if line.strip().startswith("Error in rule"):
            return " ".join(l.strip() for l in lines[i:i + 4] if l.strip())
    return f"Snakemake exited with return code {return_code}; no specific rule-level error line was found in its output."


def run_snakemake(
        runtime_config: Path,
        cores: int,
        use_conda: bool,
        logger: logging.Logger,
        snakefile: Path = Path("ctrl/Snakefile"),
        dry_run: bool = False,
        printshellcmds: bool = False,
        output_log_path: Path | None = None,
        rerun_triggers: list[str] | None = None,
        concurrent_slots: int | None = None,
        global_log_path: Path | None = None
) -> tuple[subprocess.CompletedProcess, str]:
    """
    Run Snakemake using the runtime configuration.

    Args:
        runtime_config (Path): The path to the runtime configuration file.
        cores (int): The number of cores to use.
        use_conda (bool): Whether to use conda.
        logger (logging.Logger): The logger to use.
        snakefile (Path, optional): The path to the Snakefile. Defaults to "ctrl/Snakefile".
        dry_run (bool, optional): Whether to run Snakemake in dry-run mode. Defaults to False.
        printshellcmds (bool, optional): Whether to print shell commands. Defaults to False.
        rerun_triggers (list[str] | None): Passed through to Snakemake's own --rerun-triggers. None leaves Snakemake's default (every trigger) in effect.
        output_log_path (Path | None): Where to durably save Snakemake's full captured stdout/stderr, for later inspection after a crash. When None, output is still captured and returned but not saved to its own file. Written ONCE, after the whole run finishes -- not useful for watching a live run.
        concurrent_slots (int | None): The real ceiling on concurrent GetOrganelle-class jobs (see utils.system.resources.recommended_concurrent_slots) -- sets the `getorganelle_slots` resource pool. None falls back to this machine's own live-detected recommendation, computed fresh right here rather than ever hardcoded to one specific box.
        global_log_path (Path | None): The same global log file every exec script's own logger already writes to (see logger_factory.make_logger) -- when given, Snakemake's own scheduler output (job dispatch, "N of M steps done", rule reasons) is appended here LIVE, line by line, as it's produced, not just captured for the end-of-run output_log_path. Lets a single `tail -f` on the global log show both per-sample application activity and scheduler-level dispatch together. None skips this entirely (every prior caller that doesn't pass it keeps working exactly as before).

    Returns:
        tuple[subprocess.CompletedProcess, str]: The result of the Snakemake run, and its full captured stdout+stderr (Snakemake's own output is still streamed live to the console as it runs -- capturing it doesn't silence it).
    """
    if concurrent_slots is None:
        concurrent_slots = recommended_concurrent_slots()
    # Resolved against the running interpreter's own venv (sys.executable's
    # directory) rather than a bare "snakemake" -- subprocess.Popen does not
    # activate a venv, so a bare name only resolves if the venv's bin/ is
    # already first on $PATH, which invoking this launcher via its own
    # venv's python directly does not guarantee.
    snakemake_binary = str(Path(sys.executable).parent / "snakemake")

    # This is the command syntax required for snakemake. We explicitly set entrez to 1 to prevent
    # entrez from timing out/overloading when the task is parallelized.
    # --conda-not-block-search-path-envvars is required because Snakemake
    # otherwise sanitizes $PATH for any rule carrying a conda: directive
    # even without --use-conda, breaking resolution of externally-prepared
    # tool binaries (minimap2/samtools/blastn/GetOrganelle/mitos2).
    command = [
        snakemake_binary,
        "-s",
        str(snakefile),
        "--configfile",
        str(runtime_config),
        "--cores",
        str(cores),
        "--conda-not-block-search-path-envvars",
        # Conda activation adds real per-job overhead once --use-conda
        # is actually in effect; the default 5s latency-wait is too
        # tight for a rule's declared outputs to always be visible by
        # the time Snakemake checks for them.
        "--latency-wait",
        "30",
        # A real interruption (this process killed, the machine
        # rebooted, etc.) leaves Snakemake's own job metadata unable
        # to confirm the in-flight job's output is trustworthy -- it
        # refuses to proceed at all (a real IncompleteFilesException)
        # until told explicitly it's safe to redo that one job. Since
        # this launcher's whole purpose is unattended, resumable runs,
        # that confirmation has to be automatic, not a manual step
        # left for whoever restarts it.
        "--rerun-incomplete",
        "--resources",
        "entrez=1",
        # Caps concurrent GetOrganelle-class jobs (primary assembly,
        # rescue reassembly -- see their own rules' matching
        # `resources: getorganelle_slots = 1`) regardless of how many
        # idle threads `--cores` would otherwise allow. Confirmed live,
        # 2026-11-06: resuming with a large backlog let a burst of these
        # jobs launch within seconds of each other, and the kernel
        # logged real memory pressure in that window while many
        # reassembly attempts died before writing any log at all. Was
        # hardcoded to 8 (tuned for one specific machine's real core/RAM
        # mix); now computed live from THIS machine's own real CPU count
        # and total RAM (see recommended_concurrent_slots), so the same
        # launcher sizes itself correctly on a different box too, rather
        # than silently carrying over a number tuned for Dr. Stuart's
        # machine specifically.
        f"getorganelle_slots={concurrent_slots}"
    ]
    if use_conda: command.append("--use-conda")
    if dry_run: command.append("--dry-run")
    if printshellcmds: command.append("--printshellcmds")
    if rerun_triggers: command.extend(["--rerun-triggers", *rerun_triggers])

    # Logging and running the command. Output is streamed live to this
    # process's own stdout (so watching it run interactively still
    # works exactly as before) AND captured line-by-line into a buffer,
    # so a failure can be reported with more than just a return code.
    # Output is captured line-by-line into a buffer in addition to
    # streaming live, so it survives even if the terminal scrolled
    # past it or the session wasn't attended.
    logger.info(f"Running Snakemake with command: {' '.join(command)}")
    captured_lines: list[str] = []

    # Opened once, before the subprocess starts, and appended to for the
    # entire run -- every exec script's own make_logger() FileHandler
    # already appends to this same path concurrently from its own
    # subprocess; a line-buffered append-mode handle here is safe
    # alongside those for the same reason multiple exec scripts already
    # write to it safely today (each individual line write is small
    # enough that the OS's own O_APPEND write stays atomic).
    global_log_handle = None
    if global_log_path is not None:
        global_log_path.parent.mkdir(parents = True, exist_ok = True)
        global_log_handle = global_log_path.open("a", encoding = "utf-8")

    try:
        process = subprocess.Popen(command, cwd = Path.cwd(), stdout = subprocess.PIPE,
                                   stderr = subprocess.STDOUT, text = True, bufsize = 1)
        assert process.stdout is not None
        for line in process.stdout:
            print(line, end = "")
            captured_lines.append(line)
            if global_log_handle is not None:
                global_log_handle.write(f"[SNAKEMAKE] {line}")
                global_log_handle.flush()
        process.wait()
    finally:
        if global_log_handle is not None:
            global_log_handle.close()

    output = "".join(captured_lines)

    result = subprocess.CompletedProcess(args = command, returncode = process.returncode)

    if output_log_path is not None:
        try:
            output_log_path.parent.mkdir(parents = True, exist_ok = True)
            output_log_path.write_text(output, encoding = "utf-8")
        except OSError as error:
            logger.warning(f"Could not save Snakemake's captured output to {output_log_path}: {error}.")

    if result.returncode == 0: logger.info("Snakemake run completed successfully.")
    else: logger.error(f"Snakemake execution failed with return code {result.returncode}")

    return result, output

def _rerender_run_report_with_final_status(config: dict, job: PipelineJob, runtime_manifest: Path,
                                            samples: list[Sample]) -> None:
    """
    Re-renders just the run-level report (not the whole Snakemake DAG)
    now that the job's own terminal status has actually been written
    to current_metadata.json. The Snakemake rule that first renders
    this report runs mid-DAG, before this launcher's own
    mark_completed()/mark_failed() call below -- structurally, it can
    never see the job's own real final status, since that status is
    only known once the whole Snakemake subprocess has returned. Its
    first pass honestly labels this as "Running"/"Elapsed so far"
    rather than guessing (see Development Status, 2026-10-28); this
    second, lightweight pass -- the same real exec script, invoked
    again -- corrects it to the real terminal status without
    restructuring the DAG itself. A no-op when the reporting stage
    never ran in the first place (no run_report.done).

    Args:
        config (dict): The pipeline configuration.
        job (PipelineJob): The current pipeline job, already marked completed or failed.
        runtime_manifest (Path): The runtime manifest TSV used for this run.
        samples (list[Sample]): The run's own validated samples.
    """
    report_dir = job.job_dir / "13-report_pdf"
    if not (report_dir / "run_report.done").is_file():
        return

    pandoc_config = config.get("stage_options", {}).get("report_export", {}).get("pandoc", {})
    phylogeny_model = config.get("stage_options", {}).get("phylogeny", {}).get("iqtree", {}).get("model")

    exec_script = Path.cwd() / "src" / "mitopipeline" / "exec" / "13-report_pdf" / "run_run_report.py"
    command = [
        sys.executable, str(exec_script),
        "--job-id", job.job_id,
        *[arg for sample in samples for arg in ("--sample-id", sample.sample_id)],
        "--runtime-manifest", str(runtime_manifest),
        "--job-dir", str(job.job_dir),
        "--pandoc-bin", str(pandoc_config.get("bin", "pandoc")),
        "--pdf-engine", str(pandoc_config.get("engine", "typst")),
        "--output-dir", str(report_dir),
        "--log-file", str(job.job_dir / "logs" / "13-report_pdf" / "run_report.log"),
        "--global-log-file", str(job.job_dir / "logs" / "00-job" / f"{job.job_id}.global.log"),
    ]
    if phylogeny_model:
        command += ["--phylogeny-model", str(phylogeny_model)]

    job.job_logger.info("Re-rendering the run-level report now that the job's final status is known.")
    result = subprocess.run(command, cwd = Path.cwd(), capture_output = True, text = True)
    if result.returncode != 0:
        job.job_logger.warning(
            f"Final-status run-report re-render exited {result.returncode}; the report keeps its "
            f"original mid-run status. stderr: {result.stderr[-2000:]}"
        )

def _systemd_user_scope_available(systemd_run_binary: str) -> bool:
    """
    Whether `systemd-run --user --scope` can actually launch a real
    scope on this machine right now -- false on a machine with no
    reachable user systemd/D-Bus session even if the systemd-run
    binary itself exists. Checked with a real, trivial scope rather
    than assumed from the binary's mere presence, so a genuine
    environment failure here is never confused with the later wrapped
    pipeline's own legitimate nonzero exit.

    Args:
        systemd_run_binary (str): Path to the systemd-run binary (see shutil.which).

    Returns:
        bool: True if a real scope was launched and exited cleanly.
    """
    try:
        result = subprocess.run(
            [systemd_run_binary, "--user", "--scope", "--quiet", "--", "true"],
            capture_output = True, text = True, timeout = 15
        )
    except (OSError, subprocess.TimeoutExpired):
        return False
    return result.returncode == 0

def _reexec_under_memory_cgroup() -> int | None:
    """
    Re-launches this exact invocation (same interpreter, same argv)
    inside a memory-capped `systemd --user --scope` cgroup sized by
    cgroup_memory_max_mb, so a plain
    `python pipeline_launcher.py --config ...` gets the same OOM
    containment that previously required a separate, manual
    systemd-run wrapper. Fails open -- returns None (meaning: the
    caller should proceed directly, uncontained, no OOM monitor
    started) rather than guessing when systemd-run isn't installed, no
    user systemd session is reachable, or real total RAM is
    undetectable (e.g. local macOS development) -- matching every
    other function in utils.system.resources.

    Returns:
        int | None: The wrapped run's own real exit code if re-launched, or None if containment isn't available on this machine.
    """
    ceiling_mb = cgroup_memory_max_mb()
    systemd_run_binary = shutil.which("systemd-run")
    if ceiling_mb is None or systemd_run_binary is None or not _systemd_user_scope_available(systemd_run_binary):
        logging.warning(
            "Running without memory-cgroup containment on this machine (systemd-run unavailable, no "
            "reachable user systemd session, or real total RAM undetectable) -- the proactive OOM "
            "monitor will not be started either."
        )
        return None

    command = [
        systemd_run_binary, "--user", "--scope", "--quiet", "--collect",
        "-p", f"MemoryMax={ceiling_mb}M",
        "--", sys.executable, str(Path(__file__).resolve()), *sys.argv[1:]
    ]
    logging.info(f"Re-launching inside a memory-capped cgroup (MemoryMax={ceiling_mb}M).")
    result = subprocess.run(command, cwd = Path.cwd(), env = {**os.environ, _CGROUP_WRAPPED_ENV: "1"})
    return result.returncode

def _start_oom_monitor(cgroup_path: Path, job: PipelineJob) -> subprocess.Popen | None:
    """
    Starts launcher/oom_monitor.py as a standalone background process
    watching this job's own real cgroup, scoped to this one launch --
    see oom_monitor.py's own module docstring. Returns None (logging
    why) rather than starting anything when the real proactive
    threshold can't be computed, so a caller never watches a cgroup
    with no real ceiling to react to.

    Args:
        cgroup_path (Path): This process's own real cgroup directory (see utils.system.resources.current_cgroup_path), only meaningful once already placed inside a fresh systemd-run --user --scope by _reexec_under_memory_cgroup.
        job (PipelineJob): The current pipeline job -- supplies the job root (oom_monitor's own state file lives under its sanctioned mlr_<job_id> scratch folder) and the logger.

    Returns:
        subprocess.Popen | None: The running monitor process, or None if it wasn't started.
    """
    threshold_mb = oom_monitor_threshold_mb()
    if threshold_mb is None:
        job.job_logger.warning("OOM monitor not started: real total RAM is undetectable on this machine.")
        return None

    monitor_script = Path(__file__).resolve().parent / "oom_monitor.py"
    log_file = job.job_dir / "logs" / "00-job" / f"{job.job_id}.oom_monitor.log"
    command = [
        sys.executable, str(monitor_script),
        "--cgroup-path", str(cgroup_path),
        "--threshold-mb", str(threshold_mb),
        "--job-root", str(job.job_dir),
        "--log-file", str(log_file),
    ]
    job.job_logger.info(f"Starting OOM monitor: cgroup={cgroup_path}, threshold={threshold_mb}MB.")
    return subprocess.Popen(command, cwd = Path.cwd())

def _stop_oom_monitor(process: subprocess.Popen, logger: logging.Logger) -> None:
    """
    Stops a monitor process started by _start_oom_monitor -- it runs
    forever on its own (see oom_monitor.run_monitor's own polling
    loop), so this launcher must be the one to end its lifetime once
    this one pipeline run is actually done, success or failure alike.

    Args:
        process (subprocess.Popen): The monitor process to stop.
        logger (logging.Logger): Where to log the shutdown.

    Returns:
        None
    """
    if process.poll() is not None:
        return
    logger.info("Stopping OOM monitor.")
    process.terminate()
    try:
        process.wait(timeout = 10)
    except subprocess.TimeoutExpired:
        process.kill()
        process.wait()

def resolve_concurrent_slots(config: dict) -> int:
    """
    This run's own real concurrency ceiling for GetOrganelle-class
    jobs, resolved from the pipeline's only two user-facing resource
    knobs -- `resources.cores`/`resources.ram_gb`. Both null (the
    default) means fully live-detected, never hardcoded to one box.
    `getorganelle_slots` itself is never a knob the user sets
    directly -- a third, redundant number that could silently
    disagree with cores/ram_gb -- it's purely derived from these two,
    once, here.

    Args:
        config (dict): The run's own config (job or runtime) -- reads config["resources"]["cores"]/["ram_gb"].

    Returns:
        int: The real recommended concurrent-slots ceiling (see utils.system.resources.recommended_concurrent_slots).
    """
    resources_config = config.get("resources") or {}
    configured_cores = resources_config.get("cores")
    configured_ram_gb = resources_config.get("ram_gb")
    configured_ram_mb = int(configured_ram_gb * 1024) if configured_ram_gb is not None else None
    return recommended_concurrent_slots(cpu_count_override = configured_cores,
                                        total_ram_mb_override = configured_ram_mb)

def launch_pipeline(
        config_path: Path,
        snakefile: Path = Path("ctrl/Snakefile"),
        cores: int = 1,
        use_conda: bool = True,
        dry_run: bool = False,
        printshellcmds: bool = False,
        rerun_triggers: list[str] | None = None,
        cgroup_path: Path | None = None
):
    """
    Launch a single pipeline run.

    Args:
        config_path (Path): The path to the config file.
        snakefile (Path, optional): The path to the Snakefile. Defaults to "ctrl/Snakefile".
        cores (int, optional): The number of cores to use for the job. Defaults to 1.
        use_conda (bool, optional): Whether to use conda. Defaults to True.
        dry_run (bool, optional): Whether to run the pipeline in dry-run mode. Defaults to False.
        printshellcmds (bool, optional): Whether to print shell commands. Defaults to False.
        rerun_triggers (list[str] | None, optional): Passed through to Snakemake's own --rerun-triggers. Defaults to None (Snakemake's default).
        cgroup_path (Path | None, optional): This process's own real cgroup directory, once already placed inside a memory-capped systemd scope by _reexec_under_memory_cgroup -- starts the proactive OOM monitor against it for the duration of this run. None (the default) skips the monitor entirely, e.g. when cgroup containment isn't available on this machine.

    Returns:
        None
    """
    config = load_config(config_path)
    job_configuration = _job_configuration(config)
    output_root = Path(job_configuration["output_root"])
    output_root.mkdir(parents = True, exist_ok = True)

    job_id = _configured_job_id(job_configuration.get("job_id"))
    job = PipelineJob(parent_dir = output_root, job_id = job_id)

    oom_monitor_process = _start_oom_monitor(cgroup_path, job) if cgroup_path is not None else None
    try:
        # Running the pipeline - the try block
        # creates the runtime manifest, obtains the samples, and then runs snakemake.
        # If any error occurs, the job is marked as failed.
        try:
            if job_id is not None:
                job.job_logger.info(f"Resuming configured pipeline job: {job_id}.")
            runtime_manifest, samples = prepare_runtime_manifest(config = config, job = job)
            job.set_samples([sample.sample_id for sample in samples])
            job.mark_created()

            # Resolving the real concurrency ceiling for GetOrganelle-class
            # jobs ONCE here, so the exact same number ends up on both the
            # `--resources getorganelle_slots=N` launch flag (run_snakemake)
            # AND in the written runtime_config.yaml (read by each rule as
            # `--concurrent-slots`, to size that job's own live memory
            # budget) -- never two independently-derived numbers that could
            # drift apart.
            concurrent_slots = resolve_concurrent_slots(config)
            resources_config = config.get("resources") or {}
            runtime_config_values = dict(config)
            runtime_config_values["resources"] = {**resources_config, "getorganelle_slots": concurrent_slots}

            runtime_config = write_runtime_config(config = runtime_config_values,
                                                  job = job,
                                                  runtime_manifest = runtime_manifest)
            result, snakemake_output = run_snakemake(
                snakefile = snakefile,
                runtime_config = runtime_config,
                cores = cores,
                use_conda = use_conda,
                logger = job.job_logger,
                dry_run = dry_run,
                printshellcmds = printshellcmds,
                rerun_triggers = rerun_triggers,
                output_log_path = job.job_dir / "logs" / "00-job" / f"{job.job_id}.snakemake_output.log",
                concurrent_slots = concurrent_slots,
                # The exact same path every exec script's own make_logger()
                # already appends to (see ctrl/Snakefile's GLOBAL_LOG_FILE) --
                # so a single tail on this one file shows scheduler dispatch
                # and per-sample application activity together, live.
                global_log_path = job.job_dir / "logs" / "00-job" / f"{job.job_id}.global.log"
            )
            if result.returncode == 0:
                job.mark_completed()
            else:
                # The last ~40 lines are usually enough to see the actual
                # error without dumping an entire multi-hour run's output
                # into the job's metadata file; the full output is still on
                # disk at output_log_path for anyone who needs more.
                log_tail = "\n".join(snakemake_output.splitlines()[-40:])
                job.mark_failed(
                    failure_reason = _extract_failure_reason(snakemake_output, result.returncode),
                    log_tail = log_tail
                )
            _rerender_run_report_with_final_status(config, job, runtime_manifest, samples)
            return result.returncode
        except Exception as error:
            job.job_logger.exception(f"Pipeline launch failed: {error}.")
            job.mark_failed(failure_reason = f"{type(error).__name__}: {error}")
            return 1
    finally:
        if oom_monitor_process is not None:
            _stop_oom_monitor(oom_monitor_process, job.job_logger)

def _configured_job_id(value: object) -> str | None:
    """
    Normalizes a job ID from the config.

    Args:
        value (object): The job ID to normalize.

    Returns:
        str | None: The normalized job ID, or None if the value is None.
    """
    # 1. If the value is None, return None.
    # 2. Convert the value to a string and strip whitespace.
    # 3. If the string is empty or "null", return None.
    # 4. Otherwise, return the string.
    if value is None: return None
    text = str(value).strip()
    if text == "" or text.lower() == "null": return None
    return text

def main() -> int:
    """
    Runs the pipeline launcher from the command line.

    On a real systemd/cgroup-v2 Linux machine, the first (unwrapped)
    invocation transparently re-launches itself inside a memory-capped
    `systemd --user --scope` cgroup and starts the proactive OOM
    monitor against that same cgroup -- see
    _reexec_under_memory_cgroup and _start_oom_monitor. A plain
    `python pipeline_launcher.py --config ...` is therefore enough on
    its own; no separate accessory script or manual systemd-run
    wrapper is required. Fails open to a direct, uncontained run (no
    monitor) on any machine where that containment genuinely isn't
    available (e.g. local macOS development).

    Returns:
        int: Pipeline process return code.
    """
    already_wrapped = os.environ.get(_CGROUP_WRAPPED_ENV) == "1"
    if not already_wrapped:
        wrapped_returncode = _reexec_under_memory_cgroup()
        if wrapped_returncode is not None:
            return wrapped_returncode

    args = parse_args()
    return launch_pipeline(
        config_path = Path(args.config),
        snakefile = Path(args.snakefile),
        cores = args.cores,
        use_conda = not args.no_conda,
        dry_run = args.dry_run,
        printshellcmds = args.printshellcmds,
        rerun_triggers = args.rerun_triggers,
        cgroup_path = current_cgroup_path() if already_wrapped else None
    )

# Running the main command.
if __name__ == "__main__":
    raise SystemExit(main())