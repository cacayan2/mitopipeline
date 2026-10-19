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
from pathlib import Path
import subprocess
import pandas as pd
import yaml
import sys
from mitopipeline.utils.manifest.manifest_parser import parse_sample_manifest
from mitopipeline.models.pipeline_job.pipeline_job import PipelineJob
from mitopipeline.models.sample.sample import Sample

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
    # The minimum command to run the pipeline is:
    ## python pipeline_launcher.py --config <path/to/config.yaml>
    # Recommended:
    ## python pipeline_launcher.py --config <path/to/config,yaml> --cores <number>
    parser.add_argument("--config", type = Path, required = True, help = "The path to the config file.")
    parser.add_argument("--snakefile", type = Path, default = "ctrl/Snakefile", help = "The path to the snakefile (optional).")
    parser.add_argument("--cores", type = int, default = 1, help = "The number of cores to use for the job (optional, but 1 is default so recommend changing this).")
    parser.add_argument("--no-conda", action = "store_false", help = "The flag to disable the conda environment (optional, but recommend not messing with this).")
    parser.add_argument("--dry-run", action = "store_true", help = "The flag to run the pipeline in dry-run mode (optional).")
    parser.add_argument("--printshellcmds", action = "store_true", help = "The flag to print shell commands (optional - useful for debugging).")

    # Returning the parsed arguments.
    return parser.parse_args()

def load_config(config_path: Path, job: PipelineJob) -> dict:
    """
    Loads a YAML configuration file.

    Args:
        config_path (Path): The path to the YAML configuration file.
        job (PipelineJob): The pipeline job object - used for logging.

    Returns:
        dict: The loaded configuration as a dictionary.
    """
    # Normalizing the config path and validating that it exists.
    config_path = Path(config_path)
    if not config_path.exists() or not config_path.is_file():
        job.job_logger.error(f"Config file does not exist: {config_path}.")
        sys.exit(1)

    # Loading the config file and validating that it is not empty.
    with open(config_path, "r", encoding = "utf-8") as handle: 
        config = yaml.safe_load(handle)
    if config is None:
        job.job_logger.error(f"Config file is empty: {config_path}.")
        sys.exit(1)

    # Returning the loaded config.
    return config

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

    # Creating a list of acceptable fastq file suffixes and initializing an empty dictionary for the bulk list of files.
    fastq_suffixes = [
        ".fastq.gz",
        ".fq.gz",
        ".fastq",
        ".fq"
    ]
    fastq_files = list[Path] = []

    # Iterating through the suffixes and searching through the input_dir for each one - appending to the list of files.
    # This is followed up by sorting and removing duplicates. 
    for suffix in fastq_suffixes:
        fastq_files.extend(input_dir.glob(f"*{suffix}"))
    fastq_files = sorted(set(fastq_files))

    # Now creating individual dictionaries for r1 and r2 files. The key is the sample ID, and the value is the path to the file.
    # We then iterate through each of the fastq files and determine if it is an R1 or R2 file, and extract the sample ID from each.  
    r1_files = dict[str, Path] = {}
    r2_files = dict[str, Path] = {}
    for fastq in fastq_files:
        name = fastq.name
        if "_R1" in name:
            sample_id = name.split("_R1")[0]
            r1_files[sample_id] = fastq
        elif "_R2" in name:
            sample_id = name.split("_R2")[0]
            r2_files[sample_id] = fastq
        elif "_1" in name:
            sample_id = name.split("_1")[0]
            r1_files[sample_id] = fastq
        elif "_2" in name:
            sample_id = name.split("_2")[0]
            r2_files[sample_id] = fastq

    # The finalized list of sample ids are the union of the r1 and r2 sample ids and validated that it is not empty.
    sample_ids = sorted(set(r1_files) | set(r2_files))
    if len(sample_ids) == 0:
        job.job_logger.error(f"No fastq files found in {input_dir}.")
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
            "genus": "",
            "species": "",
            "source": "discovered_fastq"
        })

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

    # Setting required columns for manifest.
    required_columns = ["sample_id", "r1", "r2"]

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
    manifest_value = config.get("manifest")

    # Obtaining and validating the input directory if provided.
    input_dir_value = config.get("input_directory")
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
        job.logger.error("Neither sample manifest nor input directory is provided. Aborting.")
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
    its captured output (Emil's request, 2026-10-19: a crash should say
    more than just "it failed"). Snakemake names the failing rule on a
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
        output_log_path: Path | None = None
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
        output_log_path (Path | None): Where to durably save Snakemake's full captured stdout/stderr, for later inspection after a crash. When None, output is still captured and returned but not saved to its own file.

    Returns:
        tuple[subprocess.CompletedProcess, str]: The result of the Snakemake run, and its full captured stdout+stderr (Snakemake's own output is still streamed live to the console as it runs -- capturing it doesn't silence it).
    """
    # This is the command syntax required for snakemake. We explicitly set entrez to 1 to prevent
    # entrez from timing out/overloading when the task is parallelized.
    command = [
        "snakemake",
        "-s",
        str(snakefile),
        "--configfile",
        str(runtime_config),
        "--cores",
        str(cores),
        "--resources",
        "entrez = 1"
    ]
    if use_conda: command.append("--use-conda")
    if dry_run: command.append("--dry-run")
    if printshellcmds: command.append("--printshellcmds")

    # Logging and running the command. Output is streamed live to this
    # process's own stdout (so watching it run interactively still
    # works exactly as before) AND captured line-by-line into a buffer,
    # so a failure can be reported with more than just a return code
    # (Emil's request, 2026-10-19) -- previously this output went
    # straight to the terminal and nowhere else, lost the moment the
    # terminal scrolled past it or the session wasn't attended.
    logger.info(f"Running Snakemake with command: {' '.join(command)}")
    captured_lines: list[str] = []
    process = subprocess.Popen(command, cwd = Path.cwd(), stdout = subprocess.PIPE,
                               stderr = subprocess.STDOUT, text = True, bufsize = 1)
    assert process.stdout is not None
    for line in process.stdout:
        print(line, end = "")
        captured_lines.append(line)
    process.wait()
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

def launch_pipeline(
        config_path: Path,
        snakefile: Path = Path("ctrl/Snakefile"),
        cores: int = 1,
        use_conda: bool = True,
        dry_run: bool = False,
        printshellcmds: bool = False
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

    Returns:
        None
    """
    config = load_config(config_path)
    output_root = Path(config["output_root"])
    output_root.mkdir(parents = True, exist_ok = True)

    job_id = _configured_job_id(config.get("job_id"))
    job = PipelineJob(parent_dir = output_root, job_id = job_id)

    # Running the pipeline - the try block 
    # creates the runtime manifest, obtains the samples, and then runs snakemake. 
    # If any error occurs, the job is marked as failed. 
    try:
        if job_id is not None:
            job.job_logger.info(f"Resuming configured pipeline job: {job_id}.")
        runtime_manifest, samples = prepare_runtime_manifest(config = config, job = job)
        job.set_samples([sample.sample_id for sample in samples])
        job.mark_created()
        runtime_config = write_runtime_config(config = config,
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
            output_log_path = job.job_dir / "logs" / "00-job" / f"{job.job_id}.snakemake_output.log"
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
        return result.returncode
    except Exception as error:
        job.job_logger.exception(f"Pipeline launch failed: {error}.")
        job.mark_failed(failure_reason = f"{type(error).__name__}: {error}")
        return 1

def main() -> int:
    """
    Runs the pipeline launcher from the command line.

    Returns:
        int: Pipeline process return code.
    """
    args = parse_args()
    return launch_pipeline(
        config_path = Path(args.config),
        snakefile = Path(args.snakefile),
        cores = args.cores,
        use_conda = not args.no_conda,
        dry_run = args.dry_run,
        printshellcmds = args.printshellcmds
    )

# Running the main command. 
if __name__ == "__main__":
    raise SystemExit(main())

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