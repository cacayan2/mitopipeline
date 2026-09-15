"""
run_fastqc_raw.py

Execution layer for FastQC on raw sequencing data.
"""

# Imports
import argparse
from pathlib import Path
from mitopipeline.api.fastqc import FastQCRunner
from mitopipeline.logging.logger_factory import make_logger
from mitopipeline.models.sample import Sample

def parse_args() -> argparse.Namespace:
    """
    Parses the command line arguments for this script.

    Returns:
        argparse.Namespace: The parsed command line arguments.
    """
    # Creating the argument parser.
    parser = argparse.ArgumentParser(description = "Runs FastQC on raw sequencing data.")

    # Adding arguments to the parser.
    parser.add_argument("--sample-id", type = str, required = True, help = "The sample id.")
    parser.add_argument("--r1", type = str, required = True, help = "The path to the R2 fastq file.")
    parser.add_argument("--r2", type = str, required = True, help = "The path to the R2 fastq file.")
    parser.add_argument("--output-dir", type = Path, required = True, help = "The path to the output directory.")
    parser.add_argument("--working-dir", type = Path, required = True, help = "The path to the working directory.")
    parser.add_argument("--log-file", type = Path, required = True, help = "The path to the log file.")
    parser.add_argument("--global-log-file", type = Path, required = True, help = "The path to the global log file.")
    parser.add_argument("--threads", type = int, required = True, help = "The number of threads to use.")

    # Parsing the arguments.
    return parser.parse_args()

def main() -> int:
    """
    Runs FastQC on raw sequencing data.

    Returns:
        int: The return code for the process.
    """
    # Parsing the command line arguments.
    args = parse_args()

    # Creating a logger for this step.
    logger = make_logger(
        name = "qc_raw",
        log_file_path = args.log_file,
        global_log_file_path = args.global_log_file
    )

    # Creating the sample object.
    sample = Sample(
        sample_id = args.sample_id,
        r1 = Path(args.r1),
        r2 = Path(args.r2)
    )

    # Creating the FastQC runner and obtaining context for logger.
    runner = FastQCRunner(
        working_dir = Path(args.working_dir),
        output_dir = Path(args.output_dir),
        sample = sample,
        r1_output_stem = "R1",
        r2_output_stem = "R2",
        threads = args.threads,
        logger = logger,
        threads = args.threads
    )
    context = runner._log_context()

    # Logging the start of the step.
    logger.info(f"{context} Initiating FastQC on raw sequencing data.")

    # Creating the sample object.
    sample = Sample(
        sample_id = args.sample_id,
        r1 = Path(args.r1),
        r2 = Path(args.r2)
    )

    # Executing FastQC.
    try: 
        result = runner.run()
        logger.info(f"{context} FastQC execution on raw sample {sample.sample_id} complete.")
        return 0 if result.success else result.return_code or 1
    except Exception as e:
        logger.error(f"{context} FastQC execution on raw sample {sample.sample_id} failed: {e}.")
        return 1

if __name__ == "__main__":
    raise SystemExit(main())
