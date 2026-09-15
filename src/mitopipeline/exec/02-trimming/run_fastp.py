"""
run_fastp.py

Execution layer connecting fastp to Snakemake.
"""

# Imports
import argparse
from logging import Logger
from pathlib import Path
from mitopipeline.api.fastp import FastpRunner
from mitopipeline.models.sample import Sample
from mitopipeline.logging.logger_factory import make_logger

def parse_args() -> argparse.Namespace: 
    """Parse command-line arguments.
    
    Returns:
        argparse.Namespace: Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description = "Process sample inputs and directory paths."
    )

    # Define arguments for files.
    parser.add_argument("--sample-id", type=str, required=True, help="The unique identifier for the sample.")
    parser.add_argument("--in1", type=Path, required=True, help="Path to the R1 input FASTQ.")
    parser.add_argument("--in2", type=Path, required=True, help="Path to the R2 input FASTQ.")
    parser.add_argument("--output-dir", type=Path, required=True, help="Directory for fastp outputs.")
    parser.add_argument("--working-dir", type=Path, default=".", help="Working directory for execution.")
    parser.add_argument("--threads", type=int, required=True, help="Number of threads.")
    parser.add_argument("--log-file", type = Path, required=True, help="Path to the log file.")
    parser.add_argument("--global-log-file", type = Path, required=True, help="Path to the global log file.")
    parser.add_argument(
    "--qualified_quality_phred",
    type=int,
    default=30,
    help="Minimum Phred quality score."
    )
    parser.add_argument(
        "--length_required",
        type=int,
        default=50,
        help="Minimum read length after filtering."
    )
    parser.add_argument(
        "--trim_front1",
        type=int,
        default=0,
        help="Trim fixed bases from the front of R1."
    )
    parser.add_argument(
        "--trim_front2",
        type=int,
        default=0,
        help="Trim fixed bases from the front of R2."
    )
    parser.add_argument(
        "--trim_tail1",
        type=int,
        default=0,
        help="Trim fixed bases from the end of R1."
    )
    parser.add_argument(
        "--trim_tail2",
        type=int,
        default=0,
        help="Trim fixed bases from the end of R2."
    )
    parser.add_argument(
        "--cut_window_size",
        type=int,
        default=4,
        help="Sliding window size."
    )
    parser.add_argument(
        "--cut_mean_quality",
        type=int,
        default=20,
        help="Minimum average quality within the sliding window."
    )
    parser.add_argument(
        "--n_base_limit",
        type=int,
        default=None,
        help="Maximum number of N bases allowed."
    )
    parser.add_argument(
        "--unqualified_percent_limit",
        type=int,
        default=None,
        help="Maximum percentage of low-quality bases."
    )
    parser.add_argument(
        "--average_qual",
        type=int,
        default=None,
        help="Minimum average read quality."
    )
    parser.add_argument(
        "--report_title",
        default=None,
        help="Title displayed in the HTML report."
    )

    parser.add_argument(
        "--adapter_sequence",
        default=None,
        help="Adapter sequence for read 1."
    )

    parser.add_argument(
        "--adapter_sequence_r2",
        default=None,
        help="Adapter sequence for read 2."
    )

    parser.add_argument(
        "--adapter_fasta",
        default=None,
        help="FASTA file containing adapter sequences."
    )
    parser.add_argument(
    "--detect_adapter_for_pe",
    action="store_true",
    help="Automatically detect adapters for paired-end reads."
    )
    parser.add_argument(
        "--cut_front",
        action="store_true",
        help="Enable front sliding-window trimming."
    )

    parser.add_argument(
        "--cut_tail",
        action="store_true",
        help="Enable tail sliding-window trimming."
    )

    parser.add_argument(
        "--cut_right",
        action="store_true",
        help="Enable right-end sliding-window trimming."
    )

    parser.add_argument(
        "--disable_quality_filtering",
        action="store_true",
        help="Disable quality filtering."
    )

    parser.add_argument(
        "--disable_length_filtering",
        action="store_true",
        help="Disable length filtering."
    )

    parser.add_argument(
        "--trim_poly_g",
        action="store_true",
        help="Enable poly-G trimming."
    )

    parser.add_argument(
        "--disable_trim_poly_g",
        action="store_true",
        help="Disable automatic poly-G trimming."
    )

    parser.add_argument(
        "--trim_poly_x",
        action="store_true",
        help="Enable poly-X trimming."
    )

    # Return the parsed arguments.
    return parser.parse_args()

def main() -> int:
    """
    Runs fastp for a single sequencing sample.

    Args:
        None

    Returns:
        int: The return code.
    """
    # Parsing the command line arguments.
    args = parse_args()

    # Creating a logger for this step.
    logger = make_logger(
        name = "fastp",
        log_file_path = args.log_file,
        global_log_file_path = args.global_log_file
    )

    # Creating the sample object.
    sample = Sample(
        sample_id = args.sample_id,
        r1 = Path(args.r1),
        r2 = Path(args.r2)
    )

    # Creating tool options dictionary.
    tool_options = {
        "qualified_quality_phred": args.qualified_quality_phred,
        "length_required": args.length_required,
        "trim_front1": args.trim_front1,
        "trim_tail1": args.trim_tail1,
        "trim_front2": args.trim_front2,
        "trim_tail2": args.trim_tail2,
        "cut_window_size": args.cut_window_size,
        "cut_mean_quality": args.cut_mean_quality,
        "n_base_limit": args.n_base_limit,
        "unqualified_percent_limit": args.unqualified_percent_limit,
        "average_qual": args.average_qual,
        "report_title": args.report_title,
        "adapter_sequence": args.adapter_sequence,
        "adapter_sequence_r2": args.adapter_sequence_r2,
        "adapter_fasta": args.adapter_fasta,
        "detect_adapter_for_pe": args.detect_adapter_for_pe,
        "cut_front": args.cut_front,
        "cut_tail": args.cut_tail,
        "cut_right": args.cut_right,
        "disable_quality_filtering": args.disable_quality_filtering,
        "disable_length_filtering": args.disable_length_filtering,
        "trim_poly_g": args.trim_poly_g,
        "disable_trim_poly_g": args.disable_trim_poly_g,
        "trim_poly_x": args.trim_poly_x,
    }

    # Creating the FastPRunner and obtaining context for logger.
    runner = FastpRunner(
        working_dir = Path(args.working_dir),
        output_dir = Path(args.output_dir),
        sample = sample,
        tool_options = tool_options,
        threads = args.threads,
        logger = logger,
        tool_name = "fastp"
    )
    context = runner._log_context()

    # Logging the start of the step.
    logger.info(f"{context} Initiating fastp on raw sequencing data.")

    # Executing fastp.
    try:
        result = runner.run()
        logger.info(f"{context} Successfully ran fastp on raw sequencing data.")
        return 0 if result.success else result.return_code or 1
    except Exception as e:
        logger.error(f"{context} Failed to run fastp on raw sequencing data: {e}")
        return 1

