"""
fastp.py

Fastp wrapper for command building and execution. 
"""

# Imports
from __future__ import annotations
from logging import Logger
from pathlib import Path
import os
import sys
from mitopipeline.models.base_tool.base_tool import BaseTool
from mitopipeline.models.sample.sample import Sample

class FastpRunner(BaseTool):
    """
    This class provides the functionality to run fastp and normalize I/O
    within a single sample directory.

    Attributes:
        working_dir (Path): The path to the directory of the input files.
        output_dir (Path): The directory to save the output files to.
        sample (Sample): The sample object that is being worked with.
        tool_options (dict): A dictionary of options to pass to the tool.
        r1_output_stem (str): The stem of the R1 file - this will be appended to with file extensions for the outputs.
        r2_output_stem (str): The stem of the R2 file - this will be appended to with file extensions for the outputs.
        logger (Logger): The logger to use for this of the pipeline.
        tool_name (str): The name of the tool (by default will be fastp).
        threads (int): The number of the threads for this process to use.
        output_r1 (Path): The path to the output R1 file.
        output_r2 (Path): The path to the output R2 file.
        output_html (Path): The path to the output HTML file.
        output_json (Path): The path to the output JSON file.
    """
    def __init__(self, 
                 working_dir: Path,
                 output_dir: Path,
                 sample: Sample,
                 tool_options: dict | None = None,
                 threads: int = 4,
                 logger: Logger | None = None,
                 tool_name: str = "fastp") -> None:
        """
        Instantiation method for the creation of the 
        FastpRunner class. Calls the BaseTool instantiation methods and
        sets a few other object variables. 
        
        Args:
            working_dir (Path): The path to the directory of the input files.
            output_dir (Path): The directory to save the output files to.
            sample (Sample): The sample object that is being worked with.
            tool_options (dict): A dictionary of options to pass to the tool.
            threads (int): The number of the threads for this process to use.
            logger (Logger): The logger to use for this of the pipeline.
            tool_name (str): The name of the tool (by default will be fastp).
        
        Returns: 
            None
        """ 
        # Using the parent class constructor and setting the output file names. 
        super().__init__(working_dir, output_dir, sample, tool_name, threads, logger, tool_options)
        self.output_r1 = Path(self.output_dir / f"{self.sample.name}_R1.fastq.gz")
        self.output_r2 = Path(self.output_dir / f"{self.sample.name}_R2.fastq.gz")
        self.output_html = Path(self.output_dir / f"{self.sample.name}_fastp.html")
        self.output_json = Path(self.output_dir / f"{self.sample.name}_fastp.json")
        self.tool_options = tool_options or {}

    def validate_inputs(self) -> None:
        """
        Validates the input files for fastp.
        This function appropriately logs the result, raises an error otherwise.

        Args:
            None

        Returns:
            None
        """
        # Adding logger context.
        context = super()._log_context()

        # Checking if input files exist.
        for path in (self.sample.r1, self.sample.r2):
            if not path.is_file():
                self.logger.error(f"{context} Input file does not exist: {path}.")
                sys.exit(1)
            if self.threads <= 0:
                self.logger.error(f"{context} Please provide a valid number of threads (passed = {self.threads}).")
                sys.exit(1)
            if self.threads > os.cpu_count():
                self.logger.error(f"{context} Number of threads is greater than number of cores. Please reduce the number of threads to {os.cpu_count()}.")
                sys.exit(1)

    def build_command(self) -> list[str]:
        """
        Constructs the command to run fastp.

        Args:
            None

        Returns:
            list[str]: The command to run fastp.
        """
        # Creating the output directory if it does not exist.
        self.output_dir.mkdir(parents = True, exist_ok = True)

        # Constructing the command.
        command = [
            "fastp",
            "--threads",
            "--in1",
            str(self.sample.r1),
            "--in2",
            str(self.sample.r2),
            "--out1",
            str(self.output_r1),
            "--out2",
            str(self.output_r2),
            "--html",
            str(self.output_html),
            "--json",
            str(self.output_json),
            "--thread",
            str(self.threads)
        ]

        # Returning the command with additional options.
        return self._add_options(command)

    def validate_outputs(self) -> None:
        """
        Validates the outputs of fastp.
        This function appropriately logs the result,
        raises an error otherwise.

        Args:
            None
        
        Returns:
            None
        """
        # Iterating through each of the expected outputs - if they do not exist or are empty,
        # then raise an error.
        for path in (
            self.output_r1,
            self.output_r2,
            self.output_html,
            self.output_json
        ):
            if not path.is_file() or path.stat().st_size == 0:
                self.logger.error(f"Missing or empty output: {path}")
                sys.exit(1)

    def _add_options(self, command: list[str]) -> list[str]:
        """
        Helper function which appends optional flags to an original command.

        Args:
            command (list[str]): The command to modify.

        Returns:
            list[str]: The modified command.
        """
        # A dictionary of options that require a value -
        # the key is the config name and the value is the 
        # true flag passed to fastp.
        value_options = {
            "qualified_quality_phred": "--qualified_quality_phred",
            "length_required": "--length_required",
            "trim_front1": "--trim_front1",
            "trim_tail1": "--trim_tail1",
            "trim_front2": "--trim_front2",
            "trim_tail2": "--trim_tail2",
            "cut_window_size": "--cut_window_size",
            "cut_mean_quality": "--cut_mean_quality",
            "n_base_limit": "--n_base_limit",
            "unqualified_percent_limit": "--unqualified_percent_limit",
            "average_qual": "--average_qual",
            "report_title": "--report_title",
            "adapter_sequence": "--adapter_sequence",
            "adapter_sequence_r2": "--adapter_sequence_r2",
            "adapter_fasta": "--adapter_fasta",
        }

        # Iterating through each of the value options and
        # appending values that are not None to the command.
        for key, flag in value_options.items():
            value = self.tool_options.get(key)
            if value is not None:
                command.extend([flag, str(value)])

        # A dictionary of options that require a boolean value -
        # the key is the config name and the value is the 
        # true flag passed to fastp.
        boolean_options = {
            "detect_adapter_for_pe": "--detect_adapter_for_pe",
            "cut_front": "--cut_front",
            "cut_tail": "--cut_tail",
            "cut_right": "--cut_right",
            "disable_quality_filtering": "--disable_quality_filtering",
            "disable_length_filtering": "--disable_length_filtering",
            "trim_poly_g": "--trim_poly_g",
            "disable_trim_poly_g": "--disable_trim_poly_g",
            "trim_poly_x": "--trim_poly_x",
        }

        # Iterating through each of the boolean options and
        # appending values that are True to the command.
        for key, flag in boolean_options.items():
            value = self.tool_options.get(key)
            if value is True:
                command.append(flag)

        # Returning the modified command.
        return command
            
