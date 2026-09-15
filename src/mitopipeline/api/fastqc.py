"""
fastqc.py

FastQC wrapper which supports parallelization using a single sample directory.
"""

# Imports
from __future__ import annotations
from logging import Logger
from pathlib import Path
import shutil
import os
import sys
from mitopipeline.models.base_tool.base_tool import BaseTool
from mitopipeline.models.sample import Sample

class FastQCRunner(BaseTool):
    """
    This class provides the functionality to run FastQC and normalize I/O 
    within a single sample directory. 

    Attributes:
        working_dir (Path): The path to the directory of the input files.
        output_dir (Path): The directory to save the output files to.
        sample (Sample): The sample object that is being worked with.
        r1_output_stem (str): The stem of the R1 file - this will be appended to with file extensions for the outputs.
        r2_output_stem (str): The stem of the R2 file - this will be appended to with file extensions for the outputs.
        logger (Logger): The logger to use for this of the pipeline.
        tool_name (str): The name of the tool (by default will be fastqc).
        threads (int): The number of the threads for this process to use.
    """
    def __init__(self,
                 working_dir: Path,
                 output_dir: Path,
                 sample: Sample,
                 r1_output_stem: str,
                 r2_output_stem: str,
                 logger: Logger | None = None, 
                 tool_name: str = "fastqc",
                 threads: int = 4) -> None:
        """
        Instantiation method for the creation of the 
        FastQCRunner class. Calls the BaseTool instantiation methods and
        sets a few other object variables.

        Args:
            working_dir (Path): The path to the directory of the input files.
            output_dir (Path): The directory to save the output files to.
            sample (Sample): The sample object that is being worked with.
            r1_output_stem (str): The stem of the R1 file - this will be appended to with file extensions for the outputs.
            r2_output_stem (str): The stem of the R2 file - this will be appended to with file extensions for the outputs.
            logger (Logger): The logger to use for this of the pipeline.
            tool_name (str): The name of the tool (by default will be fastqc).
            threads (int): The number of the threads for this process to use.
            
        Returns: 
            None
        """
        # FastQCRunner inherits from BaseTool - here we invoke BaseTool's
        # instantiation function within FastQCRunner's instantiation function. 
        super().__init__(
            tool_name = tool_name,
            workind_dir = Path(working_dir),
            logger = logger
        )        
        
        # Setting the other class variables.
        self.output_dir = Path(output_dir)
        self.sample = sample
        self.r1_output_stem = r1_output_stem
        self.r2_output_stem = r2_output_stem
        self.threads = int(threads)
        
    def validate_inputs(self) -> None:
        """
        Validates the input R1 and R2 files for FastQC.
        This function appropriately logs the result,
        raises an error otherwise. 
        
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

        # Validating that the number of threads is valid.
        if self.threads <= 0:
            self.logger.error(f"{context} Please provide a valid number of threads (passed = {self.threads}).")
            sys.exit(1)
        if self.threads > os.cpu_count():
            self.logger.error(f"{context} Number of threads is greater than number of cores. Please reduce the number of threads to {os.cpu_count()}.")

    def build_command(self) -> list[str]:
        """
        Constructs the command to be executed.

        Args:
            None

        Returns:
            list[str]: The command to be executed.
        """
        # Creating the output directory if it does not exist.
        self.output_dir.mkdir(parents = True, exist_ok = True)

        # Constructing the command.
        return [
            "fastqc", 
            "--threads", 
            str(self.threads),
            str(self.sample.r1),
            str(self.sample.r2),
            "-o",
            str(self.output_dir)
        ]

    def postprocess_outputs(self) -> None:
        """
        Further processes the outputs of fastqc, namely:
        1. Checks if the outputs exist.
        2. Checks if the outputs are the same as the expected outputs.
        3. Moves the outputs to the expected location.

        Args: 
            None

        Returns:
            None
        """
        context = super()._log_context()

        # Obtaining the FastQC outputs generated by command execution. 
        sources = [
            *self._native_outputs(self.sample.r1),
            *self._native_outputs(self.sample.r2),
        ]

        # Obtaining the expected outputs.
        targets = self._expected_fastqc_outputs()

        # Iterating through the outputs and moving them to the expected location.
        for source, target in zip(sources, targets, strict=True):
            if not source.is_file():
                self.logger.error(f"{context} Native FastQC output missing: {source}")
            if source != target:
                shutil.move(str(source), str(target))

        # Logging the result.
        self.logger.info(f"{context} Successfully processed FastQC outputs.")

    def validate_outputs(self) -> None:
        """
        Validates the outputs of fastqc. 
        This function appropriately logs the result,
        raises an error otherwise. 

        Args:
            None

        Returns:
            None
        """
        context = super()._log_context()

        for path in self._expected_fastqc_outputs():
            if not path.is_file() or path.stat().st_size == 0:
                self.logger.error(f"{context} Missing or empty output: {path}")
                sys.exit(1)
    
    def _native_outputs(self, fastq: Path) -> tuple[Path, Path]:
        """
        This function generates the expected outputs for FastQC as a tuple.

        Args: 
            fastq (Path): The path to the FastQ file.
        
        Returns:
            tuple[Path, Path]: A tuple containing the paths to the FastQC outputs.
        """
        stem = self._strip_fastq_suffix(fastq)
        return (
            self.output_dir / f"{stem}_fastqc.html",
            self.output_dir / f"{stem}_fastqc.zip"
        )

    def _expected_fastqc_outputs(self) -> list[Path]:
        """
        Generates a list of expected outputs for FastQC.

        Args: 
            None

        Returns:
            list[Path]: A list of expected outputs.
        """
        return [
            self.output_dir / f"{self.r1_output_stem}_fastqc.html",
            self.output_dir / f"{self.r1_output_stem}_fastqc.zip",
            self.output_dir / f"{self.r2_output_stem}_fastqc.html",
            self.output_dir / f"{self.r2_output_stem}_fastqc.zip"
        ]

    @staticmethod
    def _strip_fastq_suffix(path: Path) -> str:
        """
        Removes the suffix from a FastQ file name.
        This function is static to allow for unit testing.

        Args:
            path (Path): The path to the FastQ file.

        Returns:
            str: The stem of the FastQ file.
        """
        # Getting the file name.
        name = path.name

        # Iterating through the suffixes and removing them.
        for suffix in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
            if name.endswith(suffix):
                return name[:-len(suffix)]

        # If no suffix is found, return the original name.
        return path.stem