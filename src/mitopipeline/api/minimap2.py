"""
minimap2.py

minimap2 wrapper for command building and execution -- maps a sample's
reads against a reference FASTA during assembly rescue, so samtools can
sort/index/extract the reads that actually mapped (#423/#424/#425).
"""

# Imports
from __future__ import annotations
from logging import Logger
from pathlib import Path
import os
import sys
from mitopipeline.models.base_tool.base_tool import BaseTool
from mitopipeline.models.sample.sample import Sample

class Minimap2Runner(BaseTool):
    """
    Maps a sample's paired reads against a reference FASTA with
    minimap2, writing a SAM file for samtools to process next.

    Attributes:
        working_dir (Path): The path to the directory of the input files.
        output_dir (Path): The directory to save the output files to.
        sample (Sample): The sample object that is being worked with.
        reference_fasta (Path): The reference FASTA to map reads against (the rescue target -- typically the sample's own fragmented/missing assembly attempt or its seed reference).
        preset (str): minimap2's -x preset (default "sr", short-read paired-end).
        logger (Logger): The logger to use for this step of the pipeline.
        tool_name (str): The name of the tool (by default will be minimap2).
        threads (int): The number of threads for this process to use.
        output_sam (Path): The path to the mapping output SAM file.
    """
    def __init__(self,
                 working_dir: Path,
                 output_dir: Path,
                 sample: Sample,
                 reference_fasta: Path,
                 preset: str = "sr",
                 logger: Logger | None = None,
                 tool_name: str = "minimap2",
                 threads: int = 4,
                 metrics_dir: Path | None = None) -> None:
        """
        Instantiation method for the creation of the Minimap2Runner
        class. Calls the BaseTool instantiation methods and sets a few
        other object variables.

        Args:
            working_dir (Path): The path to the directory of the input files.
            output_dir (Path): The directory to save the output files to.
            sample (Sample): The sample object that is being worked with.
            reference_fasta (Path): The reference FASTA to map reads against.
            preset (str): minimap2's -x preset.
            logger (Logger): The logger to use for this step of the pipeline.
            tool_name (str): The name of the tool (by default will be minimap2).
            threads (int): The number of threads for this process to use.
            metrics_dir (Path | None): Where to persist this run's compute metrics as JSON; None skips persistence. See BaseTool.

        Returns:
            None
        """
        super().__init__(tool_name = tool_name, working_dir = working_dir, logger = logger, metrics_dir = metrics_dir)
        self.output_dir = Path(output_dir)
        self.sample = sample
        self.reference_fasta = Path(reference_fasta)
        self.preset = preset
        self.threads = int(threads)
        self.output_sam = self.output_dir / f"{self.sample.sample_id}.rescue_mapping.sam"

    def validate_inputs(self) -> None:
        """
        Validates the reference FASTA, read files, and thread count for
        minimap2. This function appropriately logs the result, exits
        otherwise.

        Args:
            None

        Returns:
            None
        """
        context = self._log_context()

        if not self.reference_fasta.is_file():
            self.logger.error(f"{context} Reference FASTA does not exist: {self.reference_fasta}.")
            sys.exit(1)

        for path in (self.sample.r1, self.sample.r2):
            if not path.is_file():
                self.logger.error(f"{context} Input file does not exist: {path}.")
                sys.exit(1)

        if self.threads <= 0:
            self.logger.error(f"{context} Please provide a valid number of threads (passed = {self.threads}).")
            sys.exit(1)
        if os.cpu_count() and self.threads > os.cpu_count():
            self.logger.error(f"{context} Number of threads is greater than number of cores. Please reduce the number of threads to {os.cpu_count()}.")
            sys.exit(1)

    def build_command(self) -> list[str]:
        """
        Constructs the command to run minimap2 against this sample's
        paired reads.

        Args:
            None

        Returns:
            list[str]: The command to run minimap2.
        """
        self.output_dir.mkdir(parents = True, exist_ok = True)
        return [
            "minimap2",
            "-ax", self.preset,
            "-t", str(self.threads),
            "-o", str(self.output_sam),
            str(self.reference_fasta),
            str(self.sample.r1),
            str(self.sample.r2),
        ]

    def validate_outputs(self) -> None:
        """
        Validates that minimap2 actually produced a mapping output.
        This function appropriately logs the result, exits otherwise.

        Args:
            None

        Returns:
            None
        """
        context = self._log_context()
        if not self.output_sam.is_file() or self.output_sam.stat().st_size == 0:
            self.logger.error(f"{context} Missing or empty mapping output: {self.output_sam}")
            sys.exit(1)
