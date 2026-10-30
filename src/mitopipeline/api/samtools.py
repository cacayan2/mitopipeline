"""
samtools.py

samtools wrapper for command building and execution -- turns
minimap2's rescue-mapping SAM into a sorted, indexed BAM, then extracts
the reads that actually mapped as rescued paired FASTQ.

One wrapper, three operations (sort/index/fastq), not three separate
classes: they share the same sample-scoped output paths and chain in a
fixed order (sort -> index -> fastq), so one class with an `operation`
switch keeps that relationship explicit instead of three classes that
would otherwise need to agree on the same paths independently.
"""

# Imports
from __future__ import annotations
from logging import Logger
from pathlib import Path
import os
import sys
from mitopipeline.models.base_tool.base_tool import BaseTool
from mitopipeline.models.sample.sample import Sample

VALID_OPERATIONS = {"sort", "index", "fastq"}

class SamtoolsRunner(BaseTool):
    """
    Sorts, indexes, or extracts rescued reads from a rescue-mapping BAM
    with samtools, depending on `operation`.

    Attributes:
        working_dir (Path): The path to the directory of the input files.
        output_dir (Path): The directory to save the output files to.
        sample (Sample): The sample object that is being worked with.
        operation (str): Which samtools step to run -- "sort", "index", or "fastq".
        input_path (Path): The file this operation reads -- the mapping SAM for "sort", the sorted BAM for "index"/"fastq".
        logger (Logger): The logger to use for this step of the pipeline.
        tool_name (str): The name of the tool (by default will be samtools).
        threads (int): The number of threads for this process to use.
        sorted_bam (Path): The path to the coordinate-sorted BAM (sort's output; index's and fastq's default input).
        bam_index (Path): The path to the BAM index (index's output).
        rescued_r1 (Path): The path to the rescued, mapped R1 FASTQ (fastq's output).
        rescued_r2 (Path): The path to the rescued, mapped R2 FASTQ (fastq's output).
    """
    def __init__(self,
                 working_dir: Path,
                 output_dir: Path,
                 sample: Sample,
                 operation: str,
                 input_path: Path,
                 logger: Logger | None = None,
                 tool_name: str = "samtools",
                 threads: int = 4,
                 metrics_dir: Path | None = None) -> None:
        """
        Instantiation method for the creation of the SamtoolsRunner
        class. Calls the BaseTool instantiation methods and sets a few
        other object variables.

        Args:
            working_dir (Path): The path to the directory of the input files.
            output_dir (Path): The directory to save the output files to.
            sample (Sample): The sample object that is being worked with.
            operation (str): Which samtools step to run -- "sort", "index", or "fastq".
            input_path (Path): The file this operation reads.
            logger (Logger): The logger to use for this step of the pipeline.
            tool_name (str): The name of the tool (by default will be samtools).
            threads (int): The number of threads for this process to use.
            metrics_dir (Path | None): Where to persist this run's compute metrics as JSON; None skips persistence. See BaseTool.

        Returns:
            None
        """
        super().__init__(tool_name = tool_name, working_dir = working_dir, logger = logger, metrics_dir = metrics_dir)
        self.output_dir = Path(output_dir)
        self.sample = sample
        self.operation = operation
        self.input_path = Path(input_path)
        self.threads = int(threads)
        self.sorted_bam = self.output_dir / f"{self.sample.sample_id}.rescue.sorted.bam"
        self.bam_index = Path(str(self.sorted_bam) + ".bai")
        # Plain .fastq, NOT .fastq.gz: confirmed against a real run that
        # this samtools version (1.3.1) writes uncompressed text to -1/-2
        # regardless of the filename given -- a .gz-named file holding
        # plain text made GetOrganelle try to gunzip it and fail on an
        # "empty file" (it extracted nothing from genuinely-uncompressed
        # data). The filename now matches what's actually on disk.
        self.rescued_r1 = self.output_dir / f"{self.sample.sample_id}.rescue.R1.fastq"
        self.rescued_r2 = self.output_dir / f"{self.sample.sample_id}.rescue.R2.fastq"

    def validate_inputs(self) -> None:
        """
        Validates the operation name, input file, and thread count for
        samtools. This function appropriately logs the result, exits
        otherwise.

        Args:
            None

        Returns:
            None
        """
        context = self._log_context()

        if self.operation not in VALID_OPERATIONS:
            self.logger.error(f"{context} Unsupported samtools operation: {self.operation}. Must be one of {sorted(VALID_OPERATIONS)}.")
            sys.exit(1)

        if not self.input_path.is_file():
            self.logger.error(f"{context} Input file does not exist: {self.input_path}.")
            sys.exit(1)

        if self.threads <= 0:
            self.logger.error(f"{context} Please provide a valid number of threads (passed = {self.threads}).")
            sys.exit(1)
        if os.cpu_count() and self.threads > os.cpu_count():
            self.logger.error(f"{context} Number of threads is greater than number of cores. Please reduce the number of threads to {os.cpu_count()}.")
            sys.exit(1)

    def version_command(self) -> list[str] | None:
        """Real, confirmed flag: `samtools --version` prints e.g. "samtools 1.24" as its first line."""
        return ["samtools", "--version"]

    def build_command(self) -> list[str]:
        """
        Constructs the samtools command for this runner's operation.

        Args:
            None

        Returns:
            list[str]: The command to run.
        """
        self.output_dir.mkdir(parents = True, exist_ok = True)

        if self.operation == "sort":
            return ["samtools", "sort", "-@", str(self.threads), "-o", str(self.sorted_bam), str(self.input_path)]
        if self.operation == "index":
            # NOT passing -@ here: confirmed against a real failure
            # ("index: invalid option -- '@'", samtools 1.3.1) that
            # unlike sort, index doesn't accept a thread count on every
            # samtools version this wrapper might run against. Index
            # is comparatively cheap regardless.
            return ["samtools", "index", str(self.input_path)]
        # "fastq": -n keeps read names intact (paired-end output needs
        # matching names); -0/-s discard unpaired/singleton reads --
        # rescue only ever wants properly paired mapped reads. Also NOT
        # passing -@, for the same real-failure reason as index above
        # ("fastq: invalid option -- '@'" on the same samtools version).
        return [
            "samtools", "fastq", "-n",
            "-1", str(self.rescued_r1), "-2", str(self.rescued_r2),
            "-0", os.devnull, "-s", os.devnull,
            str(self.input_path),
        ]

    def validate_outputs(self) -> None:
        """
        Validates that this runner's operation actually produced its
        expected output. This function appropriately logs the result,
        exits otherwise.

        Args:
            None

        Returns:
            None
        """
        context = self._log_context()
        if self.operation == "sort":
            expected = [self.sorted_bam]
        elif self.operation == "index":
            expected = [self.bam_index]
        else:
            expected = [self.rescued_r1, self.rescued_r2]

        for path in expected:
            if not path.is_file() or path.stat().st_size == 0:
                self.logger.error(f"{context} Missing or empty output for operation '{self.operation}': {path}")
                sys.exit(1)
