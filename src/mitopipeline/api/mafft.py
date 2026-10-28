"""
mafft.py

MAFFT wrapper for command building and execution -- aligns a combined
multi-sequence FASTA (#344/#345's combined dataset) for downstream
phylogenetic inference (#346).
"""

# Imports
from __future__ import annotations
from logging import Logger
from pathlib import Path
import os
import shlex
import sys
from mitopipeline.models.base_tool.base_tool import BaseTool

class MafftRunner(BaseTool):
    """
    Aligns a multi-sequence FASTA with MAFFT.

    Attributes:
        working_dir (Path): The working directory for the command.
        input_fasta (Path): The unaligned, combined multi-sequence FASTA.
        output_fasta (Path): The path to write the aligned FASTA to.
        sample_id (str): A label for this alignment (e.g. the run/job id, since this is a batch-level dataset, not one sample's own output).
        strategy (str): MAFFT's alignment strategy. "auto" (MAFFT's own --auto) is the only one with a dedicated flag here; any other value is passed through as a literal flag (e.g. "--localpair").
        threads (int): The number of threads for this process to use.
        logger (Logger): The logger to use for this step of the pipeline.
        tool_name (str): The name of the tool (by default will be mafft).
    """
    def __init__(self,
                 working_dir: Path,
                 input_fasta: Path,
                 output_fasta: Path,
                 sample_id: str,
                 strategy: str = "auto",
                 threads: int = 4,
                 logger: Logger | None = None,
                 tool_name: str = "mafft",
                 metrics_dir: Path | None = None) -> None:
        """
        Instantiation method for the creation of the MafftRunner
        class. Calls the BaseTool instantiation method and sets a few
        other object variables.

        Args:
            working_dir (Path): The working directory for the command.
            input_fasta (Path): The unaligned, combined multi-sequence FASTA.
            output_fasta (Path): The path to write the aligned FASTA to.
            sample_id (str): A label for this alignment, for logging context.
            strategy (str): MAFFT's alignment strategy. Defaults to "auto".
            threads (int): The number of threads for this process to use.
            logger (Logger): The logger to use for this step of the pipeline.
            tool_name (str): The name of the tool (by default will be mafft).
            metrics_dir (Path | None): Where to persist this run's compute metrics as JSON; None skips persistence. See BaseTool.

        Returns:
            None
        """
        super().__init__(tool_name = tool_name, working_dir = working_dir, logger = logger, metrics_dir = metrics_dir)
        self.input_fasta = Path(input_fasta)
        self.output_fasta = Path(output_fasta)
        self.sample_id = sample_id
        self.strategy = strategy
        self.threads = int(threads)

    def validate_inputs(self) -> None:
        """
        Validates the input FASTA and thread count. This function
        appropriately logs the result, exits otherwise.

        Args:
            None

        Returns:
            None
        """
        context = self._log_context()

        if not self.input_fasta.is_file() or self.input_fasta.stat().st_size == 0:
            self.logger.error(f"{context} Missing or empty input FASTA: {self.input_fasta}.")
            sys.exit(1)

        if self.threads <= 0:
            self.logger.error(f"{context} Please provide a valid number of threads (passed = {self.threads}).")
            sys.exit(1)
        if os.cpu_count() and self.threads > os.cpu_count():
            self.logger.error(f"{context} Number of threads is greater than number of cores. Please reduce the number of threads to {os.cpu_count()}.")
            sys.exit(1)

    def version_command(self) -> list[str] | None:
        """Real, confirmed flag: `mafft --version` prints e.g. "v7.526 (2024/Apr/26)" (to stderr)."""
        return ["mafft", "--version"]

    def build_command(self) -> list[str]:
        """
        Constructs the command to run MAFFT. MAFFT itself writes its
        alignment to stdout, not a file -- wrapped in a literal shell
        invocation here (not shell=True on the whole pipeline, just
        this one real redirect) so the aligned output lands at
        self.output_fasta the same way every other tool wrapper
        writes a real output file, matching BaseTool's own run()
        contract.

        Args:
            None

        Returns:
            list[str]: The command to run.
        """
        self.output_fasta.parent.mkdir(parents = True, exist_ok = True)

        strategy_flag = "--auto" if self.strategy == "auto" else self.strategy
        mafft_command = [
            "mafft", strategy_flag, "--nuc", "--thread", str(self.threads), str(self.input_fasta),
        ]
        redirect = f"{shlex.join(mafft_command)} > {shlex.quote(str(self.output_fasta))}"
        return ["bash", "-c", redirect]

    def validate_outputs(self) -> None:
        """
        Validates that MAFFT actually wrote an aligned FASTA. This
        function appropriately logs the result, exits otherwise.

        Args:
            None

        Returns:
            None
        """
        context = self._log_context()
        if not self.output_fasta.is_file() or self.output_fasta.stat().st_size == 0:
            self.logger.error(f"{context} Missing or empty MAFFT alignment output: {self.output_fasta}")
            sys.exit(1)
