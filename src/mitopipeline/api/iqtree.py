"""
iqtree.py

IQ-TREE wrapper for command building and execution -- infers a
maximum-likelihood phylogeny from an aligned multi-sequence FASTA
, for a batch-level tree across this run's accepted samples.
"""

# Imports
from __future__ import annotations
from logging import Logger
from pathlib import Path
import os
import sys
from mitopipeline.models.base_tool.base_tool import BaseTool

class IqtreeRunner(BaseTool):
    """
    Infers a maximum-likelihood phylogeny from an aligned FASTA with
    IQ-TREE.

    Attributes:
        working_dir (Path): The working directory for the command.
        alignment_fasta (Path): The aligned multi-sequence FASTA (MAFFT's own output).
        output_prefix (Path): The prefix IQ-TREE writes all of its own output files under (e.g. "{prefix}.treefile").
        sample_id (str): A label for this run, for logging context.
        sequence_type (str): IQ-TREE's -st value (e.g. "DNA").
        model (str): IQ-TREE's -m value (e.g. "MFP" for automatic model selection).
        ultrafast_bootstrap (int): The number of ultrafast bootstrap replicates (-B). IQ-TREE itself requires at least 1000 when this feature is used.
        sh_alrt (int): The number of SH-aLRT replicates (--alrt).
        random_seed (int): A fixed random seed (--seed), for reproducibility.
        threads (int): The number of threads for this process to use (-T).
        logger (Logger): The logger to use for this step of the pipeline.
        tool_name (str): The name of the tool (by default will be iqtree).
        treefile (Path): The real Newick tree IQ-TREE writes, derived from output_prefix.
    """
    def __init__(self,
                 working_dir: Path,
                 alignment_fasta: Path,
                 output_prefix: Path,
                 sample_id: str,
                 sequence_type: str = "DNA",
                 model: str = "MFP",
                 ultrafast_bootstrap: int = 1000,
                 sh_alrt: int = 1000,
                 random_seed: int = 12345,
                 threads: int = 4,
                 logger: Logger | None = None,
                 tool_name: str = "iqtree",
                 metrics_dir: Path | None = None) -> None:
        """
        Instantiation method for the creation of the IqtreeRunner
        class. Calls the BaseTool instantiation method and sets a few
        other object variables.

        Args:
            working_dir (Path): The working directory for the command.
            alignment_fasta (Path): The aligned multi-sequence FASTA.
            output_prefix (Path): The prefix IQ-TREE writes its own output under.
            sample_id (str): A label for this run, for logging context.
            sequence_type (str): IQ-TREE's -st value. Defaults to "DNA".
            model (str): IQ-TREE's -m value. Defaults to "MFP".
            ultrafast_bootstrap (int): Ultrafast bootstrap replicate count. Defaults to 1000.
            sh_alrt (int): SH-aLRT replicate count. Defaults to 1000.
            random_seed (int): A fixed random seed. Defaults to 12345.
            threads (int): The number of threads for this process to use.
            logger (Logger): The logger to use for this step of the pipeline.
            tool_name (str): The name of the tool (by default will be iqtree).
            metrics_dir (Path | None): Where to persist this run's compute metrics as JSON; None skips persistence. See BaseTool.

        Returns:
            None
        """
        super().__init__(tool_name = tool_name, working_dir = working_dir, logger = logger, metrics_dir = metrics_dir)
        self.alignment_fasta = Path(alignment_fasta)
        self.output_prefix = Path(output_prefix)
        self.sample_id = sample_id
        self.sequence_type = sequence_type
        self.model = model
        self.ultrafast_bootstrap = int(ultrafast_bootstrap)
        self.sh_alrt = int(sh_alrt)
        self.random_seed = int(random_seed)
        self.threads = int(threads)
        self.treefile = Path(f"{self.output_prefix}.treefile")

    def validate_inputs(self) -> None:
        """
        Validates the aligned input FASTA, bootstrap/SH-aLRT replicate
        counts, and thread count. This function appropriately logs
        the result, exits otherwise.

        Args:
            None

        Returns:
            None
        """
        context = self._log_context()

        if not self.alignment_fasta.is_file() or self.alignment_fasta.stat().st_size == 0:
            self.logger.error(f"{context} Missing or empty alignment FASTA: {self.alignment_fasta}.")
            sys.exit(1)

        # IQ-TREE itself requires at least 1000 replicates for ultrafast
        # bootstrap and SH-aLRT when either is actually used.
        if self.ultrafast_bootstrap < 1000:
            self.logger.error(f"{context} Ultrafast bootstrap replicates must be at least 1000 (passed = {self.ultrafast_bootstrap}).")
            sys.exit(1)
        if self.sh_alrt < 1000:
            self.logger.error(f"{context} SH-aLRT replicates must be at least 1000 (passed = {self.sh_alrt}).")
            sys.exit(1)

        if self.threads <= 0:
            self.logger.error(f"{context} Please provide a valid number of threads (passed = {self.threads}).")
            sys.exit(1)
        if os.cpu_count() and self.threads > os.cpu_count():
            self.logger.error(f"{context} Number of threads is greater than number of cores. Please reduce the number of threads to {os.cpu_count()}.")
            sys.exit(1)

    def version_command(self) -> list[str] | None:
        """Real, confirmed flag: `iqtree3 --version` prints e.g. "IQ-TREE version 3.1.2 for Linux x86 64-bit built ..." as its first line."""
        return ["iqtree3", "--version"]

    def build_command(self) -> list[str]:
        """
        Constructs the command to run IQ-TREE.

        Args:
            None

        Returns:
            list[str]: The command to run.
        """
        self.output_prefix.parent.mkdir(parents = True, exist_ok = True)

        return [
            "iqtree3",
            "-s", str(self.alignment_fasta),
            "-st", self.sequence_type,
            "-m", self.model,
            "-B", str(self.ultrafast_bootstrap),
            "--alrt", str(self.sh_alrt),
            "-T", str(self.threads),
            "--seed", str(self.random_seed),
            "--prefix", str(self.output_prefix),
            # Without this, IQ-TREE refuses to run at all when its own
            # prior checkpoint/output files already exist under the
            # same --prefix (exits in milliseconds with no useful
            # error beyond its own exit code) -- confirmed directly: a
            # second real run against the exact same alignment, same
            # prefix, failed this way even though the first run had
            # genuinely succeeded moments earlier. This pipeline always
            # wants a fresh, deterministic result from the given
            # alignment, never a resumed/cached one.
            "-redo",
        ]

    def validate_outputs(self) -> None:
        """
        Validates that IQ-TREE actually wrote its Newick tree file.
        This function appropriately logs the result, exits otherwise.

        Args:
            None

        Returns:
            None
        """
        context = self._log_context()
        if not self.treefile.is_file() or self.treefile.stat().st_size == 0:
            self.logger.error(f"{context} Missing or empty IQ-TREE output: {self.treefile}")
            sys.exit(1)
