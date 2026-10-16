"""
getorganelle.py

GetOrganelle wrapper for command building and execution.
"""

# Imports
from __future__ import annotations
from logging import Logger
from pathlib import Path
import os
import shutil
import sys
from mitopipeline.models.base_tool.base_tool import BaseTool
from mitopipeline.models.sample.sample import Sample

# The organelle types GetOrganelle itself supports via -F.
VALID_ORGANELLE_TYPES = {
    "animal_mt",
    "embplant_pt",
    "embplant_mt",
    "fungus_mt",
    "fungus_nr",
    "embplant_nr",
    "anonym",
}

class GetOrganelleRunner(BaseTool):
    """
    This class provides the functionality to run GetOrganelle and
    normalize I/O within a single sample directory.

    Attributes:
        working_dir (Path): The path to the directory of the input files.
        output_dir (Path): The directory to save the output files to.
        sample (Sample): The sample object that is being worked with.
        organelle_type (str): The GetOrganelle organelle type to assemble (e.g. "animal_mt").
        tool_options (dict): A dictionary of options to pass to the tool. A "seed_file" entry (-s) is what makes this a seeded run; its absence is a valid no-seed run, falling back to GetOrganelle's own default database for organelle_type.
        logger (Logger): The logger to use for this of the pipeline.
        tool_name (str): The name of the tool (by default will be getorganelle).
        threads (int): The number of the threads for this process to use.
        standardized_fasta (Path): The path to the standardized assembly FASTA output.
        standardized_gfa (Path): The path to the standardized assembly graph output.
    """
    def __init__(self,
                 working_dir: Path,
                 output_dir: Path,
                 sample: Sample,
                 organelle_type: str,
                 tool_options: dict | None = None,
                 threads: int = 4,
                 logger: Logger | None = None,
                 tool_name: str = "getorganelle") -> None:
        """
        Instantiation method for the creation of the
        GetOrganelleRunner class. Calls the BaseTool instantiation methods
        and sets a few other object variables.

        Args:
            working_dir (Path): The path to the directory of the input files.
            output_dir (Path): The directory to save the output files to.
            sample (Sample): The sample object that is being worked with.
            organelle_type (str): The GetOrganelle organelle type to assemble.
            tool_options (dict): A dictionary of options to pass to the tool.
            threads (int): The number of the threads for this process to use.
            logger (Logger): The logger to use for this of the pipeline.
            tool_name (str): The name of the tool (by default will be getorganelle).

        Returns:
            None
        """
        # Using the parent class constructor and setting the output file names.
        super().__init__(tool_name = tool_name, working_dir = working_dir, logger = logger)
        self.output_dir = Path(output_dir)
        self.sample = sample
        self.organelle_type = organelle_type
        self.tool_options = tool_options or {}
        self.threads = int(threads)
        self.standardized_fasta = self.output_dir / "data.fasta"
        self.standardized_gfa = self.output_dir / "graph.gfa"

    def validate_inputs(self) -> None:
        """
        Validates the organelle type, input files, and thread count for
        GetOrganelle. This function appropriately logs the result, exits
        otherwise.

        Args:
            None

        Returns:
            None
        """
        # Adding logger context.
        context = self._log_context()

        # Checking that the requested organelle type is supported.
        if self.organelle_type not in VALID_ORGANELLE_TYPES:
            self.logger.error(f"{context} Unsupported organelle type: {self.organelle_type}.")
            sys.exit(1)

        # Checking if input files exist.
        for path in (self.sample.r1, self.sample.r2):
            if not path.is_file():
                self.logger.error(f"{context} Input file does not exist: {path}.")
                sys.exit(1)

        # Checking that the requested thread count is valid for this machine.
        if self.threads <= 0:
            self.logger.error(f"{context} Please provide a valid number of threads (passed = {self.threads}).")
            sys.exit(1)
        if os.cpu_count() and self.threads > os.cpu_count():
            self.logger.error(f"{context} Number of threads is greater than number of cores. Please reduce the number of threads to {os.cpu_count()}.")
            sys.exit(1)

    def build_command(self) -> list[str]:
        """
        Constructs the command to run GetOrganelle, seeded or unseeded
        depending on tool_options.

        Args:
            None

        Returns:
            list[str]: The command to run GetOrganelle.
        """
        # Creating the output directory if it does not exist.
        self.output_dir.mkdir(parents = True, exist_ok = True)

        # Constructing the command.
        command = [
            "get_organelle_from_reads.py",
            "-1",
            str(self.sample.r1),
            "-2",
            str(self.sample.r2),
            "-o",
            str(self.output_dir),
            "-F",
            self.organelle_type,
            "-t",
            str(self.threads),
        ]

        # Returning the command with additional options (including -s, the
        # seed file, if a seed was selected for this sample).
        return self._add_options(command)

    def postprocess_outputs(self) -> None:
        """
        Standardizes GetOrganelle's own versioned output filenames to a
        fixed data.fasta/graph.gfa pair, then removes everything else
        GetOrganelle left behind in the output directory.

        Args:
            None

        Returns:
            None
        """
        # Obtaining the logging context.
        context = self._log_context()

        # Locating GetOrganelle's own native output files.
        fasta_candidates = sorted(self.output_dir.glob("*.path_sequence.fasta"))
        gfa_candidates = sorted(self.output_dir.glob("*.selected_graph.gfa"))

        # Handling an unexpected number of native assembly outputs.
        if len(fasta_candidates) != 1:
            self.logger.error(f"{context} Expected exactly one GetOrganelle path_sequence FASTA in {self.output_dir}; found {len(fasta_candidates)}.")
            sys.exit(1)
        if len(gfa_candidates) != 1:
            self.logger.error(f"{context} Expected exactly one GetOrganelle selected_graph GFA in {self.output_dir}; found {len(gfa_candidates)}.")
            sys.exit(1)

        # Standardizing the native outputs to fixed filenames.
        shutil.copy2(fasta_candidates[0], self.standardized_fasta)
        shutil.copy2(gfa_candidates[0], self.standardized_gfa)

        # Removing everything else GetOrganelle left behind.
        retained_names = {self.standardized_fasta.name, self.standardized_gfa.name, "assembly.log"}
        for child in self.output_dir.iterdir():
            if child.name in retained_names:
                continue
            if child.is_dir():
                shutil.rmtree(child)
            else:
                child.unlink()

    def validate_outputs(self) -> None:
        """
        Validates the standardized outputs of GetOrganelle. This function
        appropriately logs the result, exits otherwise.

        Args:
            None

        Returns:
            None
        """
        # Obtaining the logging context.
        context = self._log_context()

        # Iterating through each of the expected outputs - if they do not
        # exist or are empty, then log and exit.
        for path in (self.standardized_fasta, self.standardized_gfa):
            if not path.is_file() or path.stat().st_size == 0:
                self.logger.error(f"{context} Missing or empty assembly output: {path}")
                sys.exit(1)

    def _add_options(self, command: list[str]) -> list[str]:
        """
        Helper function which appends optional flags to an original command.

        Args:
            command (list[str]): The command to modify.

        Returns:
            list[str]: The modified command.
        """
        # A dictionary of options that require a value - the key is the
        # config name and the value is the true flag passed to GetOrganelle.
        # "seed_file" (-s) is the one that distinguishes a seeded run from
        # an unseeded one; its absence is a valid, intentional configuration.
        value_options = {
            "seed_file": "-s",
            "max_rounds": "-R",
            "kmer": "-k",
            "word_size": "-w",
            "pre_grouping": "-P",
            "max_reads": "--max-reads",
            "target_coverage": "--reduce-reads-for-coverage",
            "min_read_length": "--min-read-len",
            "max_read_length": "--max-read-len",
            "expected_max_size": "--expected-max-size",
            "expected_min_size": "--expected-min-size",
            "genes_file": "--genes",
            "exclude_fasta": "--exclude",
            "prefix": "--prefix",
            "round_output_prefix": "--out-per-round",
            "blast_path": "--which-blast",
            "bandage_path": "--which-bandage",
            "spades_path": "--which-spades",
            "disentangle_df": "--disentangle-df",
            "disentangle_time_limit": "--disentangle-time-limit",
        }

        # Iterating through each of the value options and appending values
        # that are not None to the command.
        for key, flag in value_options.items():
            value = self.tool_options.get(key)
            if value is not None:
                command.extend([flag, str(value)])

        # A dictionary of options that require a boolean value - the key is
        # the config name and the value is the true flag passed to GetOrganelle.
        boolean_options = {
            "overwrite": "--overwrite",
            "continue_run": "--continue",
            "fast_mode": "--fast",
            "reverse_lsc": "--reverse-lsc",
            "no_slim": "--no-slim",
            "keep_temp_files": "--keep-temp-files",
            "verbose": "--verbose",
        }

        # Iterating through each of the boolean options and appending flags
        # that are True to the command.
        for key, flag in boolean_options.items():
            if self.tool_options.get(key, False):
                command.append(flag)

        # Returning the modified command.
        return command
