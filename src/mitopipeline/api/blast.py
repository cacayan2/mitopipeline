"""
blast.py

BLAST wrapper for command building and execution.
"""

# Imports
from __future__ import annotations
from logging import Logger
from pathlib import Path
import os
import subprocess
import sys
from mitopipeline.models.base_tool.base_tool import BaseTool

# The BLAST -task values this wrapper supports.
VALID_BLAST_TASKS = {
    "blastn",
    "megablast",
    "dc-megablast",
    "blastn-short",
}

# The modes this wrapper supports: a local on-disk database, or NCBI's
# remote service.
VALID_BLAST_MODES = {"local", "remote"}

# The full set of tabular columns this wrapper always requests, so a
# BLAST run's output matches exactly what blast_hits.parse_blast_tsv
# (#400) expects to parse -- the bare "-outfmt 6" default only yields
# the first 12 of these, silently dropping qcovs/staxids/sscinames/
# stitle and breaking the parser on every row (caught porting this
# from v0.1.0, which built this whole 16-column list but never passed
# it to BLAST itself).
BLAST_OUTPUT_COLUMNS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore",
    "qcovs", "staxids", "sscinames", "stitle",
]
DEFAULT_OUTPUT_FORMAT = "6 " + " ".join(BLAST_OUTPUT_COLUMNS)

class BlastRunner(BaseTool):
    """
    This class provides the functionality to run nucleotide BLAST
    (local or remote) against a query FASTA and write tabular output.

    Attributes:
        working_dir (Path): The working directory for the command.
        query_fasta (Path): The path to the query FASTA file.
        database (str): The BLAST database name (a local db name, or a remote NCBI db such as "nt").
        output_file (Path): The path to write tabular BLAST output to.
        sample (object): The sample object that is being worked with.
        mode (str): Either "local" or "remote".
        threads (int): The number of threads for this process to use (local mode only).
        task (str): The blastn -task variant to use (e.g. "megablast").
        output_format (str): The -outfmt string passed to BLAST. Defaults to format 6 with the full explicit column list blast_hits.py expects.
        evalue (float): The maximum e-value for a reported hit.
        max_target_seqs (int): The maximum number of aligned sequences to keep.
        max_hsps (int | None): The maximum number of HSPs per subject sequence to keep, or None to use BLAST's own default.
        perc_identity (float | None): The minimum percent identity for a reported hit, or None for no filter.
        query_coverage (float | None): The minimum percent query coverage per hit, or None for no filter.
        word_size (int | None): The BLAST word size, or None to use BLAST's own default.
        logger (Logger): The logger to use for this step of the pipeline.
        tool_name (str): The name of the tool (by default will be blast).
    """
    def __init__(self,
                 working_dir: Path,
                 query_fasta: Path,
                 database: str,
                 output_file: Path,
                 sample,
                 mode: str = "local",
                 threads: int = 4,
                 task: str = "blastn",
                 output_format: str = DEFAULT_OUTPUT_FORMAT,
                 evalue: float = 1e-5,
                 max_target_seqs: int = 50,
                 max_hsps: int | None = 1,
                 perc_identity: float | None = None,
                 query_coverage: float | None = None,
                 word_size: int | None = None,
                 logger: Logger | None = None,
                 tool_name: str = "blast") -> None:
        """
        Instantiation method for the creation of the BlastRunner class.
        Calls the BaseTool instantiation method and sets a few other
        object variables.

        Args:
            working_dir (Path): The working directory for the command.
            query_fasta (Path): The path to the query FASTA file.
            database (str): The BLAST database name.
            output_file (Path): The path to write tabular BLAST output to.
            sample (object): The sample object that is being worked with.
            mode (str): Either "local" or "remote".
            threads (int): The number of threads for this process to use.
            task (str): The blastn -task variant to use.
            output_format (str): The -outfmt string passed to BLAST.
            evalue (float): The maximum e-value for a reported hit.
            max_target_seqs (int): The maximum number of aligned sequences to keep.
            max_hsps (int | None): The maximum number of HSPs per subject sequence to keep.
            perc_identity (float | None): The minimum percent identity for a reported hit.
            query_coverage (float | None): The minimum percent query coverage per hit.
            word_size (int | None): The BLAST word size.
            logger (Logger): The logger to use for this step of the pipeline.
            tool_name (str): The name of the tool (by default will be blast).

        Returns:
            None
        """
        # Using the parent class constructor and setting the other object variables.
        super().__init__(tool_name = tool_name, working_dir = working_dir, logger = logger)
        self.query_fasta = Path(query_fasta)
        self.database = database
        self.output_file = Path(output_file)
        self.sample = sample
        self.mode = mode
        self.threads = int(threads)
        self.task = task
        self.output_format = output_format
        self.evalue = evalue
        self.max_target_seqs = max_target_seqs
        self.max_hsps = max_hsps
        self.perc_identity = perc_identity
        self.query_coverage = query_coverage
        self.word_size = word_size

    def validate_inputs(self) -> None:
        """
        Validates the query FASTA, BLAST mode, task, database, and
        thread count. This function appropriately logs the result,
        exits otherwise.

        Args:
            None

        Returns:
            None
        """
        # Adding logger context.
        context = self._log_context()

        # Checking that the query FASTA exists.
        if not self.query_fasta.is_file():
            self.logger.error(f"{context} Query FASTA does not exist: {self.query_fasta}.")
            sys.exit(1)

        # Checking that the requested mode and task are supported.
        if self.mode not in VALID_BLAST_MODES:
            self.logger.error(f"{context} Unsupported BLAST mode: {self.mode}.")
            sys.exit(1)
        if self.task not in VALID_BLAST_TASKS:
            self.logger.error(f"{context} Unsupported BLAST task: {self.task}.")
            sys.exit(1)

        # Checking that a local database actually exists before running.
        if self.mode == "local":
            result = subprocess.run(
                ["blastdbcmd", "-db", self.database, "-dbtype", "nucl", "-info"],
                capture_output = True, text = True,
            )
            if result.returncode != 0:
                self.logger.error(f"{context} Invalid BLAST database: {self.database}. {result.stderr.strip()}")
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
        Constructs the command to run BLAST, local or remote depending
        on mode.

        Args:
            None

        Returns:
            list[str]: The command to run BLAST.
        """
        # Creating the output file's parent directory if it does not exist.
        self.output_file.parent.mkdir(parents = True, exist_ok = True)

        # Constructing the base command.
        command = [
            "blastn",
            "-query", str(self.query_fasta),
            "-db", self.database,
            "-out", str(self.output_file),
            "-outfmt", self.output_format,
            "-task", self.task,
            "-evalue", str(self.evalue),
            "-max_target_seqs", str(self.max_target_seqs),
        ]

        # Remote runs pass -remote instead of a thread count (NCBI's
        # service controls its own concurrency).
        if self.mode == "remote":
            command.append("-remote")
        else:
            command.extend(["-num_threads", str(self.threads)])

        # Appending the optional filters, only when actually configured.
        if self.max_hsps is not None:
            command.extend(["-max_hsps", str(self.max_hsps)])
        if self.perc_identity is not None:
            command.extend(["-perc_identity", str(self.perc_identity)])
        if self.query_coverage is not None:
            command.extend(["-qcov_hsp_perc", str(self.query_coverage)])
        if self.word_size is not None:
            command.extend(["-word_size", str(self.word_size)])

        # Returning the constructed command.
        return command

    def validate_outputs(self) -> None:
        """
        Validates that BLAST actually wrote its tabular output file.
        This function appropriately logs the result, exits otherwise.

        Args:
            None

        Returns:
            None
        """
        # Obtaining the logging context.
        context = self._log_context()

        # A remote/local run that finds zero hits still writes an
        # empty output file -- that's a valid result, not a failure,
        # so only missing (not merely empty) is checked here.
        if not self.output_file.is_file():
            self.logger.error(f"{context} BLAST output not found: {self.output_file}.")
            sys.exit(1)
