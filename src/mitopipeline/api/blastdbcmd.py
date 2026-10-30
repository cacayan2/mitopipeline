"""
blastdbcmd.py

blastdbcmd wrapper for command building and execution -- retrieves
one selected BLAST hit's reference FASTA from a local BLAST database
by accession, so later stages (contig validation, annotation
cross-checks) have the actual reference sequence on disk, not just
its tabular hit row.
"""

# Imports
from __future__ import annotations
from logging import Logger
from pathlib import Path
import sys
from mitopipeline.models.base_tool.base_tool import BaseTool

class BlastDbCmdRunner(BaseTool):
    """
    Retrieves one accession's FASTA sequence from a local nucleotide
    BLAST database with blastdbcmd.

    Attributes:
        working_dir (Path): The path to the directory of the input files.
        output_dir (Path): The directory to save the output files to.
        sample_id (str): The sample this retrieval is for (the hit this accession belongs to).
        accession (str): The selected BLAST hit's accession (sseqid) to retrieve.
        database (str): The local BLAST database name to retrieve from.
        logger (Logger): The logger to use for this step of the pipeline.
        tool_name (str): The name of the tool (by default will be blastdbcmd).
        output_fasta (Path): The path to the retrieved reference FASTA.
    """
    def __init__(self,
                 working_dir: Path,
                 output_dir: Path,
                 sample_id: str,
                 accession: str,
                 database: str,
                 logger: Logger | None = None,
                 tool_name: str = "blastdbcmd",
                 metrics_dir: Path | None = None) -> None:
        """
        Instantiation method for the creation of the BlastDbCmdRunner
        class. Calls the BaseTool instantiation method and sets a few
        other object variables.

        Args:
            working_dir (Path): The path to the directory of the input files.
            output_dir (Path): The directory to save the output files to.
            sample_id (str): The sample this retrieval is for.
            accession (str): The accession to retrieve.
            database (str): The local BLAST database name.
            logger (Logger): The logger to use for this step of the pipeline.
            tool_name (str): The name of the tool (by default will be blastdbcmd).
            metrics_dir (Path | None): Where to persist this run's compute metrics as JSON; None skips persistence. See BaseTool.

        Returns:
            None
        """
        super().__init__(tool_name = tool_name, working_dir = working_dir, logger = logger, metrics_dir = metrics_dir)
        self.output_dir = Path(output_dir)
        self.sample_id = sample_id
        self.accession = accession
        self.database = database
        # Sanitized the same way other exec stages sanitize filenames
        # derived from external, not-fully-trusted-format strings
        # (accessions can contain characters like '|' in older-style
        # NCBI identifiers) -- never embed the raw accession directly
        # in a path.
        safe_accession = "".join(c if c.isalnum() or c in ".-_" else "_" for c in accession)
        self.output_fasta = self.output_dir / f"{sample_id}.{safe_accession}.reference.fasta"

    def validate_inputs(self) -> None:
        """
        Validates the accession string and database name. This
        function appropriately logs the result, exits otherwise.

        Args:
            None

        Returns:
            None
        """
        context = self._log_context()

        if not self.accession or not self.accession.strip():
            self.logger.error(f"{context} No accession given to retrieve.")
            sys.exit(1)
        if not self.database or not self.database.strip():
            self.logger.error(f"{context} No BLAST database given to retrieve from.")
            sys.exit(1)

    def build_command(self) -> list[str]:
        """
        Constructs the command to retrieve one accession's FASTA from
        a local BLAST database.

        Args:
            None

        Returns:
            list[str]: The command to run.
        """
        self.output_dir.mkdir(parents = True, exist_ok = True)
        return [
            "blastdbcmd",
            "-db", self.database,
            "-entry", self.accession,
            "-out", str(self.output_fasta),
        ]

    def validate_outputs(self) -> None:
        """
        Validates that blastdbcmd actually retrieved a sequence. This
        function appropriately logs the result, exits otherwise -- an
        accession not present in the database is exactly this case
        (blastdbcmd itself may exit 0 but write nothing, or exit
        nonzero depending on version; checked here either way).

        Args:
            None

        Returns:
            None
        """
        context = self._log_context()
        if not self.output_fasta.is_file() or self.output_fasta.stat().st_size == 0:
            self.logger.error(f"{context} Accession '{self.accession}' was not retrieved from database "
                             f"'{self.database}' (missing or empty output): {self.output_fasta}")
            sys.exit(1)
