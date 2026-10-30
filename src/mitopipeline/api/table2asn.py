"""
table2asn.py

table2asn wrapper for command building and execution -- the same
BaseTool contract every other real external tool in this project
uses (build_command()/version_command(), real compute-metrics
persistence via run(), never a standalone subprocess.run() call
outside that contract).
"""

# Imports
from __future__ import annotations
from logging import Logger
from pathlib import Path
import sys
from mitopipeline.models.base_tool.base_tool import BaseTool


class Table2AsnRunner(BaseTool):
    """
    Runs table2asn to validate and package one sample's real GenBank
    submission (FASTA + 5-column feature table + submission
    template) into a real .sqn/.val/.stats/.dr/.gbf set.

    Attributes:
        working_dir (Path): The working directory for the command -- table2asn resolves its own output basename relative to this, so it must be the sample's own output directory.
        fasta_path (Path): The real, NCBI-defline-tagged genome FASTA (.fsa).
        feature_table_path (Path): The real 5-column feature table (.tbl), basename-matched to fasta_path.
        submission_template_path (Path): The real submission template (.sbt).
        output_basename (Path): The real output basename -- table2asn derives .sqn/.val/.stats/.dr/.gbf from this plus their own extensions.
        logger (Logger): The logger to use for this step of the pipeline.
        tool_name (str): The name of the tool (by default "table2asn").
        binary_path (str): The real table2asn executable path (absolute, never resolved via $PATH -- Snakemake sanitizes PATH for conda-tagged rules, the same gotcha already documented for typst in this project).
    """
    def __init__(self,
                 working_dir: Path,
                 fasta_path: Path,
                 feature_table_path: Path,
                 submission_template_path: Path,
                 output_basename: Path,
                 logger: Logger | None = None,
                 tool_name: str = "table2asn",
                 binary_path: str = "table2asn",
                 metrics_dir: Path | None = None) -> None:
        """
        Instantiation method for the creation of the Table2AsnRunner class.

        Args:
            working_dir (Path): The working directory for the command.
            fasta_path (Path): The real, NCBI-defline-tagged genome FASTA.
            feature_table_path (Path): The real 5-column feature table, basename-matched to fasta_path.
            submission_template_path (Path): The real submission template.
            output_basename (Path): The real output basename for table2asn's own derived outputs.
            logger (Logger | None): The logger to use for this step of the pipeline.
            tool_name (str): The name of the tool. Defaults to "table2asn".
            binary_path (str): The real table2asn executable path. Defaults to "table2asn" (PATH lookup).
            metrics_dir (Path | None): Where to persist this run's compute metrics as JSON. None skips persistence. See BaseTool.

        Returns:
            None
        """
        super().__init__(tool_name = tool_name, working_dir = working_dir, logger = logger, metrics_dir = metrics_dir)
        self.fasta_path = Path(fasta_path)
        self.feature_table_path = Path(feature_table_path)
        self.submission_template_path = Path(submission_template_path)
        self.output_basename = Path(output_basename)
        self.binary_path = binary_path

    def validate_inputs(self) -> None:
        """
        Validates that the real FASTA, feature table, and submission
        template all exist, and that the FASTA/feature table share
        the real basename table2asn requires to pair them.

        Args:
            None

        Returns:
            None
        """
        context = self._log_context()
        if not self.fasta_path.is_file():
            self.logger.error(f"{context} FASTA input does not exist: {self.fasta_path}.")
            sys.exit(1)
        if not self.feature_table_path.is_file():
            self.logger.error(f"{context} Feature table does not exist: {self.feature_table_path}.")
            sys.exit(1)
        if not self.submission_template_path.is_file():
            self.logger.error(f"{context} Submission template does not exist: {self.submission_template_path}.")
            sys.exit(1)
        if self.fasta_path.stem != self.feature_table_path.stem:
            self.logger.error(f"{context} FASTA and feature table basenames must match for table2asn "
                              f"to pair them: {self.fasta_path.name} vs {self.feature_table_path.name}.")
            sys.exit(1)

    def version_command(self) -> list[str] | None:
        """Real, confirmed flag: `table2asn -version` prints e.g. "table2asn: 1.29.324"."""
        return [self.binary_path, "-version"]

    def build_command(self) -> list[str]:
        """
        Constructs the real table2asn command: -V vb requests a
        verbose validation report (the real source of .val/.stats),
        -Z requests the discrepancy report (.dr) -- both the real
        checks this project's own validation step reads.

        Args:
            None

        Returns:
            list[str]: The command to run table2asn.
        """
        self.output_basename.parent.mkdir(parents = True, exist_ok = True)
        # Every path must be absolute, not just relative to the real
        # process cwd: BaseTool.run() executes with cwd=self.working_dir
        # (the sample's own subdirectory, matching v0.1.0's real
        # convention), so a relative -o here resolves against THAT,
        # not the project root -- confirmed directly against a real
        # Snakemake run, where the rule's own relative JOB_DIR-based
        # path doubled up underneath working_dir and table2asn failed
        # with a bare "Error in command-line arguments."
        return [
            self.binary_path,
            "-t", str(self.submission_template_path.resolve()),
            "-i", str(self.fasta_path.resolve()),
            "-f", str(self.feature_table_path.resolve()),
            "-o", str(self.output_basename.with_suffix(".sqn").resolve()),
            "-V", "vb",
            "-Z",
        ]

    def validate_outputs(self) -> None:
        """
        Validates that table2asn actually produced its own real .val
        validation report -- the one output it always writes on a
        genuine run, whether or not that run found anything to flag.
        This is a bare "did the tool actually run" sanity check, not
        a PASS/FAIL judgment: a .val full of real reject/error
        messages is still a successful, informative run, read
        separately by table2asn_validation.py's own real severity
        parsing -- never conflated with a tool that silently produced
        nothing at all.

        Args:
            None

        Returns:
            None
        """
        context = self._log_context()
        validation_file = self.output_basename.with_suffix(".val")
        if not validation_file.is_file():
            self.logger.error(f"{context} table2asn did not produce its own .val report: {validation_file}")
            sys.exit(1)
