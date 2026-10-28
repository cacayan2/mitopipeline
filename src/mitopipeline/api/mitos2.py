"""
mitos2.py

MITOS2 (runmitos.py) wrapper for command building and execution --
annotates one assembled mitogenome FASTA with its gene features
(#409), so later stages (feature parsing, biological validation
gates) have a structured annotation to read.
"""

# Imports
from __future__ import annotations
from logging import Logger
from pathlib import Path
import sys
from mitopipeline.models.base_tool.base_tool import BaseTool

class Mitos2Runner(BaseTool):
    """
    Annotates one assembled mitogenome FASTA with MITOS2.

    Attributes:
        working_dir (Path): The path to the directory of the input files.
        output_dir (Path): The directory to save the output files to.
        sample_id (str): The sample this annotation is for.
        input_fasta (Path): The assembled mitogenome FASTA to annotate.
        genetic_code (int): The genetic code to use (MITOS2 -c).
        refdir (Path): The base directory containing the reference data (MITOS2 -R).
        refseqver (str): The reference data version, a subdirectory of refdir (MITOS2 -r).
        circular (bool): Whether to treat the sequence as circular. MITOS2's own default is True; False passes --linear.
        tool_options (dict): A dictionary of boolean options to pass to the tool (noplots, best, ncbicode, zip_output).
        logger (Logger): The logger to use for this step of the pipeline.
        tool_name (str): The name of the tool (by default will be mitos2).
        result_gff (Path): The path to MITOS2's own gene-feature GFF output.
        result_fas (Path): The path to MITOS2's own per-feature nucleotide FASTA output.
        result_faa (Path): The path to MITOS2's own per-feature protein FASTA output.
        result_bed (Path): The path to MITOS2's own BED output.
        result_table (Path): The path to MITOS2's own annotation table output (named "result.seq" by the tool itself).
        result_geneorder (Path): The path to MITOS2's own gene-order output.
    """
    def __init__(self,
                 working_dir: Path,
                 output_dir: Path,
                 sample_id: str,
                 input_fasta: Path,
                 genetic_code: int,
                 refdir: Path,
                 refseqver: str,
                 circular: bool = True,
                 tool_options: dict | None = None,
                 logger: Logger | None = None,
                 tool_name: str = "mitos2",
                 metrics_dir: Path | None = None) -> None:
        """
        Instantiation method for the creation of the Mitos2Runner
        class. Calls the BaseTool instantiation method and sets a few
        other object variables.

        Args:
            working_dir (Path): The path to the directory of the input files.
            output_dir (Path): The directory to save the output files to.
            sample_id (str): The sample this annotation is for.
            input_fasta (Path): The assembled mitogenome FASTA to annotate.
            genetic_code (int): The genetic code to use.
            refdir (Path): The base directory containing the reference data.
            refseqver (str): The reference data version, a subdirectory of refdir.
            circular (bool): Whether to treat the sequence as circular. Defaults to True (MITOS2's own default).
            tool_options (dict): A dictionary of boolean options to pass to the tool.
            logger (Logger): The logger to use for this step of the pipeline.
            tool_name (str): The name of the tool (by default will be mitos2).
            metrics_dir (Path | None): Where to persist this run's compute metrics as JSON; None skips persistence. See BaseTool.

        Returns:
            None
        """
        super().__init__(tool_name = tool_name, working_dir = working_dir, logger = logger, metrics_dir = metrics_dir)
        self.output_dir = Path(output_dir)
        self.sample_id = sample_id
        self.input_fasta = Path(input_fasta)
        self.genetic_code = genetic_code
        self.refdir = Path(refdir)
        self.refseqver = refseqver
        self.circular = circular
        self.tool_options = dict(tool_options or {})

        self.result_gff = self.output_dir / f"{sample_id}.annotation.gff"
        self.result_fas = self.output_dir / f"{sample_id}.annotation.fas"
        self.result_faa = self.output_dir / f"{sample_id}.annotation.faa"
        self.result_bed = self.output_dir / f"{sample_id}.annotation.bed"
        self.result_table = self.output_dir / f"{sample_id}.annotation.seq"
        self.result_geneorder = self.output_dir / f"{sample_id}.annotation.geneorder"

        # runmitos.py's own real, fixed output filenames (the tool
        # offers no output-prefix option) -- renamed to the sample-
        # prefixed paths above in postprocess_outputs(), matching the
        # same convention already applied to every other stage's own
        # real output (fastp, FastQC, GetOrganelle).
        self._native_names = {
            "result.gff": self.result_gff,
            "result.fas": self.result_fas,
            "result.faa": self.result_faa,
            "result.bed": self.result_bed,
            "result.seq": self.result_table,
            "result.geneorder": self.result_geneorder,
            "result.mitos": self.output_dir / f"{sample_id}.annotation.mitos",
            "result.png": self.output_dir / f"{sample_id}.annotation.png",
            "ignored.mitos": self.output_dir / f"{sample_id}.annotation.ignored.mitos",
            "stst.dat": self.output_dir / f"{sample_id}.annotation.stst.dat",
        }

    def validate_inputs(self) -> None:
        """
        Validates the input FASTA, genetic code, and taxonomy
        reference directory for MITOS2. This function appropriately
        logs the result, exits otherwise.

        Args:
            None

        Returns:
            None
        """
        context = self._log_context()

        if not self.input_fasta.is_file() or self.input_fasta.stat().st_size == 0:
            self.logger.error(f"{context} Missing or empty input FASTA: {self.input_fasta}.")
            sys.exit(1)

        if not isinstance(self.genetic_code, int) or self.genetic_code <= 0:
            self.logger.error(f"{context} Invalid genetic code: {self.genetic_code!r}.")
            sys.exit(1)

        if not self.refseqver or not self.refseqver.strip():
            self.logger.error(f"{context} No reference data version (refseqver) given.")
            sys.exit(1)

        # The real runmitos.py checks exactly this combined path is a
        # directory (confirmed against the installed source: refdir =
        # os.path.join(args.refdir, args.refseqver)) -- checked here so
        # a bad taxonomy reference fails fast, not mid-run.
        reference_path = self.refdir / self.refseqver
        if not reference_path.is_dir():
            self.logger.error(f"{context} No such reference data directory: {reference_path}.")
            sys.exit(1)

    def build_command(self) -> list[str]:
        """
        Constructs the command to run MITOS2 (runmitos.py).

        Args:
            None

        Returns:
            list[str]: The command to run.
        """
        # Unlike GetOrganelle, runmitos.py requires -o to already
        # exist (confirmed against the installed source: it calls
        # parser.error("no such directory %s") when it doesn't) --
        # so the output directory is created here, not left for the
        # tool itself.
        self.output_dir.mkdir(parents = True, exist_ok = True)

        command = [
            "runmitos.py",
            "-i", str(self.input_fasta),
            "-c", str(self.genetic_code),
            "-o", str(self.output_dir),
            "-R", str(self.refdir),
            "-r", self.refseqver,
        ]

        # runmitos.py's own default is circular=True; --linear (dest
        # "circular", store_false) is the only way to flip it to
        # False (confirmed against the installed source) -- so it's
        # only ever passed for a non-circular sequence.
        if not self.circular:
            command.append("--linear")

        return self._add_options(command)

    def postprocess_outputs(self) -> None:
        """
        Renames MITOS2's own real, fixed-name output files (e.g.
        "result.gff") into sample-prefixed equivalents (e.g.
        "{sample_id}.annotation.gff") -- runmitos.py offers no
        output-prefix option of its own, so this is done directly
        against whatever it actually wrote. A native file that was
        never produced (e.g. no --zip, no --best) is skipped, not an
        error; already-renamed files (a resumed/idempotent re-run) are
        left as-is.

        Args:
            None

        Returns:
            None
        """
        for native_name, destination in self._native_names.items():
            source = self.output_dir / native_name
            if source.is_file():
                source.rename(destination)

    def validate_outputs(self) -> None:
        """
        Validates that MITOS2 actually produced its gene-feature
        outputs. This function appropriately logs the result, exits
        otherwise.

        Args:
            None

        Returns:
            None
        """
        context = self._log_context()
        for path in (self.result_gff, self.result_fas):
            if not path.is_file() or path.stat().st_size == 0:
                self.logger.error(f"{context} Missing or empty MITOS2 annotation output: {path}")
                sys.exit(1)

    def _add_options(self, command: list[str]) -> list[str]:
        """
        Helper function which appends optional boolean flags to an
        original command.

        Args:
            command (list[str]): The command to modify.

        Returns:
            list[str]: The modified command.
        """
        # A dictionary of boolean options -- the key is the config
        # name and the value is the true flag passed to MITOS2.
        boolean_options = {
            "noplots": "--noplots",
            "best": "--best",
            "ncbicode": "--ncbicode",
            "zip_output": "--zip",
        }

        for key, flag in boolean_options.items():
            if self.tool_options.get(key, False):
                command.append(flag)

        return command
