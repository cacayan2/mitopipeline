"""
getorganelle.py

GetOrganelle wrapper for command building and execution.
"""

# Imports
from __future__ import annotations
from logging import Logger
from pathlib import Path
import json
import os
import shutil
import sys
from dataclasses import asdict, dataclass
from mitopipeline.models.base_tool.base_tool import BaseTool
from mitopipeline.models.sample.sample import Sample
from mitopipeline.utils.seed.models import SeedSelectionResult
from mitopipeline.api.fastqc import FastQCRunner

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
                 tool_name: str = "getorganelle",
                 metrics_dir: Path | None = None) -> None:
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
            metrics_dir (Path | None): Where to persist this run's compute metrics (CPU, memory, runtime) as JSON; None skips persistence. See BaseTool.

        Returns:
            None
        """
        # Using the parent class constructor and setting the output file names.
        super().__init__(tool_name = tool_name, working_dir = working_dir, logger = logger, metrics_dir = metrics_dir)
        self.output_dir = Path(output_dir)
        self.sample = sample
        self.organelle_type = organelle_type

        # Copying the caller's tool_options (never mutate what was passed
        # in, matching build_runner_from_seed_result's own stated care)
        # and forcing keep_temp_files on: without GetOrganelle's own
        # --keep-temp, the extended_*_paired.fq reads it recruits (the
        # ones actually handed to SPAdes for assembly) are deleted by
        # GetOrganelle itself before this wrapper ever sees them, so the
        # recruited-reads QC step below would have nothing to QC.
        self.tool_options = dict(tool_options or {})
        self.tool_options.setdefault("keep_temp_files", True)

        self.threads = int(threads)
        self.standardized_fasta = self.output_dir / "data.fasta"
        self.standardized_gfa = self.output_dir / "graph.gfa"
        self.recruited_reads_fastqc_dir = self.output_dir / "recruited_reads_fastqc"
        self.recruited_reads_qc_record_path = self.output_dir / f"{sample.sample_id}.recruited_reads_qc.json"

        # The classification discover_assembly_artifacts() produced,
        # captured in postprocess_outputs() before it cleans up
        # GetOrganelle's native output filenames. A caller that needs
        # the circular/fragmented/missing status after run() returns
        # must read this rather than re-deriving it, since by then the
        # native *.path_sequence.fasta/*.selected_graph.gfa files
        # discover_assembly_artifacts() looks for no longer exist.
        self.assembly_status: AssemblyArtifacts | None = None

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
        # NOT pre-creating self.output_dir here: get_organelle_from_reads.py
        # refuses to run if its -o directory already exists ("Please use a
        # new output directory, or use '--continue'/'--overwrite'"), and it
        # creates that exact leaf directory itself on a successful run.
        # Pre-creating it here made every real invocation of this runner
        # fail outright.
        #
        # Its own directory creation is a bare os.mkdir(), NOT recursive --
        # confirmed directly against a real failure (FileNotFoundError on
        # a two-level-deep output_dir, e.g. a rescue attempt's
        # "assembly/attempt0" where "assembly/" didn't exist yet; the
        # single-level paths every prior real run happened to use never
        # exposed this). So the immediate PARENT must already exist, even
        # though the leaf itself must not -- hence parent.mkdir(), never
        # output_dir.mkdir().
        self.output_dir.parent.mkdir(parents = True, exist_ok = True)

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

    def _recruited_reads_paths(self) -> tuple[Path, Path]:
        """
        The paired reads GetOrganelle recruited as mitochondrial-like
        during its own extending process -- the exact reads handed to
        SPAdes for assembly, before assembly. Confirmed against the
        installed get_organelle_from_reads.py: written as
        "{prefix}extended_1_paired.fq" / "{prefix}extended_2_paired.fq"
        in the output directory, and never removed by GetOrganelle
        itself (unlike the pre-split "{prefix}extended_1/2.fq", which
        it deletes once keep_temp_files is not set).

        Args:
            None

        Returns:
            tuple[Path, Path]: The expected (r1, r2) recruited-reads paths, whether or not they actually exist yet.
        """
        prefix = str(self.tool_options.get("prefix") or "")
        return (
            self.output_dir / f"{prefix}extended_1_paired.fq",
            self.output_dir / f"{prefix}extended_2_paired.fq",
        )

    def _run_recruited_reads_qc(self) -> None:
        """
        Runs FastQC on the reads GetOrganelle recruited as
        mitochondrial-like before assembly -- a real, early look at
        whether what's about to be assembled actually looks like
        clean, usable sequence, not just trusting the final
        assembly's own quality.

        Always writes a RecruitedReadsQCRecord (#413) -- including
        the absent/empty/malformed cases -- so a caller never has to
        re-derive what happened from log lines alone; this is a
        supplementary QC artifact, not a correctness gate on the
        assembly itself, so none of these cases fail the whole run.
        Compression is never checked here: GetOrganelle's own
        extended_*_paired.fq files are always plain text (confirmed
        against the installed get_organelle_from_reads.py), never
        gzipped, so there is nothing to validate on that axis.

        Args:
            None

        Returns:
            None
        """
        context = self._log_context()
        r1, r2 = self._recruited_reads_paths()

        if not (r1.is_file() and r2.is_file()):
            self.logger.warning(f"{context} Recruited-reads QC skipped: expected {r1.name}/{r2.name} not found in {self.output_dir} (keep_temp_files may not have taken effect).")
            self._write_recruited_reads_qc_record(RecruitedReadsQCRecord(
                sample_id = self.sample.sample_id, source_output_dir = str(self.output_dir),
                r1_path = None, r2_path = None, r1_read_count = None, r2_read_count = None,
                paired = None, qc_ran = False, fastqc_dir = None,
                reason = "recruited-reads files not found",
            ))
            return

        if r1.stat().st_size == 0 or r2.stat().st_size == 0:
            self.logger.warning(f"{context} Recruited-reads QC skipped: {r1.name}/{r2.name} present but empty.")
            self._write_recruited_reads_qc_record(RecruitedReadsQCRecord(
                sample_id = self.sample.sample_id, source_output_dir = str(self.output_dir),
                r1_path = str(r1), r2_path = str(r2), r1_read_count = None, r2_read_count = None,
                paired = None, qc_ran = False, fastqc_dir = None,
                reason = "recruited-reads file(s) empty",
            ))
            return

        r1_count = _count_fastq_records(r1)
        r2_count = _count_fastq_records(r2)
        paired = None if (r1_count is None or r2_count is None) else (r1_count == r2_count)

        if r1_count is None or r2_count is None:
            self.logger.warning(f"{context} Recruited-reads {r1.name}/{r2.name} line count is not a multiple of 4 (malformed FASTQ); running QC anyway, read counts unavailable.")
        elif not paired:
            self.logger.warning(f"{context} Recruited-reads pairing mismatch: {r1.name} has {r1_count} reads, {r2.name} has {r2_count}.")

        recruited_sample = Sample(sample_id = self.sample.sample_id, r1 = r1, r2 = r2)
        qc_runner = FastQCRunner(
            working_dir = self.output_dir,
            output_dir = self.recruited_reads_fastqc_dir,
            sample = recruited_sample,
            r1_output_stem = "recruited_R1",
            r2_output_stem = "recruited_R2",
            logger = self.logger,
            threads = self.threads,
        )
        qc_ran = False
        reason = "qc completed"
        try:
            result = qc_runner.run()
            if result.success:
                qc_ran = True
            else:
                reason = f"fastqc exited {result.return_code}"
                self.logger.warning(f"{context} FastQC on recruited reads exited with return code {result.return_code}; continuing, this is supplementary QC, not an assembly correctness check.")
        except SystemExit:
            reason = "fastqc failed validation"
            self.logger.warning(f"{context} FastQC on recruited reads failed validation; continuing, this is supplementary QC, not an assembly correctness check.")

        self._write_recruited_reads_qc_record(RecruitedReadsQCRecord(
            sample_id = self.sample.sample_id, source_output_dir = str(self.output_dir),
            r1_path = str(r1), r2_path = str(r2), r1_read_count = r1_count, r2_read_count = r2_count,
            paired = paired, qc_ran = qc_ran,
            fastqc_dir = str(self.recruited_reads_fastqc_dir) if qc_ran else None,
            reason = reason,
        ))

    def _write_recruited_reads_qc_record(self, record: RecruitedReadsQCRecord) -> None:
        """
        Persists a RecruitedReadsQCRecord as JSON in this attempt's
        own output_dir -- already attempt-specific (each rescue
        attempt gets its own output_dir), so this never collides or
        gets overwritten across rescue retries.

        Args:
            record (RecruitedReadsQCRecord): The record to persist.

        Returns:
            None
        """
        with self.recruited_reads_qc_record_path.open("w", encoding = "utf-8") as handle:
            json.dump(asdict(record), handle, indent = 4)

    def postprocess_outputs(self) -> None:
        """
        Standardizes GetOrganelle's own versioned output filenames to a
        fixed data.fasta/graph.gfa pair, runs FastQC on the reads it
        recruited before assembly, then removes everything else
        GetOrganelle left behind in the output directory.

        Args:
            None

        Returns:
            None
        """
        # Obtaining the logging context.
        context = self._log_context()

        # Classifying the result BEFORE anything below renames/deletes
        # GetOrganelle's native output -- discover_assembly_artifacts()
        # globs for the same native filenames, so it must run while
        # they still exist. Callers read this back via self.assembly_status.
        self.assembly_status = discover_assembly_artifacts(self.output_dir, self.logger)

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

        # QC'ing the recruited reads BEFORE the cleanup below deletes
        # them -- they're GetOrganelle's own intermediate files, not
        # part of the standardized output, so they don't survive past
        # this point otherwise.
        self._run_recruited_reads_qc()

        # Removing everything else GetOrganelle left behind.
        retained_names = {self.standardized_fasta.name, self.standardized_gfa.name, "assembly.log",
                           self.recruited_reads_fastqc_dir.name, self.recruited_reads_qc_record_path.name}
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
            # The real GetOrganelle flag is --keep-temp (confirmed against
            # the installed get_organelle_from_reads.py's own --help and
            # argparse dest="keep_temp_files"); --keep-temp-files isn't a
            # recognized flag at all and GetOrganelle would reject it.
            "keep_temp_files": "--keep-temp",
            "verbose": "--verbose",
        }

        # Iterating through each of the boolean options and appending flags
        # that are True to the command.
        for key, flag in boolean_options.items():
            if self.tool_options.get(key, False):
                command.append(flag)

        # Returning the modified command.
        return command

# The status an assembly-artifact discovery can report.
AssemblyDiscoveryStatus = str   # one of "circular", "fragmented", "missing"

@dataclass(frozen = True, slots = True)
class AssemblyArtifacts:
    """
    The outcome of discovering GetOrganelle's native assembly output
    in one sample's output directory.

    Attributes:
        status (str): One of "circular" (a complete, circularized assembly was found), "fragmented" (at least one candidate was found, but none is marked complete/circular), or "missing" (no candidate at all).
        fasta_path (Path | None): The path to the chosen native FASTA candidate, or None when status is "missing".
        gfa_path (Path | None): The path to the chosen native GFA candidate, or None when status is "missing".
        candidate_count (int): How many native FASTA candidates GetOrganelle actually produced.
    """
    status: str
    fasta_path: Path | None
    gfa_path: Path | None
    candidate_count: int

def discover_assembly_artifacts(output_dir: Path, logger: Logger) -> AssemblyArtifacts:
    """
    Discovers GetOrganelle's own native assembly output in output_dir,
    before postprocess_outputs standardizes it to a fixed filename
    pair. GetOrganelle names a completed, circularized assembly with
    ".complete." in its path_sequence filename; anything else it
    produces is a fragmented/draft candidate. Zero candidates means
    the assembly failed outright.

    Args:
        output_dir (Path): The directory GetOrganelle wrote its native output to.
        logger (Logger): The logger to use for this step of the pipeline.

    Returns:
        AssemblyArtifacts: The discovered status and, when available, the
        single best candidate pair (a "complete" one if any exists,
        otherwise the first candidate found).
    """
    # Discovering every native FASTA/GFA candidate GetOrganelle left behind.
    fasta_candidates = sorted(output_dir.glob("*.path_sequence.fasta"))
    gfa_candidates = sorted(output_dir.glob("*.selected_graph.gfa"))

    # Handling the case where GetOrganelle produced nothing at all.
    if not fasta_candidates:
        logger.warning(f"No native assembly FASTA found in {output_dir}; assembly likely failed.")
        return AssemblyArtifacts(status = "missing", fasta_path = None, gfa_path = None, candidate_count = 0)

    # Preferring a candidate GetOrganelle itself marked complete/circular.
    complete_candidates = [path for path in fasta_candidates if ".complete." in path.name]
    chosen_fasta = complete_candidates[0] if complete_candidates else fasta_candidates[0]
    status = "circular" if complete_candidates else "fragmented"

    # Pairing the chosen FASTA with its matching GFA, by shared prefix
    # (GetOrganelle names a candidate's FASTA and GFA with the same
    # stem up to ".path_sequence."/.selected_graph.").
    stem = chosen_fasta.name.split(".path_sequence.fasta")[0]
    matching_gfa = next((path for path in gfa_candidates if path.name.startswith(stem)), None)
    chosen_gfa = matching_gfa or (gfa_candidates[0] if gfa_candidates else None)

    if status == "fragmented":
        logger.warning(f"No complete/circular GetOrganelle candidate in {output_dir}; using the best fragmented candidate ({chosen_fasta.name}).")
    else:
        logger.info(f"Found a complete GetOrganelle assembly in {output_dir} ({chosen_fasta.name}).")

    return AssemblyArtifacts(status = status, fasta_path = chosen_fasta, gfa_path = chosen_gfa,
                              candidate_count = len(fasta_candidates))

@dataclass(frozen = True, slots = True)
class RecruitedReadsQCRecord:
    """
    The outcome of recruited-reads QC for one GetOrganelle attempt
    (#413) -- the structured artifact contract a caller reads instead
    of re-deriving read counts/pairing/provenance from log lines.

    Attributes:
        sample_id (str): The sample this record is for.
        source_output_dir (str): The GetOrganelle attempt's own output_dir this record belongs to (already attempt-specific, so repeated rescue attempts never collide).
        r1_path (str | None): The recruited-reads R1 path, or None if it was never found.
        r2_path (str | None): The recruited-reads R2 path, or None if it was never found.
        r1_read_count (int | None): R1's record count, or None if the file was empty/malformed/absent.
        r2_read_count (int | None): R2's record count, or None if the file was empty/malformed/absent.
        paired (bool | None): Whether R1 and R2 have matching read counts, or None when either count is unavailable.
        qc_ran (bool): Whether FastQC actually completed successfully on these reads.
        fastqc_dir (str | None): The FastQC output directory, or None when qc_ran is False.
        reason (str): A short, human-readable explanation of this outcome.
    """
    sample_id: str
    source_output_dir: str
    r1_path: str | None
    r2_path: str | None
    r1_read_count: int | None
    r2_read_count: int | None
    paired: bool | None
    qc_ran: bool
    fastqc_dir: str | None
    reason: str

def _count_fastq_records(path: Path) -> int | None:
    """
    Counts FASTQ records in a plain-text FASTQ file (4 lines per
    record, the universal FASTQ convention). Returns None rather than
    a count when the file's line count isn't a multiple of 4 -- that
    signals a malformed/truncated FASTQ, not a valid record count.

    Args:
        path (Path): The FASTQ file to count.

    Returns:
        int | None: The number of records, or None if the file is malformed.
    """
    with path.open("r", encoding = "utf-8", errors = "replace") as handle:
        line_count = sum(1 for _ in handle)
    if line_count % 4 != 0:
        return None
    return line_count // 4

def build_runner_from_seed_result(seed_result: SeedSelectionResult,
                                   working_dir: Path,
                                   output_dir: Path,
                                   sample: Sample,
                                   organelle_type: str = "animal_mt",
                                   tool_options: dict | None = None,
                                   threads: int = 4,
                                   logger: Logger | None = None,
                                   metrics_dir: Path | None = None) -> GetOrganelleRunner:
    """
    Configures a GetOrganelleRunner from a seed-selection result: a
    resolved seed FASTA is wired in as the -s seed file; a fallback
    result (no seed found) is left unseeded, so GetOrganelle falls
    back to its own default database for organelle_type.

    Args:
        seed_result (SeedSelectionResult): The outcome of seed_selector.select_seed() for this sample.
        working_dir (Path): The working directory for the command.
        output_dir (Path): The directory to save the output files to.
        sample (Sample): The sample object that is being worked with.
        organelle_type (str): The GetOrganelle organelle type to assemble. Defaults to "animal_mt".
        tool_options (dict | None): A base dictionary of options to pass to the tool; a resolved seed is merged into a copy of this, never mutating the caller's dict.
        threads (int): The number of threads for this process to use.
        logger (Logger | None): The logger to use for this step of the pipeline.
        metrics_dir (Path | None): Where to persist this run's compute metrics as JSON; None skips persistence. See BaseTool.

    Returns:
        GetOrganelleRunner: A runner configured with the seed result wired in.
    """
    # Copying the base tool_options so the caller's own dict is never
    # mutated by this function.
    resolved_options = dict(tool_options or {})

    # Wiring in the resolved seed FASTA only for a genuine (non-fallback) result.
    if not seed_result.fallback_used and seed_result.fasta_path is not None:
        resolved_options["seed_file"] = str(seed_result.fasta_path)

    return GetOrganelleRunner(
        working_dir = working_dir,
        output_dir = output_dir,
        sample = sample,
        organelle_type = organelle_type,
        tool_options = resolved_options,
        threads = threads,
        logger = logger,
        metrics_dir = metrics_dir,
    )
