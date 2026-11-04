"""
getorganelle_database.py

Makes sure GetOrganelle's own default databases for an organelle type
(the LabelDatabase BLAST index slim_graph.py filters assembly graphs
with, and the SeedDatabase bowtie2 index an unseeded run recruits reads
with) exist AND are usable by the BLAST/bowtie2 binaries of the
environment the pipeline is about to assemble in.

Why this exists: GetOrganelle keeps these databases in one shared
directory (~/.GetOrganelle, or $GETORG_PATH) and only ever checks that
the files are present. A BLAST index built there by one BLAST version
(another environment, an older pipeline version, a system install) can
be unreadable by the BLAST of the current environment ("BLAST Database
error: Not a valid version 4 database"). GetOrganelle then logs
"slimming ... failed", quietly carries on with an UNFILTERED seed graph,
recruits most of the library, and every sample ends in "No valid
assembly graph found!" -- a whole batch lost to a one-line makeblastdb.
This module runs that check up front, once per job, and rebuilds the
index with the environment's own tools when it fails.
"""

# Imports
from __future__ import annotations
from dataclasses import dataclass, field
import logging
import os
from pathlib import Path
import shutil
import subprocess
from typing import Callable, Sequence

# The subprocess runner signature this module depends on, injectable
# for tests: given argv, returns a CompletedProcess with returncode,
# stdout and stderr.
Runner = Callable[[Sequence[str]], subprocess.CompletedProcess]

def default_config_dir() -> Path:
    """
    GetOrganelle's own database/config directory: $GETORG_PATH when set
    (GetOrganelle honours the same variable), else ~/.GetOrganelle.

    Returns:
        Path: The directory GetOrganelle reads its default databases from.
    """
    configured = os.environ.get("GETORG_PATH")
    return Path(configured).expanduser() if configured else Path.home() / ".GetOrganelle"

def _run(argv: Sequence[str]) -> subprocess.CompletedProcess:
    """The real runner: captures output, never raises on a non-zero exit."""
    return subprocess.run(list(argv), capture_output = True, text = True, check = False)

@dataclass
class GetOrganelleDatabaseResult:
    """
    What ensure() found and did.

    Attributes:
        success (bool): Whether both databases are present and usable by this environment's tools.
        config_dir (Path): The GetOrganelle config directory that was checked.
        organelle_type (str): The organelle type (e.g. "animal_mt") checked.
        blast_version (str | None): The environment's BLAST version string, if it could be read.
        downloaded (bool): Whether get_organelle_config.py had to fetch the databases.
        label_index_rebuilt (bool): Whether the LabelDatabase BLAST index had to be rebuilt.
        seed_index_rebuilt (bool): Whether the SeedDatabase bowtie2 index had to be rebuilt.
        reason (str): A one-line human explanation of the outcome.
        actions (list[str]): The commands run, in order, for the record.
    """
    success: bool
    config_dir: Path
    organelle_type: str
    blast_version: str | None = None
    downloaded: bool = False
    label_index_rebuilt: bool = False
    seed_index_rebuilt: bool = False
    reason: str = ""
    actions: list[str] = field(default_factory = list)

class GetOrganelleDatabase:
    """
    Verifies (and repairs) GetOrganelle's default databases for one
    organelle type against the current environment's own tools.
    """

    # bowtie2 writes six index files; the first forward one is enough to
    # tell whether an index exists at all (large indexes use .bt2l).
    _BOWTIE2_FIRST_INDEX = (".1.bt2", ".1.bt2l")

    def __init__(self, organelle_type: str, logger: logging.Logger,
                 config_dir: Path | None = None, runner: Runner = _run,
                 which: Callable[[str], str | None] = shutil.which) -> None:
        """
        Args:
            organelle_type (str): The GetOrganelle organelle type (e.g. "animal_mt").
            logger (logging.Logger): The logger to use.
            config_dir (Path | None): GetOrganelle's config directory. Defaults to default_config_dir().
            runner (Runner): How commands are run (injectable for tests).
            which (Callable): How executables are located on PATH (injectable for tests).
        """
        self.organelle_type = organelle_type
        self.logger = logger
        self.config_dir = Path(config_dir) if config_dir is not None else default_config_dir()
        self.runner = runner
        self.which = which

    # --- Paths ------------------------------------------------------------

    @property
    def label_prefix(self) -> Path:
        """The LabelDatabase path without extension (what -db / makeblastdb -out take)."""
        return self.config_dir / "LabelDatabase" / self.organelle_type

    @property
    def label_fasta(self) -> Path:
        return self.label_prefix.with_suffix(".fasta")

    @property
    def seed_prefix(self) -> Path:
        """The SeedDatabase path without extension (what bowtie2 -x / bowtie2-build take)."""
        return self.config_dir / "SeedDatabase" / self.organelle_type

    @property
    def seed_fasta(self) -> Path:
        return self.seed_prefix.with_suffix(".fasta")

    # --- Steps ------------------------------------------------------------

    def _call(self, argv: Sequence[str], actions: list[str]) -> subprocess.CompletedProcess:
        """Runs one command, recording it, and logs its stderr on failure."""
        actions.append(" ".join(str(a) for a in argv))
        self.logger.info(f"Running: {actions[-1]}")
        completed = self.runner(argv)
        if completed.returncode != 0:
            detail = (completed.stderr or completed.stdout or "").strip()
            self.logger.warning(f"Command exited {completed.returncode}: {argv[0]}" + (f" -- {detail}" if detail else ""))
        return completed

    def blast_version(self) -> str | None:
        """The environment's blastn version line, or None if blastn is absent."""
        if self.which("blastn") is None:
            return None
        completed = self.runner(["blastn", "-version"])
        first_line = (completed.stdout or "").strip().splitlines()
        return first_line[0] if first_line else None

    def label_index_readable(self) -> bool:
        """
        Whether this environment's BLAST can open the LabelDatabase
        index -- the exact check that fails with "Not a valid version 4
        database" when the index was built by a different BLAST.
        """
        completed = self.runner(["blastdbcmd", "-db", str(self.label_prefix), "-info"])
        return completed.returncode == 0

    def seed_index_present(self) -> bool:
        """Whether a bowtie2 index exists for the SeedDatabase FASTA."""
        return any(self.seed_prefix.with_name(self.seed_prefix.name + suffix).is_file()
                   for suffix in self._BOWTIE2_FIRST_INDEX)

    def ensure(self) -> GetOrganelleDatabaseResult:
        """
        Ensures both databases exist and are usable, downloading or
        rebuilding as needed. Never raises: the result says what
        happened, and success False means assembly cannot work in this
        environment until the stated reason is fixed.

        Returns:
            GetOrganelleDatabaseResult: The outcome.
        """
        result = GetOrganelleDatabaseResult(success = False, config_dir = self.config_dir,
                                            organelle_type = self.organelle_type)
        result.blast_version = self.blast_version()

        # The tools this check depends on must be on PATH (the assembly
        # conda environment provides all of them with getorganelle).
        missing = [tool for tool in ("blastn", "blastdbcmd", "makeblastdb", "bowtie2-build", "get_organelle_config.py")
                   if self.which(tool) is None]
        if missing:
            result.reason = f"required tool(s) not on PATH: {', '.join(missing)}"
            self.logger.error(f"GetOrganelle database check cannot run: {result.reason}.")
            return result

        # 1. Fetch the databases if GetOrganelle has never had them here.
        if not self.label_fasta.is_file() or not self.seed_fasta.is_file():
            self.logger.info(f"GetOrganelle {self.organelle_type} databases not found under {self.config_dir}; downloading.")
            completed = self._call(["get_organelle_config.py", "--add", self.organelle_type], result.actions)
            result.downloaded = True
            if completed.returncode != 0 or not self.label_fasta.is_file() or not self.seed_fasta.is_file():
                result.reason = f"get_organelle_config.py --add {self.organelle_type} did not produce {self.label_fasta} and {self.seed_fasta}"
                self.logger.error(result.reason)
                return result

        # 2. The LabelDatabase BLAST index must be readable by THIS
        #    environment's BLAST; rebuild it with this environment's
        #    makeblastdb if not.
        if not self.label_index_readable():
            self.logger.warning(f"LabelDatabase index {self.label_prefix} is missing or unreadable by this environment's BLAST "
                                f"({result.blast_version or 'unknown version'}); rebuilding it with makeblastdb.")
            completed = self._call(["makeblastdb", "-in", str(self.label_fasta), "-dbtype", "nucl", "-out", str(self.label_prefix)],
                                   result.actions)
            result.label_index_rebuilt = True
            if completed.returncode != 0 or not self.label_index_readable():
                result.reason = f"could not build a LabelDatabase BLAST index readable by this environment at {self.label_prefix}"
                self.logger.error(result.reason)
                return result

        # 3. The SeedDatabase bowtie2 index just has to exist (bowtie2's
        #    index format is stable across versions).
        if not self.seed_index_present():
            self.logger.warning(f"SeedDatabase bowtie2 index {self.seed_prefix} is missing; building it.")
            completed = self._call(["bowtie2-build", str(self.seed_fasta), str(self.seed_prefix)], result.actions)
            result.seed_index_rebuilt = True
            if completed.returncode != 0 or not self.seed_index_present():
                result.reason = f"could not build a SeedDatabase bowtie2 index at {self.seed_prefix}"
                self.logger.error(result.reason)
                return result

        result.success = True
        rebuilt = [name for name, flag in (("LabelDatabase BLAST index", result.label_index_rebuilt),
                                           ("SeedDatabase bowtie2 index", result.seed_index_rebuilt)) if flag]
        result.reason = ("databases ready" + (f" (downloaded)" if result.downloaded else "")
                         + (f"; rebuilt: {', '.join(rebuilt)}" if rebuilt else ""))
        self.logger.info(f"GetOrganelle {self.organelle_type} databases under {self.config_dir}: {result.reason}.")
        return result
