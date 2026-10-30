"""
table2asn_validation.py

Prepares table2asn's own real basename-matched (.fsa/.tbl) input
pair and interprets its real validation output (.val/.stats/.dr)
into a typed PASS/FAIL result -- the real local validation step
NCBI's own documentation calls for before any GenBank submission,
run entirely offline against files already on disk. Never submits
or contacts NCBI itself; table2asn is a local validator.
"""

# Imports
from __future__ import annotations
from dataclasses import asdict, dataclass
import re
from logging import Logger
from pathlib import Path

# table2asn's own real message severities, as they appear verbatim
# in its .val/.stats/.dr files and stdout/stderr (confirmed against
# v0.1.0's own real prototype of this exact check).
SEVERITY_PATTERN = re.compile(r"\b(REJECT|FATAL|ERROR|WARNING|INFO)\b", flags = re.IGNORECASE)


@dataclass(frozen = True, slots = True)
class Table2AsnValidationResult:
    """
    One sample's real table2asn validation outcome.

    Attributes:
        sample_id (str): The sample id.
        status (str): "PASS" or "FAIL".
        return_code (int): table2asn's own real process exit code.
        sqn_file (str): The real .sqn submission file path, if produced.
        validation_file (str): The real .val file path.
        stats_file (str): The real .stats file path.
        discrepancy_file (str): The real .dr file path.
        genbank_flatfile (str): The real .gbf GenBank flatfile path, if produced.
        reject_count (int): Real REJECT-level message count.
        fatal_count (int): Real FATAL-level message count.
        error_count (int): Real ERROR-level message count.
        warning_count (int): Real WARNING-level message count (reviewable, never a blocker on its own).
        info_count (int): Real INFO-level message count.
        blockers (str): A semicolon-joined, human-readable list of the real reasons status is FAIL, empty when PASS.
    """
    sample_id: str
    status: str
    return_code: int
    sqn_file: str
    validation_file: str
    stats_file: str
    discrepancy_file: str
    genbank_flatfile: str
    reject_count: int
    fatal_count: int
    error_count: int
    warning_count: int
    info_count: int
    blockers: str

    def to_dict(self) -> dict:
        """Serializes this result to a JSON-safe dictionary."""
        return asdict(self)


def tag_fasta_with_ncbi_defline(source_fasta: Path, destination_fasta: Path, sample_id: str, organism: str,
                                 topology: str, genome_location: str) -> None:
    """
    Copies a genome FASTA while rewriting its defline with the real
    NCBI source modifiers table2asn's own documentation requires
    ([organism=...] [topology=...] [location=...]) -- never invented,
    read back from this sample's own already-known organism/topology.

    Args:
        source_fasta (Path): The real assembled genome FASTA.
        destination_fasta (Path): Where to write the tagged copy (table2asn's own required .fsa basename).
        sample_id (str): The sample id, used as the real sequence id.
        organism (str): The real organism name.
        topology (str): "circular" or "linear".
        genome_location (str): The real genome location (e.g. "mitochondrion").

    Returns:
        None

    Raises:
        ValueError: If source_fasta has no real sequence.
    """
    lines = Path(source_fasta).read_text(encoding = "utf-8").splitlines()
    sequence_lines = [line.strip() for line in lines if line.strip() and not line.startswith(">")]
    if not sequence_lines:
        raise ValueError(f"No sequence was found in {source_fasta}.")

    defline = f">{sample_id} [organism={organism}] [topology={topology}] [location={genome_location}]"
    Path(destination_fasta).write_text(defline + "\n" + "\n".join(sequence_lines) + "\n", encoding = "utf-8")


def parse_table2asn_messages(paths: list[Path]) -> dict[str, int]:
    """
    Counts table2asn's own real messages by severity, across every
    real output file it produced (and stdout/stderr, when included).
    A missing file (table2asn didn't produce it, e.g. after a hard
    failure) is skipped, not an error -- its absence is itself
    already reflected in the caller's own blockers list.

    Args:
        paths (list[Path]): The real files to scan.

    Returns:
        dict[str, int]: {"reject", "fatal", "error", "warning", "info"} real counts.
    """
    counts = {"reject": 0, "fatal": 0, "error": 0, "warning": 0, "info": 0}
    for path in paths:
        path = Path(path)
        if not path.exists():
            continue
        text = path.read_text(encoding = "utf-8", errors = "replace")
        for match in SEVERITY_PATTERN.finditer(text):
            counts[match.group(1).lower()] += 1
    return counts


def build_validation_result(sample_id: str, return_code: int, output_basename: Path,
                             stdout_path: Path, stderr_path: Path) -> Table2AsnValidationResult:
    """
    Builds the real PASS/FAIL validation result for one sample from
    table2asn's own real output files. FAIL requires a concrete,
    named reason (a nonzero exit, a real reject/fatal/error message,
    or a missing .sqn) -- never a bare "something went wrong."
    WARNING-level messages never fail validation on their own; they
    are counted and reported for manual review, matching table2asn's
    own real severity model (this is local validation, not GenBank
    curator review -- see the module docstring).

    Args:
        sample_id (str): The sample id.
        return_code (int): table2asn's own real process exit code.
        output_basename (Path): The real output basename passed to table2asn's own -o flag.
        stdout_path (Path): Where this run's real stdout was captured.
        stderr_path (Path): Where this run's real stderr was captured.

    Returns:
        Table2AsnValidationResult: The real, typed validation outcome.
    """
    expected = {
        "sqn": output_basename.with_suffix(".sqn"), "val": output_basename.with_suffix(".val"),
        "stats": output_basename.with_suffix(".stats"), "dr": output_basename.with_suffix(".dr"),
        "gbf": output_basename.with_suffix(".gbf"),
    }
    counts = parse_table2asn_messages([expected["val"], expected["stats"], expected["dr"], stdout_path, stderr_path])

    blockers: list[str] = []
    if return_code != 0:
        blockers.append(f"table2asn exited with code {return_code}")
    if counts["reject"]:
        blockers.append(f"{counts['reject']} reject-level message(s)")
    if counts["fatal"]:
        blockers.append(f"{counts['fatal']} fatal message(s)")
    if counts["error"]:
        blockers.append(f"{counts['error']} error message(s)")
    if not expected["sqn"].exists():
        blockers.append("No .sqn submission file was generated")

    return Table2AsnValidationResult(
        sample_id = sample_id, status = "PASS" if not blockers else "FAIL", return_code = return_code,
        sqn_file = str(expected["sqn"]), validation_file = str(expected["val"]), stats_file = str(expected["stats"]),
        discrepancy_file = str(expected["dr"]), genbank_flatfile = str(expected["gbf"]),
        reject_count = counts["reject"], fatal_count = counts["fatal"], error_count = counts["error"],
        warning_count = counts["warning"], info_count = counts["info"], blockers = ";".join(blockers),
    )
