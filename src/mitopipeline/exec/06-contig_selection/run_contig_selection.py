"""
run_contig_selection.py

Execution layer that selects one sample's final mitochondrial contig
across one or more GetOrganelle assembly attempts: reads each
attempt's recorded status (the `.assembly_status.json` run_getorganelle_
assembly.py writes) and its standardized output (`data.fasta`/
`graph.gfa`, same directory), then records the winning selection as a
JSON manifest. The winning candidate's FASTA is also resolved down to
a single record here (the longest, when the native file holds more
than one disjoint scaffold) -- every later single-sequence consumer
(phylogeny dataset construction, MITOS2 annotation) gets exactly one
real sequence to work with, never a multi-record file it has to
reject.

Doesn't run GetOrganelle or rescue itself -- those are earlier stages'
jobs. This step exists so later stages (BLAST,
annotation, reporting) read one deterministic, already-decided answer
instead of re-deriving "which attempt won" themselves.
"""

# Imports
import argparse
import json
from pathlib import Path
from Bio import SeqIO
from mitopipeline.logging.logger_factory import make_logger
from mitopipeline.utils.contig_selection.models import ContigCandidate
from mitopipeline.utils.contig_selection.selector import select_mitochondrial_contig

def parse_args() -> argparse.Namespace:
    """
    Parses the command line arguments for this script.

    Returns:
        argparse.Namespace: The parsed command line arguments.
    """
    parser = argparse.ArgumentParser(description = "Selects one sample's final mitochondrial contig across its assembly attempts.")

    parser.add_argument("--sample-id", type = str, required = True, help = "The sample id.")
    parser.add_argument("--attempt-dir", type = Path, required = True, action = "append",
                        help = "An assembly attempt's output directory (data.fasta/graph.gfa/<sample_id>.assembly_status.json). "
                               "Repeat in attempt order (first = primary, then each rescue attempt); order only matters for tie-breaking.")
    parser.add_argument("--output-dir", type = Path, required = True, help = "The directory to write this step's own output (the selection JSON) to.")
    parser.add_argument("--log-file", type = Path, required = True, help = "The path to the log file.")
    parser.add_argument("--global-log-file", type = Path, required = True, help = "The path to the global log file.")

    return parser.parse_args()

def _load_candidate(attempt: int, sample_id: str, attempt_dir: Path, logger) -> ContigCandidate:
    """
    Reads one assembly attempt's recorded status and standardized
    output paths into a ContigCandidate.

    A directory with no status file at all (the assembly step never
    ran, or crashed before writing one) is treated as "missing" --
    the same outcome as an assembly that genuinely produced nothing --
    rather than raising, since a never-attempted rescue slot is not a
    distinct failure mode this step needs to react to differently.

    Args:
        attempt (int): The attempt number (0 = primary, 1+ = rescue attempts).
        sample_id (str): The sample id.
        attempt_dir (Path): The attempt's own output directory.
        logger: The logger to use.

    Returns:
        ContigCandidate: The candidate read from this attempt.
    """
    status_path = attempt_dir / f"{sample_id}.assembly_status.json"
    try:
        status = json.loads(status_path.read_text(encoding = "utf-8"))["status"]
    except (OSError, json.JSONDecodeError, KeyError) as e:
        logger.warning(f"{{{sample_id}}} attempt {attempt}: no readable status at {status_path} ({e!r}); treating as missing.")
        return ContigCandidate(attempt = attempt, status = "missing", fasta_path = None, gfa_path = None)

    if status not in ("circular", "fragmented"):
        return ContigCandidate(attempt = attempt, status = "missing", fasta_path = None, gfa_path = None)

    fasta_path = attempt_dir / "data.fasta"
    gfa_path = attempt_dir / "graph.gfa"
    if not fasta_path.is_file():
        logger.warning(f"{{{sample_id}}} attempt {attempt}: status is '{status}' but {fasta_path} is missing; treating as missing.")
        return ContigCandidate(attempt = attempt, status = "missing", fasta_path = None, gfa_path = None)

    return ContigCandidate(attempt = attempt, status = status, fasta_path = fasta_path,
                           gfa_path = gfa_path if gfa_path.is_file() else None)

def _write_resolved_fasta(source_fasta: Path, output_path: Path) -> None:
    """
    Writes the single longest record from source_fasta to output_path.

    A normal single-record candidate (the common case) round-trips
    through unchanged. A "fragmented" candidate can hold multiple
    disjoint scaffold records when GetOrganelle couldn't link them
    into one contig (confirmed live: a 5-scaffold, 374-854bp-each
    primary assembly) -- every downstream single-sequence consumer
    (phylogeny dataset construction, MITOS2 annotation) expects
    exactly one record and rejects the whole file otherwise. This
    resolves that ambiguity here, once, the same way this step
    already resolves which ATTEMPT won -- so nothing downstream has
    to special-case a multi-record FASTA itself.

    Args:
        source_fasta (Path): The selected candidate's native FASTA (one or more records).
        output_path (Path): Where to write the single resolved record.

    Returns:
        None
    """
    records = list(SeqIO.parse(source_fasta, "fasta"))
    longest = max(records, key = lambda r: len(r.seq))
    output_path.parent.mkdir(parents = True, exist_ok = True)
    SeqIO.write([longest], output_path, "fasta")

def main() -> int:
    """
    Selects and records one sample's final mitochondrial contig.

    Returns:
        int: The return code for the process.
    """
    args = parse_args()

    logger = make_logger(
        name = "contig_selection",
        log_file_path = args.log_file,
        global_log_file_path = args.global_log_file,
    )
    context = f"{{{args.sample_id}}}"

    candidates = [_load_candidate(i, args.sample_id, d, logger)
                  for i, d in enumerate(args.attempt_dir)]

    selection = select_mitochondrial_contig(candidates, logger = logger)

    args.output_dir.mkdir(parents = True, exist_ok = True)
    output_path = args.output_dir / f"{args.sample_id}.contig_selection.json"

    if selection is None:
        logger.warning(f"{context} No circular or fragmented candidate among {len(candidates)} attempt(s); nothing selected.")
        with output_path.open("w", encoding = "utf-8") as handle:
            json.dump({"sample_id": args.sample_id, "selected": False,
                      "attempts_considered": len(candidates)}, handle, indent = 4)
        return 1

    logger.info(f"{context} Selected attempt {selection.source_attempt} "
                f"({selection.status}, {selection.length_bp}bp): {selection.reason}")
    resolved_fasta = args.output_dir / f"{args.sample_id}.selected.fasta"
    _write_resolved_fasta(selection.fasta_path, resolved_fasta)
    with output_path.open("w", encoding = "utf-8") as handle:
        json.dump({
            "sample_id": args.sample_id,
            "selected": True,
            "status": selection.status,
            "fasta_path": str(resolved_fasta),
            "gfa_path": str(selection.gfa_path) if selection.gfa_path else None,
            "source_attempt": selection.source_attempt,
            "length_bp": selection.length_bp,
            "reason": selection.reason,
            "attempts_considered": len(candidates),
        }, handle, indent = 4)
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
