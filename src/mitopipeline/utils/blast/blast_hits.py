"""
blast_hits.py

Parse, rank, and select tabular BLAST hits for a sample.
"""

# Imports
import csv
from dataclasses import asdict, replace
from logging import Logger
from pathlib import Path
from mitopipeline.utils.blast.models import BlastHit

# The tabular columns this module expects, matching
# api.blast.BLAST_OUTPUT_COLUMNS exactly -- a BLAST run that didn't
# request this full column list will fail the column-count check
# below immediately, rather than silently misaligning values.
BLAST_COLUMNS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore",
    "qcovs", "staxids", "sscinames", "stitle",
]

INTEGER_COLUMNS = {"length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send"}
FLOAT_COLUMNS = {"pident", "evalue", "bitscore", "qcovs"}

def parse_blast_tsv(blast_path: str | Path, logger: Logger | None = None) -> list[BlastHit]:
    """
    Parses tabular BLAST output into typed BlastHit records.

    Blank lines are skipped (BLAST itself never emits them, but a
    hand-edited or concatenated file might). A row with the wrong
    column count, or a numeric column that doesn't parse, is a
    malformed file -- raised immediately with the offending line
    number, rather than silently skipped, since a corrupt BLAST run
    should never be mistaken for a clean one with fewer hits.

    Args:
        blast_path (str | Path): Path to tabular BLAST output.
        logger (Logger | None): Logger to use. Defaults to None.

    Returns:
        list[BlastHit]: The parsed hits, in file order.

    Raises:
        FileNotFoundError: If blast_path does not exist.
        ValueError: If blast_path is not a file, or a row is malformed.
    """
    # Normalizing the path.
    blast_path = Path(blast_path)

    # Verifying that the path exists and is a file.
    if not blast_path.exists():
        if logger is not None: logger.error(f"File {blast_path} not found.")
        raise FileNotFoundError(f"File {blast_path} not found.")
    if not blast_path.is_file():
        if logger is not None: logger.error(f"File {blast_path} is not a file.")
        raise ValueError(f"File {blast_path} is not a file.")

    # Parsing each tabular row into a BlastHit.
    hits: list[BlastHit] = []
    with blast_path.open("r", encoding = "utf-8", errors = "replace", newline = "") as handle:
        reader = csv.reader(handle, delimiter = "\t")
        for line_number, row in enumerate(reader, start = 1):
            # Skipping genuinely blank lines.
            if not row or all(value.strip() == "" for value in row):
                continue

            # Rejecting a row with the wrong number of columns.
            if len(row) != len(BLAST_COLUMNS):
                message = f"BLAST row {line_number} contains {len(row)} columns, expected {len(BLAST_COLUMNS)}."
                if logger is not None: logger.error(message)
                raise ValueError(message)

            # Mapping columns by name, then coercing the numeric ones.
            fields = dict(zip(BLAST_COLUMNS, row, strict = True))
            try:
                for column in INTEGER_COLUMNS:
                    fields[column] = int(fields[column])
                for column in FLOAT_COLUMNS:
                    fields[column] = float(fields[column])
            except ValueError as error:
                message = f"Error parsing BLAST row {line_number}: {error}"
                if logger is not None: logger.error(message)
                raise ValueError(message) from error

            hits.append(BlastHit(**fields))

    if logger is not None: logger.info(f"Parsed {len(hits)} BLAST hits from {blast_path}.")
    return hits

def _ranking_key(hit: BlastHit) -> tuple:
    """
    Returns a deterministic BLAST ranking key. Lower tuple values
    represent a stronger match.

    Args:
        hit (BlastHit): The hit to compute a ranking key for.

    Returns:
        tuple: The ranking key.
    """
    return (hit.evalue, -hit.bitscore, -hit.qcovs, -hit.pident, -hit.length, hit.sseqid)

def rank_blast_hits(hits: list[BlastHit]) -> list[BlastHit]:
    """
    Ranks BLAST hits, keeping only the single best hit per subject
    (sseqid).

    Args:
        hits (list[BlastHit]): The hits to rank.

    Returns:
        list[BlastHit]: The best hit per subject, sorted strongest first.
    """
    # Keeping only the strongest hit seen so far per subject.
    best_by_subject: dict[str, BlastHit] = {}
    for hit in hits:
        current_best = best_by_subject.get(hit.sseqid)
        if current_best is None or _ranking_key(hit) < _ranking_key(current_best):
            best_by_subject[hit.sseqid] = hit

    return sorted(best_by_subject.values(), key = _ranking_key)

def select_top_blast_hits(hits: list[BlastHit], maximum_hits: int = 6, logger: Logger | None = None) -> list[BlastHit]:
    """
    Selects the highest-ranked unique BLAST hits, numbering each one's
    rank starting at 1.

    Args:
        hits (list[BlastHit]): The hits to select from.
        maximum_hits (int): The maximum number of hits to select. Defaults to 6.
        logger (Logger | None): Logger to use. Defaults to None.

    Returns:
        list[BlastHit]: The selected hits, each with rank set.

    Raises:
        ValueError: If maximum_hits is not greater than 0.
    """
    # Verifying that maximum_hits is a sensible value.
    if maximum_hits <= 0:
        if logger is not None: logger.error("maximum_hits must be greater than 0.")
        raise ValueError("maximum_hits must be greater than 0.")

    ranked = rank_blast_hits(hits)
    selected = [replace(hit, rank = rank) for rank, hit in enumerate(ranked[:maximum_hits], start = 1)]

    if logger is not None and len(selected) < maximum_hits:
        logger.warning(f"Only {len(selected)} unique BLAST hits were available; {maximum_hits} were requested.")

    return selected

def write_top_blast_hits(hits: list[BlastHit], output_path: str | Path, sample_id: str, logger: Logger | None = None) -> None:
    """
    Writes selected BLAST hits to a TSV file, one row per hit.

    Args:
        hits (list[BlastHit]): The hits to write.
        output_path (str | Path): The path to write to.
        sample_id (str): The sample these hits belong to.
        logger (Logger | None): Logger to use. Defaults to None.

    Returns:
        None
    """
    # Normalizing the output path and ensuring its parent exists.
    output_path = Path(output_path)
    output_path.parent.mkdir(parents = True, exist_ok = True)

    # Writing one row per hit, sample_id and rank first.
    output_columns = ["sample_id", "rank", *BLAST_COLUMNS]
    with output_path.open("w", encoding = "utf-8", newline = "") as handle:
        writer = csv.DictWriter(handle, fieldnames = output_columns, delimiter = "\t")
        writer.writeheader()
        for hit in hits:
            row = {"sample_id": sample_id, **asdict(hit)}
            writer.writerow(row)

def read_top_blast_hits(path: str | Path, sample_id: str | None = None, logger: Logger | None = None) -> list[BlastHit]:
    """
    Reads back a TSV write_top_blast_hits() wrote -- the inverse of
    that function, NOT parse_blast_tsv(), since this format carries
    two extra leading columns (sample_id, rank) that raw blastn
    tabular output never has.

    Args:
        path (str | Path): The top-hit TSV to read.
        sample_id (str | None): If given, only rows for this sample are returned (the TSV may carry multiple samples' hits).
        logger (Logger | None): Logger to use. Defaults to None.

    Returns:
        list[BlastHit]: The hits read, already ranked (write_top_blast_hits itself always writes them in rank order).

    Raises:
        FileNotFoundError: If path does not exist.
    """
    path = Path(path)
    if not path.is_file():
        if logger is not None: logger.error(f"File {path} not found.")
        raise FileNotFoundError(f"File {path} not found.")

    hits: list[BlastHit] = []
    with path.open("r", encoding = "utf-8", newline = "") as handle:
        for row in csv.DictReader(handle, delimiter = "\t"):
            if sample_id is not None and row.get("sample_id") != sample_id:
                continue
            fields = {column: row[column] for column in BLAST_COLUMNS}
            for column in INTEGER_COLUMNS:
                fields[column] = int(fields[column])
            for column in FLOAT_COLUMNS:
                fields[column] = float(fields[column])
            fields["rank"] = int(row["rank"]) if row.get("rank") else None
            hits.append(BlastHit(**fields))

    if logger is not None: logger.info(f"Read {len(hits)} top BLAST hit(s) from {path}.")
    return hits

    if logger is not None: logger.info(f"Wrote {len(hits)} BLAST hits to {output_path}.")
