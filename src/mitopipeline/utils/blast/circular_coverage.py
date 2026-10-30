"""
circular_coverage.py

Corrects BLAST query coverage for a circular assembly (the biological validation gate's own
"blast_identity" check): linear BLAST can only report a contiguous
match, so a circular genome whose rotation differs from its
reference's own rotation shows an origin-spanning match as two
separate, shorter HSPs -- read as low query coverage even when the
whole genome is truly a strong match. The standard fix is to BLAST
a doubled copy of the query (sequence + sequence) so a wraparound
match becomes one unbroken HSP, then measure coverage against the
real, undoubled genome length.
"""

# Imports
from logging import Logger
from pathlib import Path
from mitopipeline.utils.blast.blast_hits import rank_blast_hits
from mitopipeline.utils.blast.models import BlastHit

def double_fasta_sequence(fasta_path: str | Path, output_path: str | Path) -> int:
    """
    Writes a doubled-sequence FASTA (the real sequence concatenated
    to itself, same header) -- a single-record FASTA is expected,
    matching contig_selection's own "one representative sequence"
    convention.

    Args:
        fasta_path (str | Path): The real, undoubled genome FASTA.
        output_path (str | Path): Where to write the doubled FASTA.

    Returns:
        int: The real, undoubled genome length.

    Raises:
        FileNotFoundError: If fasta_path doesn't exist.
        ValueError: If fasta_path has no records.
    """
    fasta_path = Path(fasta_path)
    if not fasta_path.is_file():
        raise FileNotFoundError(f"Genome FASTA not found: {fasta_path}.")

    header, sequence_lines = None, []
    with fasta_path.open("r", encoding = "utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    break
                header = line
            elif header is not None:
                sequence_lines.append(line.strip())

    if header is None:
        raise ValueError(f"Genome FASTA contains no records: {fasta_path}.")

    sequence = "".join(sequence_lines)
    if not sequence:
        raise ValueError(f"Genome FASTA record has no bases: {fasta_path}.")

    output_path = Path(output_path)
    output_path.parent.mkdir(parents = True, exist_ok = True)
    output_path.write_text(f"{header}\n{sequence}{sequence}\n", encoding = "utf-8")

    return len(sequence)

def compute_circular_query_coverage(doubled_query_hits: list[BlastHit], genome_length: int) -> tuple[float, int] | None:
    """
    Computes the real, circularity-corrected query coverage from a
    BLAST run against the doubled query -- the strongest hit's own
    alignment length against the real, undoubled genome length,
    capped at 100% (a repetitive region could otherwise double-count
    across the duplicated query).

    Args:
        doubled_query_hits (list[BlastHit]): BLAST hits from running the doubled query against the reference database.
        genome_length (int): The real, undoubled genome length.

    Returns:
        tuple[float, int] | None: (corrected coverage percent, alignment length), or None if there are no hits.

    Raises:
        ValueError: If genome_length is not greater than zero.
    """
    if genome_length <= 0:
        raise ValueError(f"genome_length must be greater than zero, got {genome_length}.")

    ranked = rank_blast_hits(doubled_query_hits)
    if not ranked:
        return None

    top_hit = ranked[0]
    coverage = min(100.0, top_hit.length / genome_length * 100.0)
    return coverage, top_hit.length
