"""
combined_fasta.py

Builds the combined, batch-level multi-sequence FASTA (#344/#345) a
run's accepted samples are aligned and compared against each other
on -- one representative sequence per sample, keyed by sample_id, in
deterministic order, with any rejected or missing sample skipped and
recorded rather than silently dropped.
"""

# Imports
from __future__ import annotations
from dataclasses import dataclass
from logging import Logger
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

@dataclass(frozen = True, slots = True)
class CombinedFastaEntry:
    """
    One sample's own candidate for the combined FASTA, as contig
    selection and biological validation already decided it.

    Attributes:
        sample_id (str): The sample id.
        fasta_path (Path | None): The sample's selected assembly FASTA, or None if nothing was selected for it.
        accepted (bool): Whether this sample's assembly was biologically accepted (its own validation didn't FAIL).
    """
    sample_id: str
    fasta_path: Path | None
    accepted: bool

@dataclass(frozen = True, slots = True)
class CombinedFastaResult:
    """
    The outcome of building the combined FASTA.

    Attributes:
        output_path (Path): Where the combined FASTA was written.
        included_sample_ids (tuple[str, ...]): Every sample actually included, in the deterministic (sorted) order they were written.
        excluded_sample_ids (tuple[str, ...]): Every sample considered but excluded, sorted.
        exclusion_reasons (dict[str, str]): Why each excluded sample was excluded, keyed by sample_id.
    """
    output_path: Path
    included_sample_ids: tuple[str, ...]
    excluded_sample_ids: tuple[str, ...]
    exclusion_reasons: dict[str, str]

def build_combined_fasta(entries: list[CombinedFastaEntry], output_path: str | Path,
                          logger: Logger | None = None) -> CombinedFastaResult:
    """
    Builds the combined, batch-level FASTA from every sample entry
    considered, in deterministic (sorted-by-sample_id) order --
    regardless of the order entries was given in, so the combined
    dataset (and everything downstream of it) never depends on
    incidental manifest/filesystem ordering.

    A rejected sample (accepted is False), a sample with no selected
    FASTA at all, or a FASTA that doesn't actually exist/is empty are
    all skipped and recorded, never silently dropped or treated as a
    crash. A sample's FASTA with more than one record (a fragmented
    assembly with multiple scaffolds) contributes only its first
    (longest-ranked, per contig_selection's own ordering) record --
    the combined dataset is one representative sequence per sample.

    Args:
        entries (list[CombinedFastaEntry]): Every sample considered for this batch.
        output_path (str | Path): Where to write the combined FASTA.
        logger (Logger | None): Logger to use. Defaults to None.

    Returns:
        CombinedFastaResult: The outcome, including which samples were included/excluded and why.

    Raises:
        ValueError: If entries contains duplicate sample_ids, or if zero samples end up included.
    """
    output_path = Path(output_path)

    sample_ids = [entry.sample_id for entry in entries]
    if len(sample_ids) != len(set(sample_ids)):
        duplicates = sorted({sid for sid in sample_ids if sample_ids.count(sid) > 1})
        raise ValueError(f"Duplicate sample_id(s) in combined FASTA input: {duplicates}")

    exclusion_reasons: dict[str, str] = {}
    included_records: list[tuple[str, SeqRecord]] = []

    for entry in sorted(entries, key = lambda e: e.sample_id):
        if not entry.accepted:
            exclusion_reasons[entry.sample_id] = "assembly not biologically accepted"
            continue
        if entry.fasta_path is None:
            exclusion_reasons[entry.sample_id] = "no selected assembly FASTA"
            continue
        fasta_path = Path(entry.fasta_path)
        if not fasta_path.is_file() or fasta_path.stat().st_size == 0:
            exclusion_reasons[entry.sample_id] = f"selected FASTA missing or empty: {fasta_path}"
            continue

        records = list(SeqIO.parse(fasta_path, "fasta"))
        if not records:
            exclusion_reasons[entry.sample_id] = f"selected FASTA contains no records: {fasta_path}"
            continue
        if len(records) > 1 and logger is not None:
            logger.warning(f"{{{entry.sample_id}}} {fasta_path} has {len(records)} records; "
                           f"using only the first ({records[0].id}) for the combined dataset.")

        included_records.append((entry.sample_id, records[0]))

    if not included_records:
        raise ValueError("No samples were included in the combined FASTA -- nothing to align or build a tree from.")

    output_path.parent.mkdir(parents = True, exist_ok = True)
    with output_path.open("w", encoding = "utf-8") as handle:
        for sample_id, record in included_records:
            handle.write(f">{sample_id}\n{str(record.seq)}\n")

    included_sample_ids = tuple(sample_id for sample_id, _ in included_records)
    excluded_sample_ids = tuple(sorted(exclusion_reasons))

    if logger is not None:
        logger.info(f"Combined FASTA: {len(included_sample_ids)} sample(s) included, "
                    f"{len(excluded_sample_ids)} excluded. Written to {output_path}.")

    return CombinedFastaResult(
        output_path = output_path,
        included_sample_ids = included_sample_ids,
        excluded_sample_ids = excluded_sample_ids,
        exclusion_reasons = exclusion_reasons,
    )
