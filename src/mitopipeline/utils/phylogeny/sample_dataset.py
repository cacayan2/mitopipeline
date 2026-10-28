"""
sample_dataset.py

Builds one sample's own real phylogenetic dataset: its assembled
genome plus the real sequences for its own top-ranked BLAST reference
hits, fetched live from NCBI Entrez -- the same per-sample design
v0.1.0 used (ported here, not re-derived), so a report can always
carry a real tree for any sample with a usable assembly and at least
one BLAST hit, independent of how many other samples are in the same
batch.
"""

# Imports
import re
from logging import Logger
from pathlib import Path
from typing import Callable
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
from mitopipeline.utils.blast.blast_hits import read_top_blast_hits
from mitopipeline.utils.blast.models import BlastHit

_DATABASE_PREFIXES = {"ref", "gb", "emb", "dbj", "tpg", "tpe", "tpd"}

def normalize_blast_accession(subject_id: str) -> str:
    """
    Extracts a bare nucleotide accession from a BLAST subject
    identifier, which may carry a database-prefix wrapper (e.g.
    "ref|NC_002333.2|") or already be bare, depending on how the
    local BLAST database was built.

    Args:
        subject_id (str): The BLAST hit's own sseqid value.

    Returns:
        str: The bare accession.

    Raises:
        ValueError: If no accession can be extracted.
    """
    subject_id = subject_id.strip()
    if not subject_id:
        raise ValueError("BLAST subject identifier cannot be empty.")

    parts = [part.strip() for part in subject_id.split("|") if part.strip()]
    if len(parts) >= 2 and parts[0].lower() in _DATABASE_PREFIXES:
        accession = parts[1]
    else:
        accession = parts[0]

    if not accession:
        raise ValueError(f"Could not extract an accession from BLAST subject identifier: {subject_id}.")
    return accession

def standardize_sequence_name(value: str) -> str:
    """Standardizes free text for safe use in a FASTA record identifier."""
    value = value.strip()
    value = re.sub(r"\s+", "_", value)
    value = re.sub(r"[^a-zA-Z0-9_.-]", "_", value)
    value = re.sub(r"_+", "_", value)
    return value.strip("_") or "unknown"

def build_assembled_sequence_name(sample_id: str) -> str:
    """Builds the standardized identifier for a sample's own assembled genome record."""
    return f"{standardize_sequence_name(sample_id)}|assembled"

def build_reference_sequence_name(accession: str, scientific_name: str) -> str:
    """Builds a readable, unique identifier for one reference sequence's tree tip."""
    return f"{standardize_sequence_name(scientific_name)}__{standardize_sequence_name(accession)}"

def _scientific_name_is_missing(value: object) -> bool:
    """Returns whether a scientific name is absent or a known placeholder."""
    return str(value or "").strip().lower() in {"", "n/a", "na", "none", "null", "unknown", "not available"}

def _infer_scientific_name_from_record(record: SeqRecord, accession: str) -> str:
    """Infers a scientific name from an NCBI record's own annotations/description when BLAST's own sscinames is missing."""
    organism = str(record.annotations.get("organism", "")).strip()
    if not _scientific_name_is_missing(organism):
        return organism

    description = record.description.strip()
    if description.startswith(accession):
        description = description[len(accession):].strip()

    match = re.match(r"^([A-Z][A-Za-z.-]+)\s+([a-z][A-Za-z.-]+)", description)
    return f"{match.group(1)} {match.group(2)}" if match is not None else "unknown_species"

def read_assembled_genome(assembly_fasta: str | Path, sample_id: str, logger: Logger | None = None) -> SeqRecord:
    """
    Reads and standardizes one sample's own assembled mitochondrial
    genome for inclusion in its phylogenetic dataset.

    Args:
        assembly_fasta (str | Path): The sample's own selected/resolved assembly FASTA.
        sample_id (str): The sample id.
        logger (Logger | None): Logger to use. Defaults to None.

    Returns:
        SeqRecord: The standardized assembled-genome record.

    Raises:
        FileNotFoundError: If assembly_fasta does not exist.
        ValueError: If the FASTA is empty, carries more than one record, or the sequence itself is empty.
    """
    assembly_fasta = Path(assembly_fasta)
    if not assembly_fasta.is_file():
        raise FileNotFoundError(f"Assembly FASTA not found: {assembly_fasta}.")

    records = list(SeqIO.parse(assembly_fasta, "fasta"))
    if not records:
        raise ValueError(f"No sequences found in assembly FASTA: {assembly_fasta}.")
    if len(records) != 1:
        raise ValueError(f"Expected exactly one assembled sequence in {assembly_fasta}, but found {len(records)}.")

    record = records[0]
    if len(record.seq) == 0:
        raise ValueError(f"Assembly sequence is empty: {assembly_fasta}.")

    record.id = build_assembled_sequence_name(sample_id)
    record.name = record.id
    record.description = ""

    if logger is not None: logger.info(f"Read assembled genome containing {len(record.seq)} bases from {assembly_fasta}.")
    return record

def retrieve_reference_sequences(accessions: list[str], entrez_email: str, entrez_api_key: str | None = None,
                                  logger: Logger | None = None) -> list[SeqRecord]:
    """
    Retrieves real nucleotide reference sequences from NCBI Entrez
    for a sample's own selected BLAST hit accessions -- a live
    network call, never fabricated or cached locally.

    Args:
        accessions (list[str]): The bare NCBI accessions to fetch.
        entrez_email (str): The contact email NCBI requires for Entrez requests.
        entrez_api_key (str | None): An optional NCBI API key.
        logger (Logger | None): Logger to use. Defaults to None.

    Returns:
        list[SeqRecord]: The retrieved records, in NCBI's own returned order.

    Raises:
        ValueError: If accessions is empty or entrez_email is blank.
    """
    if not accessions:
        raise ValueError("At least one accession is required.")
    if not entrez_email.strip():
        raise ValueError("An Entrez email address is required.")

    Entrez.email = entrez_email.strip()
    Entrez.api_key = entrez_api_key.strip() if entrez_api_key else None

    if logger is not None: logger.info(f"Retrieving {len(accessions)} reference sequence(s) from NCBI Entrez.")
    with Entrez.efetch(db = "nucleotide", id = ",".join(accessions), rettype = "fasta", retmode = "text") as handle:
        records = list(SeqIO.parse(handle, "fasta"))
    if logger is not None: logger.info(f"Retrieved {len(records)} reference sequence(s) from NCBI Entrez.")
    return records

def order_and_standardize_reference_sequences(hits: list[BlastHit], retrieved_records: list[SeqRecord],
                                               logger: Logger | None = None) -> list[SeqRecord]:
    """
    Orders Entrez's returned records by the sample's own BLAST rank
    and assigns each a readable, unique tree-tip identifier.

    Args:
        hits (list[BlastHit]): The sample's own selected BLAST hits, in rank order.
        retrieved_records (list[SeqRecord]): The records Entrez returned for those hits' accessions.
        logger (Logger | None): Logger to use. Defaults to None.

    Returns:
        list[SeqRecord]: The standardized reference records, in the same order as hits.

    Raises:
        ValueError: If a hit's own accession was not among the retrieved records, or a retrieved sequence is empty.
    """
    records_by_accession: dict[str, SeqRecord] = {}
    for record in retrieved_records:
        records_by_accession[record.id] = record
        records_by_accession[record.id.split(".")[0]] = record

    ordered_records: list[SeqRecord] = []
    for hit in hits:
        accession = normalize_blast_accession(hit.sseqid)
        record = records_by_accession.get(accession) or records_by_accession.get(accession.split(".")[0])
        if record is None:
            raise ValueError(f"NCBI did not return the expected reference sequence: {accession}.")

        record = record[:]
        scientific_name = str(hit.sscinames or "").strip()
        if _scientific_name_is_missing(scientific_name):
            scientific_name = _infer_scientific_name_from_record(record, accession)

        record.id = build_reference_sequence_name(accession, scientific_name)
        record.name = record.id
        record.description = ""

        if len(record.seq) == 0:
            raise ValueError(f"Retrieved reference sequence is empty: {accession}.")
        ordered_records.append(record)

    if logger is not None: logger.info(f"Ordered and standardized {len(ordered_records)} reference sequence(s).")
    return ordered_records

def validate_sample_dataset(records: list[SeqRecord], expected_reference_count: int) -> None:
    """
    Validates a sample's combined phylogenetic dataset before it's
    written to disk.

    Args:
        records (list[SeqRecord]): The combined records (assembled genome first, references after).
        expected_reference_count (int): How many reference records should be present.

    Raises:
        ValueError: If the total count is wrong, the assembled genome isn't first, identifiers collide, or any sequence is empty.
    """
    expected_total = expected_reference_count + 1
    if len(records) != expected_total:
        raise ValueError(f"Expected {expected_total} total sequence(s) but found {len(records)}.")
    if not records[0].id.endswith("|assembled"):
        raise ValueError("The assembled genome must be the first record.")

    identifiers = [record.id for record in records]
    if len(identifiers) != len(set(identifiers)):
        raise ValueError("Duplicate sequence identifiers were generated.")
    for record in records:
        if len(record.seq) == 0:
            raise ValueError(f"Sequence is empty: {record.id}.")

def build_sample_phylogeny_dataset(sample_id: str, assembly_fasta: str | Path, top_hits_tsv: str | Path,
                                    output_fasta: str | Path, entrez_email: str, entrez_api_key: str | None = None,
                                    logger: Logger | None = None,
                                    retrieve_function: Callable[..., list[SeqRecord]] | None = None) -> int:
    """
    Builds and writes one sample's own real phylogenetic dataset: its
    assembled genome plus the real sequences for its own top-ranked
    BLAST hits, fetched live from NCBI Entrez.

    Args:
        sample_id (str): The sample id.
        assembly_fasta (str | Path): The sample's own selected/resolved assembly FASTA.
        top_hits_tsv (str | Path): The sample's own ranked BLAST top-hit TSV (write_top_blast_hits()'s own format).
        output_fasta (str | Path): Where to write the combined FASTA.
        entrez_email (str): The contact email NCBI requires for Entrez requests.
        entrez_api_key (str | None): An optional NCBI API key.
        logger (Logger | None): Logger to use. Defaults to None.
        retrieve_function (Callable[..., list[SeqRecord]] | None): Overrides retrieve_reference_sequences -- real callers never pass this; tests inject a fake to avoid a live network call.

    Returns:
        int: The number of reference sequences included.

    Raises:
        FileNotFoundError: If assembly_fasta or top_hits_tsv does not exist.
        ValueError: If the dataset can't be built into a valid combined FASTA.
    """
    output_fasta = Path(output_fasta)
    if retrieve_function is None:
        retrieve_function = retrieve_reference_sequences

    assembly_record = read_assembled_genome(assembly_fasta, sample_id, logger = logger)
    hits = read_top_blast_hits(top_hits_tsv, sample_id = sample_id, logger = logger)
    if not hits:
        raise ValueError(f"No selected BLAST hits found in {top_hits_tsv}.")

    accessions = [normalize_blast_accession(hit.sseqid) for hit in hits]
    retrieved_records = retrieve_function(accessions = accessions, entrez_email = entrez_email,
                                          entrez_api_key = entrez_api_key, logger = logger)
    reference_records = order_and_standardize_reference_sequences(hits, retrieved_records, logger = logger)

    combined_records = [assembly_record, *reference_records]
    validate_sample_dataset(combined_records, expected_reference_count = len(hits))

    output_fasta.parent.mkdir(parents = True, exist_ok = True)
    SeqIO.write(combined_records, output_fasta, "fasta")

    if logger is not None: logger.info(f"Wrote {len(combined_records)} sequence(s) to {output_fasta}.")
    return len(reference_records)
