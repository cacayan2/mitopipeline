"""
gene_extraction.py

Extracts a named mitochondrial gene region from a reference sequence
record, handling common gene-name synonyms and strand orientation, and
writes the extracted region to a seed FASTA file.
"""

# Imports
from __future__ import annotations
import logging
import sys
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# Known alternate names for the same mitochondrial locus, grouped under
# one canonical key. Matching is case-insensitive and ignores surrounding
# whitespace -- see _canonical_gene().
GENE_SYNONYMS: dict[str, frozenset[str]] = {
    "12S_rRNA": frozenset({"12s rrna", "12s ribosomal rna", "s-rrna", "rrns"}),
    "16S_rRNA": frozenset({"16s rrna", "16s ribosomal rna", "l-rrna", "rrnl"}),
}

# Feature types that may carry a gene's name/product qualifiers in a
# GenBank-style record.
_GENE_FEATURE_TYPES = ("gene", "rRNA", "CDS")

def _canonical_gene(name: str) -> str | None:
    """
    Resolves a gene name or synonym to its canonical locus name.

    Args:
        name (str): The gene name or synonym to resolve.

    Returns:
        str | None: The canonical locus name, or None if unrecognized.
    """
    normalized = name.strip().lower()
    for canonical, synonyms in GENE_SYNONYMS.items():
        if normalized == canonical.lower() or normalized in synonyms:
            return canonical
    return None

def extract_gene(record: SeqRecord, gene: str, logger: logging.Logger) -> SeqRecord:
    """
    Extracts the named gene region from a reference record.

    Strand orientation is handled by Bio.SeqFeature's own location.extract(),
    which reverse-complements the extracted sequence automatically when the
    feature is annotated on the minus strand -- this function does not
    reimplement that logic.

    Args:
        record (SeqRecord): The reference sequence record to extract from.
        gene (str): The requested gene name or synonym.
        logger (logging.Logger): The logger to use for logging messages.

    Returns:
        SeqRecord: The extracted gene sequence, correctly oriented.
    """
    # Resolving the requested gene name to its canonical locus.
    canonical = _canonical_gene(gene)
    if canonical is None:
        logger.error(f"Unrecognized gene name or synonym: {gene}.")
        sys.exit(1)

    # Searching the record's features for a matching gene/product qualifier.
    for feature in record.features:
        if feature.type not in _GENE_FEATURE_TYPES:
            continue

        qualifier_values = [value for values in feature.qualifiers.values() for value in values]
        if not any(_canonical_gene(str(value)) == canonical for value in qualifier_values):
            continue

        # Extracting the gene region; location.extract() reverse-complements
        # automatically for a minus-strand feature.
        extracted_seq = feature.location.extract(record.seq)
        logger.info(f"Extracted {canonical} from {record.id} (strand {feature.location.strand}).")
        return SeqRecord(
            extracted_seq,
            id = f"{record.id}_{canonical}",
            description = f"{canonical} extracted from {record.id}",
        )

    # Handling a reference record missing the requested locus entirely.
    logger.error(f"Gene {gene} (canonical: {canonical}) not found in reference {record.id}.")
    sys.exit(1)

def write_gene_fasta(seq_record: SeqRecord, output_path: Path) -> Path:
    """
    Writes an extracted gene SeqRecord to a FASTA file.

    Args:
        seq_record (SeqRecord): The gene sequence to write.
        output_path (Path): The path to write the FASTA file to.

    Returns:
        Path: The path the FASTA file was written to.
    """
    # Creating the parent directory if it does not exist.
    output_path.parent.mkdir(parents = True, exist_ok = True)

    # Writing the sequence record to disk as a FASTA file.
    SeqIO.write(seq_record, output_path, "fasta")

    return output_path
