"""
seed_selector.py

Orchestrates mitochondrial seed selection for one sample: looks up a
MitoFish reference, fetches its sequence, and writes it as a seed FASTA
with full provenance. Never raises or exits over a single sample's
lookup failure -- any failure (unknown organism, ambiguous match,
corrupt/unreachable database) is reported as a typed fallback result
instead, so the caller can proceed with GetOrganelle's own default
database for the organelle type (see api/getorganelle.py's "no seed"
run, #393) rather than the whole pipeline crashing over one sample.
"""

# Imports
import hashlib
import logging
from pathlib import Path
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from mitopipeline.utils.seed.gene_extraction import write_gene_fasta
from mitopipeline.utils.seed.mitofish import MitoFishReferenceLookup
from mitopipeline.utils.seed.models import DEFAULT_SEED_LOCUS, SeedSelectionResult

# The accession recorded on a fallback result, where no real seed
# accession exists.
FALLBACK_ACCESSION = "GETORGANELLE_DEFAULT"
FALLBACK_SOURCE = "fallback_default_database"

def _checksum(sequence: str) -> str:
    """
    Computes a SHA-256 checksum of a sequence, for provenance.

    Args:
        sequence (str): The sequence to checksum.

    Returns:
        str: The hex-encoded SHA-256 checksum.
    """
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()

def _fallback(sample_id: str, scientific_name: str, gene: str, reason: str) -> SeedSelectionResult:
    """
    Builds the typed "no seed" result for a sample, recording why.

    Args:
        sample_id (str): The sample this seed selection was performed for.
        scientific_name (str): The scientific name the selection was queried against.
        gene (str): The mitochondrial locus that was requested.
        reason (str): Why a fallback was needed.

    Returns:
        SeedSelectionResult: A result with fallback_used True and no FASTA.
    """
    return SeedSelectionResult(
        sample_id = sample_id,
        scientific_name = scientific_name,
        gene = gene,
        accession = FALLBACK_ACCESSION,
        source = FALLBACK_SOURCE,
        fallback_used = True,
        fallback_reason = reason,
    )

def select_seed(sample_id: str,
                 scientific_name: str,
                 database_dir: Path,
                 output_dir: Path,
                 logger: logging.Logger,
                 gene: str = DEFAULT_SEED_LOCUS) -> SeedSelectionResult:
    """
    Selects and writes a mitochondrial seed sequence for one sample.

    Args:
        sample_id (str): The sample to select a seed for.
        scientific_name (str): The sample's scientific name.
        database_dir (Path): The path to the local MitoFish reference database.
        output_dir (Path): The directory to write the selected seed FASTA to.
        logger (logging.Logger): The logger to use for logging messages.
        gene (str): The mitochondrial locus to select. Defaults to DEFAULT_SEED_LOCUS.

    Returns:
        SeedSelectionResult: Always returned, never raises -- fallback_used
        is True and fasta_path is None when no seed could be selected.
    """
    context = f"{{{sample_id}}}"
    lookup = MitoFishReferenceLookup(database_dir = database_dir, logger = logger, sample_id = sample_id)

    # Looking up MitoFish candidates. Any failure here (a corrupt or
    # unreachable local database, missing metadata tables -- the same
    # failure a transient network error downloading the database would
    # eventually surface as) is a fallback, not a crash.
    try:
        candidates = lookup.lookup(scientific_name = scientific_name, gene = gene)
    except SystemExit:
        logger.warning(f"{context} MitoFish lookup failed for {scientific_name}; falling back to GetOrganelle's default database.")
        return _fallback(sample_id, scientific_name, gene, "mitofish_lookup_failed")

    # Handling an organism absent from MitoFish, or the locus missing for it.
    if not candidates:
        logger.warning(f"{context} No MitoFish candidate for {scientific_name} ({gene}); falling back.")
        return _fallback(sample_id, scientific_name, gene, "no_mitofish_candidate")

    # Rejecting an ambiguous match rather than silently guessing one.
    if len(candidates) > 1:
        logger.warning(f"{context} {len(candidates)} ambiguous MitoFish candidates for {scientific_name} ({gene}); falling back rather than guessing.")
        return _fallback(sample_id, scientific_name, gene, "ambiguous_candidates")

    candidate = candidates[0]

    # Fetching the candidate's actual sequence.
    try:
        sequence = lookup.fetch_sequence(candidate.accession)
    except SystemExit:
        logger.warning(f"{context} Failed to fetch sequence for accession {candidate.accession}; falling back.")
        return _fallback(sample_id, scientific_name, gene, "sequence_fetch_failed")

    # Writing the seed FASTA.
    seq_record = SeqRecord(Seq(sequence), id = candidate.accession,
                           description = f"{gene} seed for {sample_id}, {scientific_name}")
    fasta_path = output_dir / f"{sample_id}.seed.fasta"
    write_gene_fasta(seq_record, fasta_path)

    logger.info(f"{context} Selected MitoFish seed {candidate.accession} for {scientific_name} ({gene}); written to {fasta_path}.")

    return SeedSelectionResult(
        sample_id = sample_id,
        scientific_name = scientific_name,
        gene = gene,
        accession = candidate.accession,
        source = "mitofish",
        fasta_path = fasta_path,
        fallback_used = False,
        match_level = "exact",
        checksum = _checksum(sequence),
        ncbi_taxon_id = candidate.ncbi_taxon_id,
        sequence_length = len(sequence),
    )
