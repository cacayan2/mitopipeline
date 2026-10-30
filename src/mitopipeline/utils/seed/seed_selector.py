"""
seed_selector.py

Orchestrates mitochondrial seed selection for one sample: tries each
locus in a documented-density-ordered chain against MitoFish, fetches
the first clean match's sequence, and writes it as a seed FASTA with
full provenance. Never raises or exits over a single sample's lookup
failure -- any failure (unknown organism, no clean match at any locus,
corrupt/unreachable database) is reported as a typed fallback result
instead, so the caller can proceed with GetOrganelle's own default
database for the organelle type (see api/getorganelle.py's "no seed"
run) rather than the whole pipeline crashing over one sample.
"""

# Imports
import hashlib
import logging
from pathlib import Path
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from mitopipeline.utils.seed.gene_extraction import write_gene_fasta
from mitopipeline.utils.seed.mitofish import MitoFishReferenceLookup
from mitopipeline.utils.seed.models import DEFAULT_SEED_LOCUS, LOCUS_FALLBACK_ORDER, SeedSelectionResult

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
        gene (str): The mitochondrial locus that was requested (or the first locus tried, for a whole-chain exhaustion).
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
                 loci: tuple[str, ...] = LOCUS_FALLBACK_ORDER) -> SeedSelectionResult:
    """
    Selects and writes a mitochondrial seed sequence for one sample,
    trying each locus in loci in order until one yields a MitoFish
    candidate with a fetchable sequence. Multiple candidates at a
    locus are resolved deterministically (longest available sequence
    record first, accession as a final tie-breaker) rather than
    skipped as "ambiguous" -- a well-documented organism (e.g. Danio
    rerio, with 16 complete MitoFish genomes) still deserves a real,
    reproducible seed, not an unseeded fallback just because more
    than one real reference happens to exist. Falls back to no seed
    (GetOrganelle's own default database) only once every locus in
    the chain has been tried and none had a fetchable candidate.

    Args:
        sample_id (str): The sample to select a seed for.
        scientific_name (str): The sample's scientific name.
        database_dir (Path): The path to the local MitoFish reference database.
        output_dir (Path): The directory to write the selected seed FASTA to.
        logger (logging.Logger): The logger to use for logging messages.
        loci (tuple[str, ...]): The mitochondrial loci to try, in order, most densely documented first. Defaults to LOCUS_FALLBACK_ORDER.

    Returns:
        SeedSelectionResult: Always returned, never raises -- fallback_used
        is True and fasta_path is None when no locus could be resolved.
    """
    context = f"{{{sample_id}}}"
    lookup = MitoFishReferenceLookup(database_dir = database_dir, logger = logger, sample_id = sample_id)

    for gene in loci:
        # Looking up MitoFish candidates for this locus. A SystemExit
        # here means the local database itself is missing/corrupt --
        # the identical failure every other locus would hit too, so
        # this is a hard stop, not a reason to keep trying.
        try:
            candidates = lookup.lookup(scientific_name = scientific_name, gene = gene)
        except SystemExit:
            logger.warning(f"{context} MitoFish lookup failed for {scientific_name}; falling back to GetOrganelle's default database.")
            return _fallback(sample_id, scientific_name, gene, "mitofish_lookup_failed")

        # No candidate at this locus -- trying the next one in the chain.
        if not candidates:
            logger.info(f"{context} No MitoFish candidate for {scientific_name} ({gene}); trying the next locus.")
            continue

        multiple_candidates = len(candidates) > 1
        if multiple_candidates:
            logger.info(f"{context} {len(candidates)} MitoFish candidates for {scientific_name} ({gene}); selecting the best-documented one deterministically.")
            candidates = sorted(candidates, key = lambda c: (-(c.sequence_length or 0), c.accession))

        # Trying each candidate at this locus, best (longest, most
        # complete record) first -- a missing accession (metadata/
        # sequence-file sync issue) is specific to one candidate, not
        # the whole locus, so the next-best candidate at the SAME
        # locus is tried before giving up on this locus entirely.
        for candidate in candidates:
            try:
                sequence = lookup.fetch_sequence(candidate.accession)
            except SystemExit:
                logger.info(f"{context} Failed to fetch sequence for accession {candidate.accession} ({gene}); trying the next candidate.")
                continue

            # Writing the seed FASTA and returning the successful result.
            seq_record = SeqRecord(Seq(sequence), id = candidate.accession,
                                   description = f"{gene} seed for {sample_id}, {scientific_name}")
            fasta_path = output_dir / f"{sample_id}.seed.fasta"
            write_gene_fasta(seq_record, fasta_path)

            logger.info(f"{context} Selected MitoFish seed {candidate.accession} ({gene}) for {scientific_name}; written to {fasta_path}.")

            return SeedSelectionResult(
                sample_id = sample_id,
                scientific_name = scientific_name,
                gene = gene,
                accession = candidate.accession,
                source = "mitofish",
                fasta_path = fasta_path,
                fallback_used = False,
                match_level = "best_of_multiple" if multiple_candidates else "exact",
                checksum = _checksum(sequence),
                ncbi_taxon_id = candidate.ncbi_taxon_id,
                sequence_length = len(sequence),
            )

        logger.info(f"{context} No fetchable MitoFish candidate for {scientific_name} ({gene}); trying the next locus.")

    # Every locus in the chain was exhausted without a fetchable match.
    logger.warning(f"{context} No usable MitoFish match for {scientific_name} across {len(loci)} loci; falling back to GetOrganelle's default database.")
    return _fallback(sample_id, scientific_name, loci[0] if loci else DEFAULT_SEED_LOCUS, "no_fetchable_mitofish_match")
