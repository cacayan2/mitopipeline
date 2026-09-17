"""
mitofish.py

Lookup of mitochondrial refernece candidates from MitoFish.
"""

# Imports
from pathlib import Path
import logging
import sys
import pandas as pd
from mitopipeline.utils.seed.models import SeedReferenceCandidate

class MitoFishReferenceLookup:
    """
    Queries a local MitoFish reference database. 
    """
    def __init__(self, database_dir: Path, logger: logging.Logger, sample_id: str | None = None) -> None:
        """
        Initializes the MitoFishReferenceLookup with the path to the local database.

        Args:
            database_dir (Path): The path to the local MitoFish reference database.
            logger (logging.Logger): The logger to use for logging messages.

        Returns:
            None
        """
        self.database_dir = Path(database_dir)
        self.logger = logger
        self.sample_id = sample_id

    def lookup(self, scientific_name: str, gene: str = "12S_rRNA") -> list[SeedReferenceCandidate]:
        """
        Returns the candidate references for an organism and gene.

        Args:
            scientific_name (str): The scientific name of the organism.
            gene (str): The mitochondrial gene to query. Defaults to "12S_rRNA".

        Returns:
            list[SeedReferenceCandidate]: A list of candidate reference sequences. 
        """
        # Obtaining the context (sample) for the logger.
        context = self._log_context()

        # Normalizing the scientific name and gene name. 
        scientific_name = scientific_name.strip()
        gene = gene.strip()

        # Validating that the scientific name and gene are not empty.
        if not scientific_name: 
            self.logger.error(f"{context} Scientific name is empty.")
            sys.exit(1)
        if not gene: 
            self.logger.error(f"{context} Gene name is empty.")
            sys.exit(1)

        # Logging the start of the lookup.
        if self.logger is not None:
            self.logger.info(f"{context} Initiating MitoFish reference lookup for {scientific_name}, gene {gene}.")

        # Looking up the organism's NCBI taxonomy IDs.
        taxon_ids = self._lookup_taxon_ids(scientific_name)

        # Handling organisms not represented in MitoFish.
        if not taxon_ids:
            self.logger.warning(f"{context} No MitoFish taxonomy entry found for {scientific_name}.")
            return []

        # Looking up mitochondrial accessions associated with the organism.
        accessions = self._lookup_accessions(taxon_ids)

        # Handling organisms without mitochondrial sequence records.
        if not accessions:
            self.logger.warning(f"{context} No mitochondrial accessions found for {scientific_name}.")
            return []

        # Building candidates matching the requested gene.
        candidates = self._build_candidates(accessions, scientific_name, gene)

        # Logging an unsuccessful gene lookup.
        if not candidates:
            self.logger.warning(f"{context} No MitoFish reference candidatees found for {scientific_name}, gene {gene}.")
            return []

        # Logging successful lookup completion. 
        self.logger.info(f"{context} Found {len(candidates)} MitoFish reference candidate(s) for {scientific_name}, gene {gene}.")
        return candidates

    def _log_context(self) -> str:
        """
        Adds a sample label for log messages.
        
        Returns:
            str: The sample label for log messages.
        """
        # Attempting to access a direct sample ID.
        sample_id = getattr(self, "sample_id", None)

        # Returning the formatted sample label.
        if sample_id: 
            return f"{{{sample_id}}}"

        return ""

    def _lookup_taxon_ids(self, scientific_name: str) -> list[int]:
        """
        Matches an exact scientific name to NCBI taxonomy ID's using MitoFish.

        Args:
            scientific_name (str): The scientific name of the organism.

        Returns:
            list[int]: A list of NCBI taxonomy ID's. 
        """
        # Obtaining the logging context.
        context = self._log_context()

        # Defining the path to the MitoFish tasxonomy name table.
        taxonomy_path = self.database_dir / "taxonid_name.parquet"

        # Validating the taxonomy name table path.
        if not taxonomy_path.is_file():
            self.logger.error(f"{context} MitoFish taxonomy name table not found at {taxonomy_path}.")
            sys.exit(1)

        # Loading the MitoFish taxonomy name table.
        try:
            taxonomy = pd.read_parquet(taxonomy_path)
        except Exception as e:
            self.logger.error(f"{context} Failed to load MitoFish taxonomy name table from {taxonomy_path}: {e}")
            sys.exit(1)

        # Resolving the exact normalized scientific name.
        matches = taxonomy[taxonomy["name"] == scientific_name]

        # Handling names absent from the MitoFish taxonomy table.
        if matches.empty:
            self.logger.warning(f"{context} No MitoFish taxonomy entry found for {scientific_name}.")
            return []

        # Extracting unique matching NCBI taxonomy IDs.
        taxon_ids = (matches["taxon_id"].dropna().astype(int).drop_duplicates().tolist())

        # Logging successful taxonomy ID lookup.
        self.logger.info(f"{context} Found {len(taxon_ids)} MitoFish taxonomy ID(s) for {scientific_name}.")
        self.logger.debug(f"{context} MitoFish taxonomy ID(s) for {scientific_name}: {taxon_ids}")

        return taxon_ids

    def _lookup_accessions(self, taxon_ids: list[int]) -> list[str]:
        """
        Matches NCBI taxonomy ID's to mitochondrial accessions using MitoFish.

        Args:
            taxon_ids (list[int]): A list of NCBI taxonomy ID's.

        Returns:
            list[str]: A list of mitochondrial accessions.         
        """
        # Obtaining the logging context.
        context = self._log_context()

        # Defining the path to the MitoFish sequence taxonomy table.
        sequence_taxonomy_path = self.database_dir / "seq_taxonid.parquet"

        # Validating the sequence taxonomy table path.
        if not sequence_taxonomy_path.is_file():
            self.logger.error(f"{context} MitoFish sequence taxonomy table not found at {sequence_taxonomy_path}.")
            sys.exit(1)

        # Loading the MitoFish sequence taxonomy table.
        try:
            sequence_taxonomy = pd.read_parquet(sequence_taxonomy_path)
        except Exception as e:
            self.logger.error(f"{context} Failed to load MitoFish sequence taxonomy table from {sequence_taxonomy_path}: {e}")
            sys.exit(1)

        # Resolving taxonomy IDs against MitoFish sequence metadata.
        matches = sequence_taxonomy[sequence_taxonomy["taxon_id"].isin(taxon_ids)]

        # Handling taxonomy IDs without sequence records.
        if matches.empty:
            self.logger.warning(f"{context} No mitochondrial accessions found for taxonomy ID(s): {taxon_ids}.")
            return []

        # Extracting unique sequence accessions.
        accessions = (matches["accession"].dropna().drop_duplicates().tolist())

        # Logigng successful accession lookup.
        self.logger.info(f"{context} Found {len(accessions)} mitochondrial accession(s) for taxonomy ID(s): {taxon_ids}.")
        self.logger.debug(f"{context} Mitochondrial accession(s) for taxonomy ID(s) {taxon_ids}: {accessions}")

        return accessions

    def _build_candidates(self, accessions: list[str], scientific_name: str, gene: str) -> list[SeedReferenceCandidate]:
        """
        Builds seed reference candidates from MitoFish annotations.

        Args:
            accessions (list[str]): MitoFish mitochondrial accessions associated with the organism.
            scientific_name (str): The scientific name of the organism.
            gene (str): The mitochondrial gene to query.
        
        Returns:
            list[SeedReferenceCandidate]: A list of seed reference candidates.
        """
        # Obtaining the logging context.
        context = self._log_context()

        # Defining the path to the MitoFish sequence annotation table.
        annotation_path = self.database_dir / "seq_annotation.parquet"

        # Validating the sequence annotation table path.
        if not annotation_path.is_file():
            self.logger.error(f"{context} Mitofish sequence annotation table not found at {annotation_path}.")
            sys.exit(1)

        # Loading the MitoFish sequence annotation table.
        try: annotations = pd.read_parquet(annotation_path)
        except Exception as e:
            self.logger.error(f"{context} Failed to load MitoFish sequence annotation table from {annotation_path}: {e}")
            sys.exit(1)

        # Restricting annotations to accessions associated with the organism of interest.
        matches = annotations[annotations["accession"].isin(accessions)]

        # Handling accessions without any records.
        if matches.empty:
            self.logger.warning(f"{context} No MitoFish annotation records found for {scientific_name}.")
            return []



        # Restricting annotation records to the requested mitochondrial gene.
        gene_matches = matches[matches["gene"] == gene]

        # Handling accessions without the requested gene.
        if gene_matches.empty:
            self.logger.warning(f"{context} No MitoFish annotation records found for {scientific_name}, gene {gene}.")
            return []

        # Building one candidate for each matching reference record.
        candidates: list[SeedReferenceCandidate] = []

        # Iterating through the gene matches and building the candidates.
        for _, record in gene_matches.iterrows():
            # Extracting the sequence length.
            sequence_length = int(record["length"]) if pd.notna(record["length"]) else None

            # Creating the candidate.
            candidate = SeedReferenceCandidate(
                accession = str(record["accession"]),
                scientific_name = scientific_name, 
                ncbi_taxon_id = int(record["taxon_id"]),
                gene = gene,
                sequence_length = sequence_length
            )

            # Adding the candidate to the list.
            candidates.append(candidate)

        # Remocing duplicate candidates and preserving order.
        candidates = list(dict.fromkeys(candidates))

        # Logigng successful candidate construction.
        self.logger.info(f"{context} Built {len(candidates)} MitoFish reference candidate(s) for {scientific_name}, gene {gene}.")

        return candidates