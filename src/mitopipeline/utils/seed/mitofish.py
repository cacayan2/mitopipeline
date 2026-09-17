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
        accession_taxon_ids = self._lookup_accessions(taxon_ids)

        # Handling organisms without mitochondrial sequence records.
        if not accession_taxon_ids:
            self.logger.warning(f"{context} No mitochondrial accessions found for {scientific_name}.")
            return []

        # Building candidates matching the requested gene.
        candidates = self._build_candidates(accession_taxon_ids, scientific_name, gene)

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

    def _lookup_taxon_ids(self, scientific_name: str) -> list[str]:
        """
        Matches an exact scientific name to NCBI taxonomy ID's using MitoFish.

        Args:
            scientific_name (str): The scientific name of the organism.

        Returns:
            list[str]: A list of matching NCBI taxonomy ID's. 
        """
        # Obtaining the logging context.
        context = self._log_context()

        # Defining paths to the required MitoFish taxonomy tables.
        taxonomy_name_path = self.database_dir / "taxonid_name.parquet"
        species_lineage_path = self.database_dir / "speciesid_lineageid.parquet"
        taxon_species_path = self.database_dir / "taxonid_speciesid.parquet"

        # Validating the required taxonomy tables.
        required_tables = [taxonomy_name_path, species_lineage_path, taxon_species_path]
        for path in required_tables:
            if not path.is_file():
                self.logger.error(f"{context} Required MitoFish taxonomy table not found at {path}.")
                sys.exit(1)

        # Loading the required MitoFish taxonomy tables.
        try:
            taxonomy_names = pd.read_parquet(taxonomy_name_path)
            species_lineages = pd.read_parquet(species_lineage_path)
            taxon_species = pd.read_parquet(taxon_species_path)
        except Exception as e:
            self.logger.error(f"{context} Failed to load MitoFish taxonomy tables: {e}")
            sys.exit(1)

        # Resolving the exact scientific name to MitoFish lineage IDs.
        name_matches = taxonomy_names[taxonomy_names["lineage_name"] == scientific_name]

        # Handling scientific anmes absent from MitoFish.
        if name_matches.empty:
            self.logger.warning(f"{context} Scientific name {scientific_name} not found in MitoFish.")
            return []

        # Extracting unique matching lineage IDs.
        lineage_ids = name_matches["lineage_id"].dropna().drop_duplicates().tolist()

        # Resolving lineage IDs to MitoFish species IDs. 
        lineage_matches = species_lineages[species_lineages["lineage_id"].isin(lineage_ids)]

        # Handling lineage IDs without associated species IDs.
        if lineage_matches.empty:
            self.logger.warning(f"{context} Lineage IDs {lineage_ids} not found in MitoFish.")
            return []

        # Extracting unique matching species IDs.
        species_ids = (lineage_matches["species_id"].dropna().drop_duplicates().tolist())

        # Resolving MitoFish species IDs to NCBi taxonomy IDs.
        taxon_matches = taxon_species[taxon_species["species_id"].isin(species_ids)]

        # Handling species IDs without associated taxonomy IDs.
        if taxon_matches.empty:
            self.logger.warning(f"{context} Species IDs {species_ids} not found in MitoFish.")
            return []

        # Extracting unique matching taxonomy IDs.
        taxon_ids = (taxon_matches["taxon_id"].dropna().drop_duplicates().tolist())

        # Logging successful lookup completion. 
        self.logger.info(f"{context} Found {len(taxon_ids)} MitoFish taxonomy ID(s) for {scientific_name}.")
        self.logger.debug(f"{context} Taxon IDs for {scientific_name}: {taxon_ids}")
        return taxon_ids

    def _lookup_accessions(self, taxon_ids: list[str]) -> dict[str, str]:
        """
        Matches NCBI taxonomy ID's to mitochondrial accessions using MitoFish.

        Args:
            taxon_ids (list[str]): A list of NCBI taxonomy ID's.

        Returns:
            dict[str, str]: A dictionary of mitochondrial accessions keyed by NCBI taxonomy ID.    
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
            return {}

        # Removing records without an accession or taxonomy ID.
        matches = matches.dropna(subset = ["accession", "taxon_id"])

        # Detecting accessions associated with multiple taxonomy IDs.
        taxon_counts = matches[["accession", "taxon_id"]].drop_duplicates().groupby("accession")["taxon_id"].nunique()
        ambiguous_accessions = taxon_counts[taxon_counts > 1].index.tolist()
        if ambiguous_accessions:
            self.logger.error(f"{context} Multiple taxonomy IDs found for mitochondrial accessions: {ambiguous_accessions}.")
            sys.exit(1)

        # Building the accession-to-taxonomy-ID mapping.
        accession_taxon_ids = matches[["accession", "taxon_id"]].drop_duplicates(subset = ["accession"]).set_index("accession")["taxon_id"].to_dict()

        # Logigng successful accession lookup.
        self.logger.info(f"{context} Found {len(accession_taxon_ids)} mitochondrial accession(s) for taxonomy ID(s): {taxon_ids}.")
        self.logger.debug(f"{context} Mitochondrial accession(s) for taxonomy ID(s) {taxon_ids}: {accession_taxon_ids}")

        return accession_taxon_ids

    def _build_candidates(self, accession_taxon_ids: dict[str, str], scientific_name: str, gene: str) -> list[SeedReferenceCandidate]:
        """
        Builds seed reference candidates from MitoFish annotations.

        Args:
            accession_taxon_ids (dict[str, str]): A mapping of mitochondrial accessions to NCBI taxonomy IDs.
            scientific_name (str): The scientific name of the organism.
            gene (str): The mitochondrial gene to query.

        Returns:
            list[SeedReferenceCandidate]: A list of candidate reference sequences.
        """
        # Obtaining the logging context.
        context = self._log_context()

        # Defining the path to the MitoFish sequence annotation table.
        annotation_path = self.database_dir / "seq_annotation.parquet"

        # Validating the sequence annotation table path.
        if not annotation_path.is_file():
            self.logger.error(f"{context} MitoFish sequence annotation table not found at {annotation_path}.")
            sys.exit(1)

        # Loading the MitoFish sequence annotation table.
        try:
            annotations = pd.read_parquet(annotation_path)
        except Exception as e:
            self.logger.error(f"{context} Failed to load MitoFish sequence annotation table from {annotation_path}: {e}")
            sys.exit(1)

        # Restricting annotations to accessions associated with the organism.
        matches = annotations[annotations["accession"].isin(accession_taxon_ids)]

        # Handling accessions without annotation records.
        if matches.empty:
            self.logger.warning(f"{context} No MitoFish annotation records found for {scientific_name}.")
            return []

        # Restricting annotation records to the requested mitochondrial gene.
        gene_matches = matches[matches["gene"] == gene]

        # Handling accessions without the requested gene.
        if gene_matches.empty:
            self.logger.warning(f"{context} No MitoFish annotation records found for {gene} in {scientific_name}.")
            return []

        # Building one candidate for each matching reference record.
        candidates: list[SeedReferenceCandidate] = []

        # Iterating through each matching reference record and populating the candidate list.
        for _, record in gene_matches.iterrows():
            accession = str(record["accession"])

            # Creating the candidate from the MitoFish annotation record.
            candidate = SeedReferenceCandidate(
                accession = accession,
                scientific_name = scientific_name,
                ncbi_taxon_id = int(accession_taxon_ids[accession]),
                gene = gene,
                sequence_length = int(record["length"]) if pd.notna(record["length"]) else None
            )

            candidates.append(candidate)

        # Removing duplicate candidates and preserving order.
        candidates = list(dict.fromkeys(candidates))

        # Logging succcessful candidate construction.
        self.logger.info(f"{context} Found {len(candidates)} candidate reference(s) for {scientific_name}, gene {gene}.")
        return candidates