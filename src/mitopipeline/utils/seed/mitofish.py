"""
mitofish.py

Lookup of mitochondrial refernece candidates from MitoFish.
"""

# Imports
import gzip
from pathlib import Path
import logging
import sys
import pandas as pd
from Bio import SeqIO
from mitopipeline.utils.seed.models import DEFAULT_SEED_LOCUS, SeedReferenceCandidate

class MitoFishReferenceLookup:
    """
    Queries a local MitoFish reference database.
    """
    # The bgzipped multi-FASTA of MitoFish's own already gene-extracted
    # sequences, keyed by accession -- see MitoFishDatabase.SEQUENCE_FILE.
    SEQUENCE_FILE = "mitofishdb.fa.gz"

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

    def lookup(self, scientific_name: str, gene: str = DEFAULT_SEED_LOCUS) -> list[SeedReferenceCandidate]:
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
        coordinate_columns = self._coordinate_columns(annotations.columns)

        # Iterating through each matching reference record and populating the candidate list.
        for _, record in gene_matches.iterrows():
            accession = str(record["accession"])
            gene_start, gene_end, gene_strand = self._coordinates_from_row(record, coordinate_columns)

            # Creating the candidate from the MitoFish annotation record.
            candidate = SeedReferenceCandidate(
                accession = accession,
                scientific_name = scientific_name,
                ncbi_taxon_id = int(accession_taxon_ids[accession]),
                gene = gene,
                sequence_length = int(record["length"]) if pd.notna(record["length"]) else None,
                gene_start = gene_start,
                gene_end = gene_end,
                gene_strand = gene_strand,
            )

            candidates.append(candidate)

        # Removing duplicate candidates and preserving order.
        candidates = list(dict.fromkeys(candidates))

        # Logging succcessful candidate construction.
        self.logger.info(f"{context} Found {len(candidates)} candidate reference(s) for {scientific_name}, gene {gene}.")
        return candidates

    # Column names the MitoFish annotation table may use for a gene's
    # coordinates within its record. Matched case-insensitively; the
    # first present name wins. Without a start AND an end column the
    # table carries no coordinates and every candidate is returned
    # without them (seed_selector's mitofish_gene mode then declines to
    # seed from the whole record rather than guess).
    _START_COLUMNS = ("start", "gene_start", "start_pos", "start_position", "from", "begin")
    _END_COLUMNS = ("end", "gene_end", "end_pos", "end_position", "to", "stop")
    _STRAND_COLUMNS = ("strand", "direction", "orientation")

    @classmethod
    def _coordinate_columns(cls, columns) -> tuple[str, str, str | None] | None:
        """
        Finds the annotation table's (start, end, strand) column names.

        Args:
            columns: The annotation table's column index.

        Returns:
            tuple[str, str, str | None] | None: The real (start, end, strand-or-None) column names, or None when the table has no start/end pair.
        """
        by_lower = {str(column).lower(): str(column) for column in columns}
        start = next((by_lower[name] for name in cls._START_COLUMNS if name in by_lower), None)
        end = next((by_lower[name] for name in cls._END_COLUMNS if name in by_lower), None)
        strand = next((by_lower[name] for name in cls._STRAND_COLUMNS if name in by_lower), None)
        if start is None or end is None:
            return None
        return start, end, strand

    @staticmethod
    def _coordinates_from_row(record, coordinate_columns: tuple[str, str, str | None] | None) -> tuple[int | None, int | None, str | None]:
        """
        Reads one annotation row's gene coordinates, tolerating blanks
        and non-numeric values (either leaves the candidate without
        coordinates rather than failing the lookup).

        Args:
            record: One row of the annotation table.
            coordinate_columns (tuple | None): From _coordinate_columns().

        Returns:
            tuple[int | None, int | None, str | None]: (start, end, strand) with strand normalised to "+"/"-" or None.
        """
        if coordinate_columns is None:
            return None, None, None
        start_col, end_col, strand_col = coordinate_columns
        try:
            start = int(float(record[start_col])) if pd.notna(record[start_col]) else None
            end = int(float(record[end_col])) if pd.notna(record[end_col]) else None
        except (TypeError, ValueError):
            return None, None, None
        strand: str | None = None
        if strand_col is not None and pd.notna(record[strand_col]):
            raw = str(record[strand_col]).strip().lower()
            if raw in ("-", "-1", "minus", "reverse", "rev", "neg", "negative", "complement"):
                strand = "-"
            elif raw in ("+", "1", "+1", "plus", "forward", "fwd", "pos", "positive"):
                strand = "+"
        return start, end, strand

    def fetch_sequence(self, accession: str) -> str:
        """
        Retrieves the sequence for a MitoFish accession from the local
        sequence database. The entries in MitoFish's sequence FASTA are
        the full records (a complete mitogenome for a complete-genome
        accession), NOT the gene region the annotation row describes;
        cutting the gene out is seed_selector's job, using the
        candidate's annotation coordinates.

        Args:
            accession (str): The MitoFish accession to fetch.

        Returns:
            str: The sequence for the requested accession.
        """
        # Obtaining the logging context.
        context = self._log_context()

        # Validating the sequence database file.
        sequence_path = self.database_dir / self.SEQUENCE_FILE
        if not sequence_path.is_file():
            self.logger.error(f"{context} MitoFish sequence database not found at {sequence_path}.")
            sys.exit(1)

        # Scanning the bgzipped multi-FASTA for the requested accession.
        try:
            with gzip.open(sequence_path, "rt") as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    if record.id == accession or accession in record.id.split("|"):
                        return str(record.seq)
        except Exception as e:
            self.logger.error(f"{context} Failed to read MitoFish sequence database at {sequence_path}: {e}")
            sys.exit(1)

        # Handling an accession absent from the sequence database.
        self.logger.error(f"{context} Accession {accession} not found in MitoFish sequence database.")
        sys.exit(1)