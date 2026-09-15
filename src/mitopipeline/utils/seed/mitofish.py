"""
mitofish.py

Lookup of mitochondrial refernece candidates from MitoFish.
"""

# Imports
from pathlib import Path
from mitopipeline.utils.seed.models import SeedReferenceCandidate

class MitoFishReferenceLookup:
    """
    Queries a local MitoFish reference database. 
    """
    def __init__(self, database_dir: Path) -> None:
        """
        Initializes the MitoFishReferenceLookup with the path to the local database.

        Args:
            database_dir (Path): The path to the local MitoFish reference database.
        
        Returns:
            None
        """
        self.database_dir = database_dir

    def lookup(self, scientific_name: str, gene: str = "12S") -> list[SeedReferenceCandidate]:
        """
        Returns the candidate references for an organism and gene.

        Args:
            scientific_name (str): The scientific name of the organism.
            gene (str): The mitochondrial gene to query. Defaults to "12S".

        Returns:
            list[SeedReferenceCandidate]: A list of candidate reference sequences. 
        """
        raise NotImplementedError("MitoFish reference lookup is not yet implemented.")