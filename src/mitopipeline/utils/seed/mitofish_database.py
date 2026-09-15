"""
mitofish_database.py

Downloads and manages the local MitoFish reference database.
"""

# Imports
from __future__ import annotations
import logging
import sys
import shutil
import tarfile
import tempfile
import urllib.request
from pathlib import Path

class MitoFishDatabase:
    """
    Manages downloading and validating a local copy of the MitoFish reference database.

    Attributes:
        database_dir (Path): The path to the local MitoFish reference database.
        version (str): The version of the MitoFish database to download.
        tables_url (str): The URL to download the MitoFish database tables.
        sequences_url (str): The URL to download the MitoFish database sequences.
        logger (logging.Logger): The logger to use for logging messages.
    """

    REQUIRED_METADATA_FILES = (
        "taxonid_name.parquet",
        "seq_taxonid.parquet",
        "seq_annotation.parquet"
    )

    SEQUENCE_FILE = "mitofishdb.fa.gz"

    def __init__(self, database_dir: Path, version: str, tables_url: str, sequences_url: str, logger: logging.Logger) -> None:
        """
        Initializes the MitoFishDatabase with the path to the local database and the version to download.

        Args:
            database_dir (Path): The path to the local MitoFish reference database.
            version (str): The version of the MitoFish database to download.
            tables_url (str): The URL to download the MitoFish database tables.
            sequences_url (str): The URL to download the MitoFish database sequences.
            logger (logging.Logger): The logger to use for logging messages.

        Returns:
            None
        """
        # Storing member variables for class. 
        self.database_root = Path(database_dir)
        self.version = version
        self.tables_url = tables_url
        self.sequences_url = sequences_url
        self.logger = logger

        # Defining the verison-specific database directory.
        self.database_dir = self.database_root / self.version

    @property
    def required_files(self) -> tuple[Path, ...]:
        """
        Returns the required files for the MitoFish database.

        Returns:
            tuple[Path, ...]: A tuple of Paths to the required files.
        """
        # Defining the required metadata files and sequence file.
        metadata_files = tuple(self.database_dir / filename for filename in self.REQUIRED_METADATA_FILES)

        # Returning the required files.
        return (*metadata_files, self.database_dir / self.SEQUENCE_FILE)

    def is_available(self) -> bool:
        """
        Checks whether the required MitoFish database files exist.
        
        Args:
            None
        
        Returns:
            bool: True if all required files exist, False otherwise.
        """
        return all(path.is_file() and path.stat().st_size > 0 for path in self.required_files)

    def validate(self) -> None:
        """
        Validates the MitoFish database by checking for the existence of required files.

        Args:
            None

        Returns:
            None
        """
        # Identifying missing or empty database files.
        missing_files = [path for path in self.required files if not path.is_file() or path.stat().st_size == 0]

        # Logging an error and exiting if any required files are missing or empty.
        if missing_files:
            missing = ", ".join(str(path) for path in missing_files)
            self.logger.error(f"MitoFish database validation failed. Missing or empty files: {missing}")
            sys.exit(1)

    def prepare(self) -> Path:
        """
        Prepares the MitoFish database by downloading and extracting the required files.

        Args:
            None

        Returns:
            Path: The path to the prepared MitoFish database.
        """
        # Reusing 

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