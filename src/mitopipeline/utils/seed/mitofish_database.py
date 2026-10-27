"""
mitofish_database.py

Downloads and manages the local MitoFish reference database.
"""

# Imports
from __future__ import annotations
import json
import logging
import sys
import shutil
import tarfile
import tempfile
import urllib.request
from pathlib import Path

# MitoFish's own site builds its bulk-download links from the latest
# release's "YYYY.MM" version string (confirmed directly against the
# real site: https://mitofish.aori.u-tokyo.ac.jp/app/api/stat/latest
# returns {"version": {"species": "YYYY.MM", "sequences": "YYYY.MM"},
# ...}, and the Download page's own real download buttons build their
# href as /download/tables/{year}/{month}/tables.tar and
# /download/fullseq/{year}/{month}/mitofishdb.fa.gz from that same
# string) -- there is no stable, version-independent URL.
MITOFISH_STAT_URL = "https://mitofish.aori.u-tokyo.ac.jp/app/api/stat/latest"
MITOFISH_TABLES_URL_TEMPLATE = "https://mitofish.aori.u-tokyo.ac.jp/download/tables/{year}/{month}/tables.tar"
MITOFISH_SEQUENCES_URL_TEMPLATE = "https://mitofish.aori.u-tokyo.ac.jp/download/fullseq/{year}/{month}/mitofishdb.fa.gz"

def mitofish_urls_for_version(version: str) -> tuple[str, str]:
    """
    Builds the real MitoFish tables/sequences download URLs for a
    given "YYYY.MM" release version string.

    Args:
        version (str): The MitoFish release version, formatted "YYYY.MM" (e.g. "2026.09").

    Returns:
        tuple[str, str]: The (tables_url, sequences_url) pair for this version.

    Raises:
        ValueError: If version is not in "YYYY.MM" form.
    """
    parts = version.split(".")
    if len(parts) != 2 or not all(part.isdigit() for part in parts):
        raise ValueError(f"MitoFish version must be formatted 'YYYY.MM', got: {version!r}.")
    year, month = parts
    return (
        MITOFISH_TABLES_URL_TEMPLATE.format(year = year, month = month.zfill(2)),
        MITOFISH_SEQUENCES_URL_TEMPLATE.format(year = year, month = month.zfill(2)),
    )

def resolve_latest_mitofish_version(logger: logging.Logger, stat_url: str = MITOFISH_STAT_URL) -> str:
    """
    Queries MitoFish's own real stat API for the current release
    version, so a fresh install always fetches the latest database
    without a maintainer having to update a hardcoded version string
    every time MitoFish publishes a new monthly release.

    Args:
        logger (logging.Logger): The logger to use for logging messages.
        stat_url (str): The MitoFish stat API URL. Defaults to the real endpoint.

    Returns:
        str: The latest release version, formatted "YYYY.MM".

    Raises:
        OSError: If the stat API can't be reached.
        (KeyError, ValueError): If the response doesn't carry the expected version field.
    """
    with urllib.request.urlopen(stat_url) as response:
        payload = json.loads(response.read())
    version = payload["version"]["sequences"]
    logger.info(f"Resolved latest MitoFish release version: {version}.")
    return version

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

    # The exact 5 tables mitofish.py's own MitoFishReferenceLookup reads
    # (taxonid_name/speciesid_lineageid/taxonid_speciesid for the
    # scientific-name resolution chain, seq_taxonid/seq_annotation for
    # the per-sequence gene lookup) -- confirmed against its real
    # column access, not guessed from the archive's own larger
    # contents (which also includes seq_description/seq_heterospecific,
    # neither ever read here).
    REQUIRED_METADATA_FILES = (
        "taxonid_name.parquet",
        "speciesid_lineageid.parquet",
        "taxonid_speciesid.parquet",
        "seq_taxonid.parquet",
        "seq_annotation.parquet",
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
        missing_files = [path for path in self.required_files if not path.is_file() or path.stat().st_size == 0]

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
        # Reusing an existing database when possible.
        if self.is_available():
            self.logger.info(f"Mitofish database {self.version} is already available at {self.database_dir}")
            return self.database_dir

        # Creating the database directory.
        self.database_dir.mkdir(parents = True, exist_ok = True)

        # Logging the start of the database creation process.
        self.logger.info(f"Preparing MitoFish database {self.version} at {self.database_dir}.")

        # Downloading required resources.
        self._download_metadata()
        self._download_sequences()

        # Validating the completed database.
        self.validate()
        self.logger.info(f"Mitofish database version {self.version} is ready at {self.database_dir}.")

        # Returning the database.
        return self.database_dir

    def _download_metadata(self) -> None:
        """
        Downloads and extracts the required Mitofish metadata files. 
        
        Args:
            None
        
        Returns:
            None
        """
        # Logging and downloading the temporary MitoFish table tarball. 
        self.logger.info(f"Downloading MitoFish metadata for version {self.version}.")
        with tempfile.TemporaryDirectory() as tmpdir:
            # Downloading the MitoFish table tarball.
            temp_dir = Path(tmpdir)
            archive_path = temp_dir / "tables.tar"
            self._download_file(url = self.tables_url, destination = archive_path)

            # Extracting the MitoFish table tarball.
            with tarfile.open(archive_path, "r") as archive:
                archive.extractall(path = temp_dir, filter = "data")

            # Copying the required metadata files.
            for filename in self.REQUIRED_METADATA_FILES:
                destination = self.database_dir / filename

                # Skipping if the destination file already exists and is not empty.
                if destination.is_file() and destination.stat().st_size > 0:
                    continue

                # Identifying the file to copy.
                matches = list(temp_dir.rglob(filename))

                # Validating the number of matches.
                if len(matches) != 1:
                    self.logger.error(f"Expected one match for {filename} in the extracted archive, found {len(matches)}.")
                    sys.exit(1)

                # Copying the file.
                shutil.copy2(matches[0], destination)

                # Logging the successful copy.
                self.logger.info(f"Successfully copied {filename} to {destination}.")

        # Logging the successful download.
        self.logger.info(f"Successfully downloaded MitoFish metadata for version {self.version}.")

    def _download_sequences(self) -> None:
        """
        Downloads the MitoFish sequence FASTA.
        
        Args:
            None
        
        Returns:
            None
        """
        # Setting the destination for the downloaded file.
        destination = self.database_dir / self.SEQUENCE_FILE

        # Skipping download if the destination file already exists and is not empty.
        if destination.is_file() and destination.stat().st_size > 0: return

        # Logging the start of the download process.
        self.logger.info(f"Downloading MitoFish sequence database for version {self.version}.")

        # Downloading the sequence FASTA.
        self._download_file(url = self.sequences_url, destination = destination)

        # Logging the successful download.
        self.logger.info(f"Successfully downloaded MitoFish sequence database for version {self.version}.")

    def _download_file(self, url: str, destination: Path) -> None:
        """
        Downloads a file from the specified URL to the specified destination.

        Args:
            url (str): The URL of the file to download.
            destination (Path): The destination path to save the downloaded file.

        Returns:
            None
        """
        # Setting up a temporary file for the download.
        temporary_path = destination.with_suffix(destination.suffix + ".tmp")

        # Attempting the download and handling exceptions.
        try:
            # Downloading the file. 
            with urllib.request.urlopen(url) as response:
                # Saving the downloaded file.
                with temporary_path.open("wb") as output:
                    shutil.copyfileobj(response, output)
            # Replacing the temporary file with the destination file.
            temporary_path.replace(destination)
        except Exception as e:
            #If the temporary file exists, delete it.
            if temporary_path.exists(): temporary_path.unlink()
            # Logging the error and exiting.
            self.logger.error(f"Failed to download MitoFish resource from {url}: {e}")
            sys.exit(1)

        # Logging the successful download.
        self.logger.info(f"Successfully downloaded MitoFish resource from {url} to {destination}.")