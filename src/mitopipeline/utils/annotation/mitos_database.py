"""
mitos_database.py

Downloads and manages the local MITOS2 reference data (BLAST/HMM/
covariance-model databases MITOS2 itself reads via -R/-r), published
on Zenodo by MITOS2's own authors.
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

ZENODO_RECORD_API_TEMPLATE = "https://zenodo.org/api/records/{record_id}"

def resolve_mitos_download_url(record_id: str, refseqver: str, logger: logging.Logger,
                                 api_template: str = ZENODO_RECORD_API_TEMPLATE) -> str:
    """
    Queries Zenodo's own real API for one record and finds the direct
    download URL for a specific reference version's archive -- never
    a hardcoded URL, since Zenodo's own file storage paths are not
    predictable from the record id and version string alone.

    Args:
        record_id (str): The Zenodo record id (e.g. "4284483", MITOS2's own published reference data).
        refseqver (str): The reference version to resolve, matching one of the record's own archive names (e.g. "refseq89m" for "refseq89m.tar.bz2").
        logger (logging.Logger): The logger to use for logging messages.
        api_template (str): The Zenodo record API URL template. Defaults to the real endpoint.

    Returns:
        str: The direct download URL for refseqver's own archive.

    Raises:
        OSError: If the Zenodo API can't be reached.
        ValueError: If the record has no file matching refseqver.
    """
    url = api_template.format(record_id = record_id)
    with urllib.request.urlopen(url) as response:
        payload = json.loads(response.read())

    filename = f"{refseqver}.tar.bz2"
    for file_entry in payload.get("files", []):
        if file_entry.get("key") == filename:
            download_url = file_entry["links"]["self"]
            logger.info(f"Resolved MITOS2 reference data download URL for {refseqver}: {download_url}.")
            return download_url

    available = ", ".join(sorted(f.get("key", "") for f in payload.get("files", [])))
    raise ValueError(f"Zenodo record {record_id} has no file named {filename!r}. Available files: {available}.")

class MitosReferenceData:
    """
    Manages downloading and validating one refseqver's local MITOS2
    reference data.

    Attributes:
        refdir (Path): The directory MITOS2's own -R flag points at (holds one subdirectory per refseqver).
        refseqver (str): The reference version to ensure is available (e.g. "refseq89m").
        record_id (str): The Zenodo record id to download from.
        logger (logging.Logger): The logger to use for logging messages.
        version_dir (Path): refdir/refseqver -- the exact directory MITOS2's own -r flag selects.
    """

    # Confirmed directly against a real, already-extracted refseq89m
    # archive -- every real refseqver archive shares this same
    # top-level layout (MITOS2's own published structure, not specific
    # to one version).
    REQUIRED_SUBDIRS = ("featureNuc", "featureProt", "ncRNA")
    REQUIRED_FILES = ("auxinfo.json",)

    def __init__(self, refdir: Path, refseqver: str, record_id: str, logger: logging.Logger) -> None:
        """
        Initializes the MitosReferenceData with the path to the local
        reference directory and which version to download.

        Args:
            refdir (Path): The directory MITOS2's own -R flag points at.
            refseqver (str): The reference version to ensure is available.
            record_id (str): The Zenodo record id to download from.
            logger (logging.Logger): The logger to use for logging messages.

        Returns:
            None
        """
        self.refdir = Path(refdir)
        self.refseqver = refseqver
        self.record_id = record_id
        self.logger = logger
        self.version_dir = self.refdir / self.refseqver

    def is_available(self) -> bool:
        """
        Checks whether this refseqver's required reference data
        already exists under refdir.

        Args:
            None

        Returns:
            bool: True if every required subdirectory/file exists, False otherwise.
        """
        if not all((self.version_dir / name).is_dir() for name in self.REQUIRED_SUBDIRS):
            return False
        return all((self.version_dir / name).is_file() and (self.version_dir / name).stat().st_size > 0
                   for name in self.REQUIRED_FILES)

    def validate(self) -> None:
        """
        Validates the local reference data by checking for the
        existence of its required subdirectories/files.

        Args:
            None

        Returns:
            None
        """
        if not self.is_available():
            self.logger.error(f"MITOS2 reference data validation failed for {self.refseqver} at {self.version_dir}.")
            sys.exit(1)

    def prepare(self) -> Path:
        """
        Prepares this refseqver's local reference data by downloading
        and extracting it from Zenodo, reusing an already-available
        copy when possible.

        Args:
            None

        Returns:
            Path: The path to the prepared reference directory (refdir/refseqver).
        """
        if self.is_available():
            self.logger.info(f"MITOS2 reference data {self.refseqver} is already available at {self.version_dir}.")
            return self.version_dir

        self.refdir.mkdir(parents = True, exist_ok = True)
        self.logger.info(f"Preparing MITOS2 reference data {self.refseqver} at {self.version_dir}.")

        download_url = resolve_mitos_download_url(self.record_id, self.refseqver, self.logger)
        with tempfile.TemporaryDirectory() as tmpdir:
            archive_path = Path(tmpdir) / f"{self.refseqver}.tar.bz2"
            self._download_file(download_url, archive_path)

            # The real archive's own top-level entry is refseqver/
            # itself (confirmed directly against the real
            # refseq89m.tar.bz2), so extracting straight into refdir
            # reproduces refdir/refseqver/... with no renaming needed.
            with tarfile.open(archive_path, "r:bz2") as archive:
                archive.extractall(path = self.refdir, filter = "data")

        self.validate()
        self.logger.info(f"MITOS2 reference data {self.refseqver} is ready at {self.version_dir}.")
        return self.version_dir

    def _download_file(self, url: str, destination: Path) -> None:
        """
        Downloads a file from the specified URL to the specified destination.

        Args:
            url (str): The URL of the file to download.
            destination (Path): The destination path to save the downloaded file.

        Returns:
            None
        """
        temporary_path = destination.with_suffix(destination.suffix + ".tmp")
        try:
            with urllib.request.urlopen(url) as response:
                with temporary_path.open("wb") as output:
                    shutil.copyfileobj(response, output)
            temporary_path.replace(destination)
        except Exception as e:
            if temporary_path.exists(): temporary_path.unlink()
            self.logger.error(f"Failed to download MITOS2 reference data from {url}: {e}")
            sys.exit(1)

        self.logger.info(f"Successfully downloaded MITOS2 reference data from {url} to {destination}.")
