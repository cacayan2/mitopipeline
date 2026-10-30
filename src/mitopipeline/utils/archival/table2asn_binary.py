"""
table2asn_binary.py

Downloads and manages the local table2asn binary -- NCBI's own real
GenBank submission-validation tool, published as a static binary on
NCBI's public FTP, not packaged for this system (confirmed directly:
neither apt nor a conda env exists here for it). Same contract as
this project's other self-installing external resources (MitoFish,
the BLAST reference database, MITOS2's Zenodo reference data): never
a manual one-off download a fresh checkout has to repeat by hand.
"""

# Imports
from __future__ import annotations
import gzip
import logging
import shutil
import stat
import subprocess
import sys
import tempfile
import urllib.request
from pathlib import Path

# NCBI's own real, stable distribution URL for the Linux x86-64
# static binary (confirmed by a real download and version check:
# "table2asn: 1.29.324").
DEFAULT_DOWNLOAD_URL = "https://ftp.ncbi.nlm.nih.gov/asn1-converters/by_program/table2asn/linux64.table2asn.gz"


class Table2AsnBinary:
    """
    Manages downloading and validating the local table2asn binary.

    Attributes:
        binary_path (Path): Where the real executable should live once prepared.
        download_url (str): NCBI's own real download URL for the binary.
        logger (logging.Logger): The logger to use for logging messages.
    """

    def __init__(self, binary_path: Path, logger: logging.Logger,
                 download_url: str = DEFAULT_DOWNLOAD_URL) -> None:
        """
        Initializes the Table2AsnBinary with the path to install to.

        Args:
            binary_path (Path): Where the real executable should live once prepared.
            logger (logging.Logger): The logger to use for logging messages.
            download_url (str): NCBI's own real download URL for the binary. Defaults to the confirmed real one.

        Returns:
            None
        """
        self.binary_path = Path(binary_path)
        self.download_url = download_url
        self.logger = logger

    def is_available(self) -> bool:
        """
        Checks whether a real, executable table2asn binary already
        exists at binary_path.

        Args:
            None

        Returns:
            bool: True if binary_path exists, is non-empty, and is executable.
        """
        return bool(self.binary_path.is_file() and self.binary_path.stat().st_size > 0
                    and self.binary_path.stat().st_mode & stat.S_IXUSR)

    def validate(self) -> None:
        """
        Validates the local binary by actually running it and
        confirming it reports a real version -- existing and
        executable isn't proof it isn't a truncated/corrupt download.

        Args:
            None

        Returns:
            None
        """
        if not self.is_available():
            self.logger.error(f"table2asn binary validation failed: not found/executable at {self.binary_path}.")
            sys.exit(1)
        result = subprocess.run([str(self.binary_path), "-version"], capture_output = True, text = True, timeout = 10)
        if result.returncode != 0 or not result.stdout.strip():
            self.logger.error(f"table2asn binary validation failed: {self.binary_path} did not report a version.")
            sys.exit(1)

    def prepare(self) -> Path:
        """
        Prepares the local table2asn binary by downloading it from
        NCBI's own real distribution and making it executable,
        reusing an already-available copy when possible.

        Args:
            None

        Returns:
            Path: The path to the prepared, validated binary.
        """
        if self.is_available():
            self.logger.info(f"table2asn binary is already available at {self.binary_path}.")
            return self.binary_path

        self.binary_path.parent.mkdir(parents = True, exist_ok = True)
        self.logger.info(f"Downloading table2asn binary from {self.download_url}.")

        with tempfile.TemporaryDirectory() as tmpdir:
            archive_path = Path(tmpdir) / "table2asn.gz"
            self._download_file(self.download_url, archive_path)
            with gzip.open(archive_path, "rb") as compressed, self.binary_path.open("wb") as destination:
                shutil.copyfileobj(compressed, destination)

        self.binary_path.chmod(self.binary_path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
        self.validate()
        self.logger.info(f"table2asn binary is ready at {self.binary_path}.")
        return self.binary_path

    def _download_file(self, url: str, destination: Path) -> None:
        """
        Downloads a file from the specified URL to the specified
        destination.

        Args:
            url (str): The URL of the file to download.
            destination (Path): The destination path to save the downloaded file.

        Returns:
            None
        """
        with urllib.request.urlopen(url, timeout = 60) as response, destination.open("wb") as handle:
            shutil.copyfileobj(response, handle)
