"""
typst_binary.py

Downloads and manages a local, standalone typst binary -- the real
PDF-rendering engine pandoc uses for every sample/run report. Never
relies on $PATH or an activated conda env's own bin/ at runtime:
confirmed live, 2026-11-19, that Snakemake's own conda-env activation
for a "localrule" (run_report/sample_report both are) does not
reliably propagate the activated env's bin/ onto $PATH on every
invocation within the same scheduler session -- two consecutive real
runs, same env, same machine, one found "typst" and one didn't. A
real, project-relative, self-installed binary at a fixed absolute
path sidesteps that race entirely, same contract as this project's
other self-installing external resources (table2asn, MitoFish, the
BLAST reference database, MITOS2's Zenodo reference data): never a
manual one-off download a fresh checkout has to repeat by hand.
"""

# Imports
from __future__ import annotations
import logging
import shutil
import stat
import subprocess
import sys
import tarfile
import tempfile
import urllib.request
from pathlib import Path

# typst's own official GitHub release distribution for the Linux
# x86-64 static (musl, libc-independent) build (confirmed by a real
# download and version check: "typst 0.15.1 (6734eafd)").
DEFAULT_DOWNLOAD_URL = "https://github.com/typst/typst/releases/download/v0.15.1/typst-x86_64-unknown-linux-musl.tar.xz"


class TypstBinary:
    """
    Manages downloading and validating the local typst binary.

    Attributes:
        binary_path (Path): Where the real executable should live once prepared.
        download_url (str): typst's own real GitHub release download URL for the binary.
        logger (logging.Logger): The logger to use for logging messages.
    """

    def __init__(self, binary_path: Path, logger: logging.Logger,
                 download_url: str = DEFAULT_DOWNLOAD_URL) -> None:
        """
        Initializes the TypstBinary with the path to install to.

        Args:
            binary_path (Path): Where the real executable should live once prepared.
            logger (logging.Logger): The logger to use for logging messages.
            download_url (str): typst's own real release download URL for the binary. Defaults to the confirmed real one.

        Returns:
            None
        """
        self.binary_path = Path(binary_path)
        self.download_url = download_url
        self.logger = logger

    def is_available(self) -> bool:
        """
        Checks whether a real, executable typst binary already
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
            self.logger.error(f"typst binary validation failed: not found/executable at {self.binary_path}.")
            sys.exit(1)
        result = subprocess.run([str(self.binary_path), "--version"], capture_output = True, text = True, timeout = 10)
        if result.returncode != 0 or not result.stdout.strip():
            self.logger.error(f"typst binary validation failed: {self.binary_path} did not report a version.")
            sys.exit(1)

    def prepare(self) -> Path:
        """
        Prepares the local typst binary by downloading it from
        typst's own real release distribution and making it
        executable, reusing an already-available copy when possible.

        Args:
            None

        Returns:
            Path: The path to the prepared, validated binary.
        """
        if self.is_available():
            self.logger.info(f"typst binary is already available at {self.binary_path}.")
            return self.binary_path

        self.binary_path.parent.mkdir(parents = True, exist_ok = True)
        self.logger.info(f"Downloading typst binary from {self.download_url}.")

        with tempfile.TemporaryDirectory() as tmpdir:
            archive_path = Path(tmpdir) / "typst.tar.xz"
            self._download_file(self.download_url, archive_path)
            # The real release asset is a tarball with the binary
            # nested one level down (e.g.
            # typst-x86_64-unknown-linux-musl/typst), never the bare
            # binary itself -- confirmed directly against a real
            # download.
            with tarfile.open(archive_path, mode = "r:xz") as archive:
                member = next((m for m in archive.getmembers() if Path(m.name).name == "typst" and m.isfile()), None)
                if member is None:
                    raise FileNotFoundError(f"No 'typst' binary found inside the downloaded archive from {self.download_url}.")
                extracted = archive.extractfile(member)
                if extracted is None:
                    raise FileNotFoundError(f"Could not read the 'typst' member inside the downloaded archive from {self.download_url}.")
                with extracted, self.binary_path.open("wb") as destination:
                    shutil.copyfileobj(extracted, destination)

        self.binary_path.chmod(self.binary_path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
        self.validate()
        self.logger.info(f"typst binary is ready at {self.binary_path}.")
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
