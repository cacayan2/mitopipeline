"""
test_mitofish_database.py

Unit tests for the MitoFish database download/setup machinery:
version resolution, URL construction, and MitoFishDatabase's own
availability/validation/preparation behavior.

Regression context: REQUIRED_METADATA_FILES used to list only 3 of
the 5 real tables mitofish.py actually reads (speciesid_lineageid.parquet
and taxonid_speciesid.parquet were missing), so even a successful
download silently left seed selection's own name-to-accession
resolution chain unable to run -- confirmed against a real prior
incident where those exact 2 files had to be copied in by hand.
"""

# Imports
import json
import logging
import tarfile
from pathlib import Path
from unittest.mock import patch
import pytest
from mitopipeline.utils.seed.mitofish_database import (
    MitoFishDatabase, mitofish_urls_for_version, resolve_latest_mitofish_version,
)

def test_required_metadata_files_covers_every_table_mitofish_py_reads() -> None:
    """Regression: this exact 5-file set is what mitofish.py's
    MitoFishReferenceLookup reads from database_dir -- confirmed
    directly against its own real column access, not guessed."""
    assert set(MitoFishDatabase.REQUIRED_METADATA_FILES) == {
        "taxonid_name.parquet",
        "speciesid_lineageid.parquet",
        "taxonid_speciesid.parquet",
        "seq_taxonid.parquet",
        "seq_annotation.parquet",
    }

def test_mitofish_urls_for_version_builds_the_real_url_pattern() -> None:
    tables_url, sequences_url = mitofish_urls_for_version("2026.09")
    assert tables_url == "https://mitofish.aori.u-tokyo.ac.jp/download/tables/2026/09/tables.tar"
    assert sequences_url == "https://mitofish.aori.u-tokyo.ac.jp/download/fullseq/2026/09/mitofishdb.fa.gz"

def test_mitofish_urls_for_version_zero_pads_a_single_digit_month() -> None:
    tables_url, _ = mitofish_urls_for_version("2026.9")
    assert "/2026/09/" in tables_url

def test_mitofish_urls_for_version_rejects_a_malformed_version() -> None:
    with pytest.raises(ValueError):
        mitofish_urls_for_version("not-a-version")
    with pytest.raises(ValueError):
        mitofish_urls_for_version("2026")

def test_resolve_latest_mitofish_version_reads_the_real_stat_api_shape() -> None:
    """Confirmed directly against the real response shape from
    https://mitofish.aori.u-tokyo.ac.jp/app/api/stat/latest."""
    payload = json.dumps({"version": {"species": "2026.09", "sequences": "2026.09"}}).encode()

    class _FakeResponse:
        def __enter__(self): return self
        def __exit__(self, *a): return False
        def read(self): return payload

    with patch("mitopipeline.utils.seed.mitofish_database.urllib.request.urlopen", return_value = _FakeResponse()):
        version = resolve_latest_mitofish_version(logging.getLogger("test"))
    assert version == "2026.09"

def _make_tables_tar(path: Path) -> None:
    """Builds a real tarball containing all 7 real MitoFish table
    files (the 5 required ones plus 2 real ones never required),
    matching the real archive's own contents."""
    with tarfile.open(path, "w") as archive:
        for filename in (*MitoFishDatabase.REQUIRED_METADATA_FILES, "seq_description.parquet", "seq_heterospecific.parquet"):
            source = path.parent / filename
            source.write_text(f"contents of {filename}")
            archive.add(source, arcname = filename)

@pytest.fixture
def database(tmp_path: Path) -> MitoFishDatabase:
    return MitoFishDatabase(
        database_dir = tmp_path / "mitofish", version = "2026.09",
        tables_url = "https://example.invalid/tables.tar",
        sequences_url = "https://example.invalid/mitofishdb.fa.gz",
        logger = logging.getLogger("test"),
    )

def test_is_available_false_when_nothing_downloaded(database: MitoFishDatabase) -> None:
    assert database.is_available() is False

def test_prepare_skips_download_when_already_available(database: MitoFishDatabase) -> None:
    database.database_dir.mkdir(parents = True)
    for filename in database.required_files:
        filename.write_text("already here")

    with patch.object(database, "_download_metadata") as mock_metadata, \
         patch.object(database, "_download_sequences") as mock_sequences:
        result = database.prepare()
        mock_metadata.assert_not_called()
        mock_sequences.assert_not_called()
    assert result == database.database_dir
    assert database.is_available() is True

def test_prepare_downloads_and_copies_every_required_table(tmp_path: Path, database: MitoFishDatabase) -> None:
    """Regression: a real download must leave all 5 required tables
    in database_dir, not just the 3 the old REQUIRED_METADATA_FILES
    bug would have copied."""
    tar_path = tmp_path / "tables.tar"
    _make_tables_tar(tar_path)

    def _fake_download_file(url: str, destination: Path) -> None:
        if destination.name == "tables.tar":
            destination.write_bytes(tar_path.read_bytes())
        else:
            destination.write_text("fake sequence data")

    with patch.object(database, "_download_file", side_effect = _fake_download_file):
        result = database.prepare()

    for filename in MitoFishDatabase.REQUIRED_METADATA_FILES:
        assert (result / filename).is_file()
    assert (result / MitoFishDatabase.SEQUENCE_FILE).is_file()
    assert database.is_available() is True

def test_prepare_exits_when_a_required_table_is_missing_from_the_archive(tmp_path: Path, database: MitoFishDatabase) -> None:
    """A real-world case this guards against: MitoFish changing its
    archive layout, or a partial/corrupt download -- validate() must
    catch it rather than silently leaving seed selection to fail
    later with a confusing file-not-found."""
    incomplete_tar = tmp_path / "incomplete.tar"
    with tarfile.open(incomplete_tar, "w") as archive:
        source = tmp_path / "taxonid_name.parquet"
        source.write_text("contents")
        archive.add(source, arcname = "taxonid_name.parquet")

    def _fake_download_file(url: str, destination: Path) -> None:
        if destination.name == "tables.tar":
            destination.write_bytes(incomplete_tar.read_bytes())
        else:
            destination.write_text("fake sequence data")

    with patch.object(database, "_download_file", side_effect = _fake_download_file):
        with pytest.raises(SystemExit):
            database.prepare()
