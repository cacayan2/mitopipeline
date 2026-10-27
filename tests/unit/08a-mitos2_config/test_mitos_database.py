"""
test_mitos_database.py

Unit tests for the MITOS2 reference data download/setup machinery:
Zenodo URL resolution and MitosReferenceData's own availability/
validation/preparation behavior.
"""

# Imports
import json
import logging
import tarfile
from pathlib import Path
from unittest.mock import patch
import pytest
from mitopipeline.utils.annotation.mitos_database import MitosReferenceData, resolve_mitos_download_url

def test_resolve_mitos_download_url_finds_the_matching_file() -> None:
    """Confirmed directly against the real response shape from
    https://zenodo.org/api/records/4284483."""
    payload = json.dumps({
        "files": [
            {"key": "refseq39.tar.bz2", "links": {"self": "https://zenodo.org/.../refseq39.tar.bz2/content"}},
            {"key": "refseq89m.tar.bz2", "links": {"self": "https://zenodo.org/.../refseq89m.tar.bz2/content"}},
        ],
    }).encode()

    class _FakeResponse:
        def __enter__(self): return self
        def __exit__(self, *a): return False
        def read(self): return payload

    with patch("mitopipeline.utils.annotation.mitos_database.urllib.request.urlopen", return_value = _FakeResponse()):
        url = resolve_mitos_download_url("4284483", "refseq89m", logging.getLogger("test"))
    assert url == "https://zenodo.org/.../refseq89m.tar.bz2/content"

def test_resolve_mitos_download_url_raises_when_version_not_found() -> None:
    payload = json.dumps({"files": [{"key": "refseq39.tar.bz2", "links": {"self": "x"}}]}).encode()

    class _FakeResponse:
        def __enter__(self): return self
        def __exit__(self, *a): return False
        def read(self): return payload

    with patch("mitopipeline.utils.annotation.mitos_database.urllib.request.urlopen", return_value = _FakeResponse()):
        with pytest.raises(ValueError):
            resolve_mitos_download_url("4284483", "refseq999x", logging.getLogger("test"))

@pytest.fixture
def reference_data(tmp_path: Path) -> MitosReferenceData:
    return MitosReferenceData(
        refdir = tmp_path / "mitos", refseqver = "refseq89m",
        record_id = "4284483", logger = logging.getLogger("test"),
    )

def test_is_available_false_when_nothing_downloaded(reference_data: MitosReferenceData) -> None:
    assert reference_data.is_available() is False

def _populate_real_structure(version_dir: Path) -> None:
    for subdir in MitosReferenceData.REQUIRED_SUBDIRS:
        (version_dir / subdir).mkdir(parents = True)
    for filename in MitosReferenceData.REQUIRED_FILES:
        (version_dir / filename).write_text("{}")

def test_prepare_skips_download_when_already_available(reference_data: MitosReferenceData) -> None:
    _populate_real_structure(reference_data.version_dir)

    with patch.object(reference_data, "_download_file") as mock_download:
        result = reference_data.prepare()
        mock_download.assert_not_called()
    assert result == reference_data.version_dir

def _make_real_archive(path: Path, refseqver: str) -> None:
    """Builds a real tar.bz2 archive matching the real MITOS2 Zenodo
    archive's own layout -- refseqver/ as its own top-level entry."""
    staging = path.parent / "staging"
    version_dir = staging / refseqver
    _populate_real_structure(version_dir)
    with tarfile.open(path, "w:bz2") as archive:
        archive.add(version_dir, arcname = refseqver)

def test_prepare_downloads_and_extracts_the_real_archive_layout(tmp_path: Path, reference_data: MitosReferenceData) -> None:
    archive_path = tmp_path / "refseq89m.tar.bz2"
    _make_real_archive(archive_path, "refseq89m")

    def _fake_download_file(url: str, destination: Path) -> None:
        destination.write_bytes(archive_path.read_bytes())

    with patch("mitopipeline.utils.annotation.mitos_database.resolve_mitos_download_url", return_value = "https://example.invalid/refseq89m.tar.bz2"), \
         patch.object(reference_data, "_download_file", side_effect = _fake_download_file):
        result = reference_data.prepare()

    assert result == reference_data.version_dir
    for subdir in MitosReferenceData.REQUIRED_SUBDIRS:
        assert (result / subdir).is_dir()
    for filename in MitosReferenceData.REQUIRED_FILES:
        assert (result / filename).is_file()
    assert reference_data.is_available() is True

def test_prepare_exits_when_the_archive_is_missing_a_required_subdir(tmp_path: Path, reference_data: MitosReferenceData) -> None:
    incomplete_dir = tmp_path / "incomplete" / "refseq89m"
    incomplete_dir.mkdir(parents = True)
    (incomplete_dir / "auxinfo.json").write_text("{}")   # missing featureNuc/featureProt/ncRNA
    archive_path = tmp_path / "refseq89m.tar.bz2"
    with tarfile.open(archive_path, "w:bz2") as archive:
        archive.add(incomplete_dir, arcname = "refseq89m")

    def _fake_download_file(url: str, destination: Path) -> None:
        destination.write_bytes(archive_path.read_bytes())

    with patch("mitopipeline.utils.annotation.mitos_database.resolve_mitos_download_url", return_value = "https://example.invalid/refseq89m.tar.bz2"), \
         patch.object(reference_data, "_download_file", side_effect = _fake_download_file):
        with pytest.raises(SystemExit):
            reference_data.prepare()
