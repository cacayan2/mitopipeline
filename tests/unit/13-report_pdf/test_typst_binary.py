"""
test_typst_binary.py

Unit tests for the typst binary auto-installer: availability/
validation checks, and a real download+tar.xz-extract+chmod round
trip against a faked (never real network) response.
"""

# Imports
import io
import logging
import tarfile
from pathlib import Path
from unittest.mock import patch
import pytest
from mitopipeline.utils.reporting.typst_binary import TypstBinary

@pytest.fixture
def binary(tmp_path: Path) -> TypstBinary:
    return TypstBinary(binary_path = tmp_path / "bin" / "typst", logger = logging.getLogger("test"))

def test_is_available_false_when_nothing_downloaded(binary: TypstBinary) -> None:
    assert binary.is_available() is False

def test_is_available_true_for_a_real_executable_file(binary: TypstBinary) -> None:
    binary.binary_path.parent.mkdir(parents = True)
    binary.binary_path.write_bytes(b"fake binary")
    binary.binary_path.chmod(0o755)
    assert binary.is_available() is True

def test_validate_exits_when_not_available(binary: TypstBinary) -> None:
    with pytest.raises(SystemExit):
        binary.validate()

def _fake_release_tarball_bytes(content: bytes) -> bytes:
    # Mirrors the real release asset's own layout: the binary sits
    # one level down inside a named directory, alongside README/LICENSE
    # files -- never the bare binary itself (confirmed directly
    # against a real download).
    buf = io.BytesIO()
    with tarfile.open(fileobj = buf, mode = "w:xz") as archive:
        info = tarfile.TarInfo(name = "typst-x86_64-unknown-linux-musl/typst")
        info.size = len(content)
        archive.addfile(info, io.BytesIO(content))
        for extra in ("README.md", "LICENSE"):
            extra_info = tarfile.TarInfo(name = f"typst-x86_64-unknown-linux-musl/{extra}")
            extra_info.size = 0
            archive.addfile(extra_info, io.BytesIO(b""))
    return buf.getvalue()

def test_prepare_downloads_extracts_and_makes_executable(tmp_path: Path, binary: TypstBinary) -> None:
    # A tiny real shell script standing in for the real binary --
    # this proves the download/extract/chmod/validate pipeline end to
    # end without a real network call or the real multi-MB binary.
    fake_binary_content = b"#!/bin/sh\necho 'typst 9.9.9'\n"
    archive_bytes = _fake_release_tarball_bytes(fake_binary_content)

    class _FakeResponse:
        def __init__(self, data: bytes) -> None:
            self._data, self._pos = data, 0

        def __enter__(self): return self
        def __exit__(self, *a): return False

        def read(self, n = -1):
            # shutil.copyfileobj calls read(bufsize) repeatedly until it
            # gets an empty result -- a fake that doesn't advance its
            # own position (returning the same prefix every call) loops
            # it forever instead of a real stream's natural end.
            if n in (-1, None):
                chunk, self._pos = self._data[self._pos:], len(self._data)
            else:
                chunk = self._data[self._pos:self._pos + n]
                self._pos += len(chunk)
            return chunk

    with patch("mitopipeline.utils.reporting.typst_binary.urllib.request.urlopen",
               return_value = _FakeResponse(archive_bytes)):
        result_path = binary.prepare()

    assert result_path == binary.binary_path
    assert binary.binary_path.is_file()
    assert binary.binary_path.read_bytes() == fake_binary_content
    assert binary.is_available() is True

def test_prepare_reuses_an_already_available_binary(tmp_path: Path, binary: TypstBinary) -> None:
    binary.binary_path.parent.mkdir(parents = True)
    binary.binary_path.write_text("#!/bin/sh\necho 'typst 1.0.0'\n")
    binary.binary_path.chmod(0o755)

    with patch("mitopipeline.utils.reporting.typst_binary.urllib.request.urlopen") as mock_urlopen:
        binary.prepare()

    mock_urlopen.assert_not_called()

def test_prepare_raises_when_archive_has_no_typst_member(tmp_path: Path, binary: TypstBinary) -> None:
    buf = io.BytesIO()
    with tarfile.open(fileobj = buf, mode = "w:xz") as archive:
        info = tarfile.TarInfo(name = "typst-x86_64-unknown-linux-musl/README.md")
        info.size = 0
        archive.addfile(info, io.BytesIO(b""))
    archive_bytes = buf.getvalue()

    class _FakeResponse:
        def __init__(self, data: bytes) -> None:
            self._data, self._pos = data, 0
        def __enter__(self): return self
        def __exit__(self, *a): return False
        def read(self, n = -1):
            if n in (-1, None):
                chunk, self._pos = self._data[self._pos:], len(self._data)
            else:
                chunk = self._data[self._pos:self._pos + n]
                self._pos += len(chunk)
            return chunk

    with patch("mitopipeline.utils.reporting.typst_binary.urllib.request.urlopen",
               return_value = _FakeResponse(archive_bytes)):
        with pytest.raises(FileNotFoundError):
            binary.prepare()
