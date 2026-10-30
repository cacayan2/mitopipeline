"""
test_table2asn_binary.py

Unit tests for the table2asn binary auto-installer: availability/
validation checks, and a real download+gunzip+chmod round trip
against a faked (never real network) response.
"""

# Imports
import gzip
import logging
from pathlib import Path
from unittest.mock import patch
import pytest
from mitopipeline.utils.archival.table2asn_binary import Table2AsnBinary

@pytest.fixture
def binary(tmp_path: Path) -> Table2AsnBinary:
    return Table2AsnBinary(binary_path = tmp_path / "bin" / "table2asn", logger = logging.getLogger("test"))

def test_is_available_false_when_nothing_downloaded(binary: Table2AsnBinary) -> None:
    assert binary.is_available() is False

def test_is_available_true_for_a_real_executable_file(binary: Table2AsnBinary) -> None:
    binary.binary_path.parent.mkdir(parents = True)
    binary.binary_path.write_bytes(b"fake binary")
    binary.binary_path.chmod(0o755)
    assert binary.is_available() is True

def test_validate_exits_when_not_available(binary: Table2AsnBinary) -> None:
    with pytest.raises(SystemExit):
        binary.validate()

def _fake_gz_bytes(content: bytes) -> bytes:
    import io
    buf = io.BytesIO()
    with gzip.GzipFile(fileobj = buf, mode = "wb") as handle:
        handle.write(content)
    return buf.getvalue()

def test_prepare_downloads_gunzips_and_makes_executable(tmp_path: Path, binary: Table2AsnBinary) -> None:
    # A tiny real shell script standing in for the real binary --
    # this proves the download/gunzip/chmod/validate pipeline end to
    # end without a real network call or the real multi-MB binary.
    fake_binary_content = b"#!/bin/sh\necho 'table2asn: 9.9.9'\n"
    compressed = _fake_gz_bytes(fake_binary_content)

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

    with patch("mitopipeline.utils.archival.table2asn_binary.urllib.request.urlopen",
               return_value = _FakeResponse(compressed)):
        result_path = binary.prepare()

    assert result_path == binary.binary_path
    assert binary.binary_path.is_file()
    assert binary.binary_path.read_bytes() == fake_binary_content
    assert binary.is_available() is True

def test_prepare_reuses_an_already_available_binary(tmp_path: Path, binary: Table2AsnBinary) -> None:
    binary.binary_path.parent.mkdir(parents = True)
    binary.binary_path.write_text("#!/bin/sh\necho 'table2asn: 1.0.0'\n")
    binary.binary_path.chmod(0o755)

    with patch("mitopipeline.utils.archival.table2asn_binary.urllib.request.urlopen") as mock_urlopen:
        binary.prepare()

    mock_urlopen.assert_not_called()
