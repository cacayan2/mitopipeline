"""
test_sra_metadata.py

Unit tests for SRA metadata generation (#368/#371/#447): manifest
values win over configured defaults, missing required fields are
tracked rather than fabricated, and paired-layout samples require a
real R2 file.
"""

# Imports
from pathlib import Path
import pytest
from mitopipeline.utils.archival.sra_metadata import build_sra_metadata, write_sra_metadata

DEFAULTS = {
    "library_strategy": "WGS", "library_source": "GENOMIC", "library_selection": "RANDOM",
    "library_layout": "PAIRED", "platform": "ILLUMINA", "instrument_model": "Illumina NovaSeq 6000",
    "filetype": "fastq", "title_template": "{sample_id} raw sequencing reads",
}

def _write_manifest(tmp_path: Path, rows: list[str], header: str = "sample_id\tr1\tr2\tgenus\tspecies") -> Path:
    r1 = tmp_path / "s1_R1.fastq.gz"
    r2 = tmp_path / "s1_R2.fastq.gz"
    r1.write_bytes(b"")
    r2.write_bytes(b"")
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text(header + "\n" + "\n".join(rows) + "\n")
    return manifest

def test_build_sra_metadata_uses_defaults_when_manifest_has_no_value(tmp_path: Path) -> None:
    manifest = _write_manifest(tmp_path, [f"sample1\t{tmp_path / 's1_R1.fastq.gz'}\t{tmp_path / 's1_R2.fastq.gz'}\tDanio\trerio"])

    records = build_sra_metadata(manifest, defaults = DEFAULTS)

    assert len(records) == 1
    record = records[0]
    assert record["sample_name"] == "sample1"
    assert record["library_strategy"] == "WGS"
    assert record["platform"] == "ILLUMINA"
    assert record["organism"] == "Danio rerio"
    assert record["review_status"] == "READY"
    assert record["missing_fields"] == ""

def test_manifest_values_win_over_defaults(tmp_path: Path) -> None:
    manifest = _write_manifest(
        tmp_path,
        [f"sample1\t{tmp_path / 's1_R1.fastq.gz'}\t{tmp_path / 's1_R2.fastq.gz'}\tDanio\trerio\tOXFORD_NANOPORE"],
        header = "sample_id\tr1\tr2\tgenus\tspecies\tplatform",
    )
    records = build_sra_metadata(manifest, defaults = DEFAULTS)
    assert records[0]["platform"] == "OXFORD_NANOPORE"

def test_missing_required_field_is_tracked_not_fabricated(tmp_path: Path) -> None:
    manifest = _write_manifest(tmp_path, [f"sample1\t{tmp_path / 's1_R1.fastq.gz'}\t{tmp_path / 's1_R2.fastq.gz'}\tDanio\trerio"])
    records = build_sra_metadata(manifest, defaults = {})

    record = records[0]
    assert record["review_status"] == "MANUAL_REVIEW_REQUIRED"
    assert "library_strategy" in record["missing_fields"]
    assert "platform" in record["missing_fields"]

def test_paired_layout_without_r2_is_flagged(tmp_path: Path) -> None:
    r1 = tmp_path / "single_R1.fastq.gz"
    r1.write_bytes(b"")
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text(f"sample_id\tr1\tr2\nsingle_sample\t{r1}\t\n")

    records = build_sra_metadata(manifest, defaults = DEFAULTS)
    assert "filename2" in records[0]["missing_fields"]
    assert records[0]["review_status"] == "MANUAL_REVIEW_REQUIRED"

def test_blank_sample_id_raises(tmp_path: Path) -> None:
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text("sample_id\tr1\tr2\n\ta\tb\n")
    with pytest.raises(ValueError):
        build_sra_metadata(manifest)

def test_missing_manifest_raises(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        build_sra_metadata(tmp_path / "missing.tsv")

def test_write_sra_metadata_writes_a_real_tsv(tmp_path: Path) -> None:
    manifest = _write_manifest(tmp_path, [f"sample1\t{tmp_path / 's1_R1.fastq.gz'}\t{tmp_path / 's1_R2.fastq.gz'}\tDanio\trerio"])
    output_path = tmp_path / "sra_metadata.tsv"

    records = write_sra_metadata(manifest, output_path, defaults = DEFAULTS)

    assert output_path.is_file()
    text = output_path.read_text()
    assert "sample_name" in text.splitlines()[0]
    assert len(records) == 1
