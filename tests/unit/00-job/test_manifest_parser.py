"""
test_manifest_parser.py

Unit tests for sample manifest parsing: required-column/duplicate-id
validation, path resolution, and that every manifest column actually
reaches the resulting Sample object -- regression coverage for a
real `speices` typo that made every manifest with a species column
raise a TypeError constructing Sample.
"""

# Imports
from pathlib import Path
import pytest
from mitopipeline.utils.manifest.manifest_parser import parse_sample_manifest

def _write_fastqs(tmp_path: Path, sample_id: str) -> tuple[Path, Path]:
    r1 = tmp_path / f"{sample_id}_R1.fastq.gz"
    r2 = tmp_path / f"{sample_id}_R2.fastq.gz"
    r1.write_bytes(b"")
    r2.write_bytes(b"")
    return r1, r2

def test_parses_genus_species_source_onto_the_sample(tmp_path: Path) -> None:
    r1, r2 = _write_fastqs(tmp_path, "sample1")
    manifest_path = tmp_path / "manifest.tsv"
    manifest_path.write_text(
        "sample_id\tr1\tr2\tgenus\tspecies\tsource\n"
        f"sample1\t{r1.name}\t{r2.name}\tDanio\trerio\tgargammel_simulated\n"
    )

    samples = parse_sample_manifest(manifest_path)
    assert len(samples) == 1
    sample = samples[0]
    assert sample.sample_id == "sample1"
    assert sample.genus == "Danio"
    assert sample.species == "rerio"
    assert sample.source == "gargammel_simulated"
    assert sample.r1 == r1
    assert sample.r2 == r2

def test_missing_required_column_raises(tmp_path: Path) -> None:
    manifest_path = tmp_path / "manifest.tsv"
    manifest_path.write_text("sample_id\tr1\n" "sample1\tfoo.fastq.gz\n")

    with pytest.raises(ValueError):
        parse_sample_manifest(manifest_path)

def test_duplicate_sample_ids_raise(tmp_path: Path) -> None:
    r1, r2 = _write_fastqs(tmp_path, "sample1")
    manifest_path = tmp_path / "manifest.tsv"
    manifest_path.write_text(
        "sample_id\tr1\tr2\n"
        f"sample1\t{r1.name}\t{r2.name}\n"
        f"sample1\t{r1.name}\t{r2.name}\n"
    )

    with pytest.raises(ValueError):
        parse_sample_manifest(manifest_path)

def test_missing_fastq_file_raises(tmp_path: Path) -> None:
    manifest_path = tmp_path / "manifest.tsv"
    manifest_path.write_text("sample_id\tr1\tr2\n" "sample1\tdoes_not_exist_R1.fastq.gz\tdoes_not_exist_R2.fastq.gz\n")

    with pytest.raises(FileNotFoundError):
        parse_sample_manifest(manifest_path)

def test_missing_manifest_file_raises(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        parse_sample_manifest(tmp_path / "does_not_exist.tsv")

def test_empty_manifest_file_raises(tmp_path: Path) -> None:
    manifest_path = tmp_path / "manifest.tsv"
    manifest_path.write_text("")

    with pytest.raises(ValueError):
        parse_sample_manifest(manifest_path)
