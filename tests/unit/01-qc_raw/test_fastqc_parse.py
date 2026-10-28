"""
test_fastqc_parse.py

Unit tests for parsing FastQC's own real zip output format
(fastqc_data.txt/summary.txt/Images/*.png) -- the real raw text
format confirmed directly against a genuine FastQC 0.12.1 zip run on
this project's own gargammel_e2e fixtures.
"""

# Imports
import zipfile
from pathlib import Path
import pytest
from mitopipeline.utils.fastqc.parse import extract_fastqc_image, parse_fastqc_zip

REAL_FASTQC_DATA = (
    "##FastQC\t0.13.0\n"
    ">>Basic Statistics\tpass\n"
    "#Measure\tValue\n"
    "Filename\tsample_R1.fastq.gz\n"
    "File type\tConventional base calls\n"
    "Encoding\tSanger / Illumina 1.9\n"
    "Total Sequences\t12507\n"
    "Total Bases\t1.2 Mbp\n"
    "Sequence length\t100\n"
    "Mean Length\t100\n"
    "Median Length\t100\n"
    "%GC\t40\n"
    ">>END_MODULE\n"
    ">>Per base sequence quality\tpass\n"
    "#Base\tMean\n"
    "1\t32.8\n"
    ">>END_MODULE\n"
)

REAL_FASTQC_SUMMARY = (
    "PASS\tBasic Statistics\tsample_R1.fastq.gz\n"
    "PASS\tPer base sequence quality\tsample_R1.fastq.gz\n"
    "WARN\tPer sequence GC content\tsample_R1.fastq.gz\n"
    "PASS\tOverrepresented sequences\tsample_R1.fastq.gz\n"
)

def _write_fastqc_zip(zip_path: Path, internal_dir: str = "sample_R1_fastqc",
                       data_text: str = REAL_FASTQC_DATA, summary_text: str = REAL_FASTQC_SUMMARY,
                       include_image: bool = True) -> Path:
    with zipfile.ZipFile(zip_path, "w") as archive:
        archive.writestr(f"{internal_dir}/fastqc_data.txt", data_text)
        archive.writestr(f"{internal_dir}/summary.txt", summary_text)
        if include_image:
            archive.writestr(f"{internal_dir}/Images/per_base_quality.png", b"\x89PNG\r\n\x1a\nfake-png-bytes")
    return zip_path

def test_parse_fastqc_zip_reads_basic_statistics(tmp_path: Path) -> None:
    zip_path = _write_fastqc_zip(tmp_path / "sample.R1_fastqc.zip")
    summary = parse_fastqc_zip(zip_path)
    assert summary is not None
    assert summary.total_sequences == 12507
    assert summary.sequence_length == "100"
    assert summary.percent_gc == pytest.approx(40.0)

def test_parse_fastqc_zip_reads_module_verdicts_and_overall() -> None:
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        zip_path = _write_fastqc_zip(Path(tmp) / "sample.R1_fastqc.zip")
        summary = parse_fastqc_zip(zip_path)
    assert summary.module_verdicts["Per sequence GC content"] == "WARN"
    assert summary.module_verdicts["Basic Statistics"] == "PASS"
    # Worst verdict across every module wins -- one real WARN among PASSes.
    assert summary.overall_verdict == "WARN"

def test_parse_fastqc_zip_overall_verdict_is_fail_when_any_module_fails(tmp_path: Path) -> None:
    summary_with_fail = REAL_FASTQC_SUMMARY + "FAIL\tAdapter Content\tsample_R1.fastq.gz\n"
    zip_path = _write_fastqc_zip(tmp_path / "sample.R1_fastqc.zip", summary_text = summary_with_fail)
    summary = parse_fastqc_zip(zip_path)
    assert summary.overall_verdict == "FAIL"

def test_parse_fastqc_zip_handles_internal_dir_name_not_matching_external_filename(tmp_path: Path) -> None:
    """Regression: FastQC's own internal zip directory name is derived
    from the input FASTQ filename it was given, not the external zip
    filename a caller later renamed it to (confirmed directly against
    a real run: run_fastqc_raw.py renames gargammel_circular_R1_fastqc.zip
    to gargammel_circular.R1_fastqc.zip, but the zip's own internal
    gargammel_circular_R1_fastqc/ directory name is untouched)."""
    zip_path = _write_fastqc_zip(tmp_path / "gargammel_circular.R1_fastqc.zip", internal_dir = "gargammel_circular_R1_fastqc")
    summary = parse_fastqc_zip(zip_path)
    assert summary is not None
    assert summary.total_sequences == 12507

def test_parse_fastqc_zip_missing_file_returns_none(tmp_path: Path) -> None:
    assert parse_fastqc_zip(tmp_path / "does_not_exist.zip") is None

def test_parse_fastqc_zip_malformed_zip_returns_none(tmp_path: Path) -> None:
    bad_zip = tmp_path / "bad.zip"
    bad_zip.write_bytes(b"not a real zip file")
    assert parse_fastqc_zip(bad_zip) is None

def test_extract_fastqc_image_writes_a_real_file(tmp_path: Path) -> None:
    zip_path = _write_fastqc_zip(tmp_path / "sample.R1_fastqc.zip")
    output_path = tmp_path / "extracted.png"

    result = extract_fastqc_image(zip_path, "per_base_quality.png", output_path)

    assert result == output_path
    assert output_path.is_file()
    assert output_path.read_bytes().startswith(b"\x89PNG")

def test_extract_fastqc_image_missing_image_returns_none(tmp_path: Path) -> None:
    zip_path = _write_fastqc_zip(tmp_path / "sample.R1_fastqc.zip", include_image = False)
    assert extract_fastqc_image(zip_path, "per_base_quality.png", tmp_path / "out.png") is None
