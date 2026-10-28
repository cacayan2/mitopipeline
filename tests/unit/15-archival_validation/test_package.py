"""
test_package.py

Unit tests for archival packaging (#375/#450): real, reproducible MD5
checksums over real files, and a submission checklist that correctly
reflects every real sample's own readiness and blockers.
"""

# Imports
import hashlib
from pathlib import Path
from mitopipeline.utils.archival.package import build_checksum_manifest, render_submission_checklist, write_checksum_manifest, write_submission_checklist

def test_checksum_manifest_is_real_and_reproducible(tmp_path: Path) -> None:
    file_a = tmp_path / "a.fasta"
    file_a.write_text(">s\nACGT\n")
    expected_md5 = hashlib.md5(file_a.read_bytes()).hexdigest()

    records = build_checksum_manifest([file_a], base_directory = tmp_path)

    assert len(records) == 1
    assert records[0]["path"] == "a.fasta"
    assert records[0]["md5"] == expected_md5

def test_checksum_manifest_skips_missing_files(tmp_path: Path) -> None:
    records = build_checksum_manifest([tmp_path / "does_not_exist.fasta"], base_directory = tmp_path)
    assert records == []

def test_checksum_manifest_is_sorted_by_path(tmp_path: Path) -> None:
    file_b = tmp_path / "b.fasta"
    file_a = tmp_path / "a.fasta"
    file_b.write_text("b")
    file_a.write_text("a")

    records = build_checksum_manifest([file_b, file_a], base_directory = tmp_path)
    assert [r["path"] for r in records] == ["a.fasta", "b.fasta"]

def test_write_checksum_manifest_writes_a_real_tsv(tmp_path: Path) -> None:
    file_a = tmp_path / "a.fasta"
    file_a.write_text("a")
    records = build_checksum_manifest([file_a], base_directory = tmp_path)
    output_path = tmp_path / "md5_checksums.tsv"

    write_checksum_manifest(records, output_path)
    assert output_path.is_file()

def test_submission_checklist_reflects_real_readiness() -> None:
    records = [
        {"sample_id": "sample1", "organism": "Danio rerio", "raw_fastq_ready": "True",
         "assembly_ready": "True", "annotation_ready": "True", "sra_ready": "True",
         "biosample_ready": "True", "genbank_ready": "True", "ready_for_review": "True",
         "submission_blockers": ""},
        {"sample_id": "sample2", "organism": "", "raw_fastq_ready": "False",
         "assembly_ready": "False", "annotation_ready": "False", "sra_ready": "False",
         "biosample_ready": "False", "genbank_ready": "False", "ready_for_review": "False",
         "submission_blockers": "raw FASTQ files missing; no SRA record"},
    ]

    checklist = render_submission_checklist(records, "job1")

    assert "1 of 2 sample(s) ready" in checklist
    assert "sample1" in checklist and "sample2" in checklist
    assert "raw FASTQ files missing; no SRA record" in checklist
    assert "never submits to NCBI directly" in checklist

def test_submission_checklist_attaches_a_hint_per_missing_field() -> None:
    records = [
        {"sample_id": "sample1", "organism": "Danio rerio", "raw_fastq_ready": "True",
         "assembly_ready": "True", "annotation_ready": "True", "sra_ready": "False",
         "biosample_ready": "False", "genbank_ready": "True", "ready_for_review": "False",
         "missing_metadata": "instrument_model;collection_date",
         "submission_blockers": "SRA fields need review: instrument_model; BioSample fields need review: collection_date"},
    ]

    checklist = render_submission_checklist(records, "job1")

    assert "**How to complete:**" in checklist
    assert "`instrument_model`: " in checklist
    assert "Illumina NovaSeq 6000" in checklist
    assert "`collection_date`: " in checklist
    assert "YYYY-MM-DD" in checklist

def test_submission_checklist_omits_how_to_complete_when_nothing_missing() -> None:
    records = [
        {"sample_id": "sample1", "organism": "Danio rerio", "raw_fastq_ready": "True",
         "assembly_ready": "True", "annotation_ready": "True", "sra_ready": "True",
         "biosample_ready": "True", "genbank_ready": "True", "ready_for_review": "True",
         "missing_metadata": "", "submission_blockers": ""},
    ]
    checklist = render_submission_checklist(records, "job1")
    assert "**How to complete:**" not in checklist

def test_write_submission_checklist_writes_a_real_file(tmp_path: Path) -> None:
    records = [{"sample_id": "sample1", "organism": "Danio rerio", "raw_fastq_ready": "True",
               "assembly_ready": "True", "annotation_ready": "True", "sra_ready": "True",
               "biosample_ready": "True", "genbank_ready": "True", "ready_for_review": "True",
               "submission_blockers": ""}]
    output_path = tmp_path / "submission_checklist.md"

    write_submission_checklist(records, "job1", output_path)
    assert output_path.is_file()
