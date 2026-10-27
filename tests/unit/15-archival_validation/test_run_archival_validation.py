"""
test_run_archival_validation.py

Fixture tests for the archival validation+packaging exec (#374/#375/
#450): combines real per-sample GenBank records, cross-checks every
table, and produces real checksums and a real submission checklist.
"""

# Imports
import importlib
import json
import sys
from pathlib import Path
import pytest

run_archival_validation = importlib.import_module("mitopipeline.exec.15-archival_validation.run_archival_validation")

def test_full_validation_and_packaging(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    r1, r2 = tmp_path / "r1.fq", tmp_path / "r2.fq"
    r1.write_text("x"); r2.write_text("y")
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text(f"sample_id\tr1\tr2\nsample1\t{r1}\t{r2}\n")

    sra_tsv = tmp_path / "sra_metadata.tsv"
    sra_tsv.write_text("sample_name\treview_status\tmissing_fields\nsample1\tREADY\t\n")
    biosample_tsv = tmp_path / "biosample_metadata.tsv"
    biosample_tsv.write_text("sample_name\torganism\treview_status\tmissing_fields\nsample1\tDanio rerio\tREADY\t\n")

    fasta = tmp_path / "data.fasta"
    fasta.write_text(">s\nACGT\n")
    genbank_records_dir = tmp_path / "genbank_records"
    genbank_records_dir.mkdir()
    (genbank_records_dir / "sample1.genbank_record.json").write_text(json.dumps({
        "sample_id": "sample1", "organism": "Danio rerio", "submission_folder": str(tmp_path / "gb"),
        "assembly_fasta": str(fasta), "annotation_gff": "", "sequence_length_bp": "4",
        "molecule_type": "DNA", "genome_location": "mitochondrion", "topology": "circular",
        "assembly_method": "GetOrganelle", "annotation_method": "MITOS2",
        "missing_files": "annotation_gff", "missing_metadata": "", "review_status": "MANUAL_REVIEW_REQUIRED",
    }))

    output_dir = tmp_path / "output"
    monkeypatch.setattr(sys, "argv", [
        "run_archival_validation.py", "--job-id", "job1", "--runtime-manifest", str(manifest),
        "--sra-metadata-tsv", str(sra_tsv), "--biosample-metadata-tsv", str(biosample_tsv),
        "--genbank-records-dir", str(genbank_records_dir), "--job-dir", str(tmp_path),
        "--output-dir", str(output_dir), "--log-file", str(tmp_path / "logs/s.log"),
        "--global-log-file", str(tmp_path / "logs/g.log"),
    ])

    assert run_archival_validation.main() == 0
    assert (output_dir / "genbank_metadata.tsv").is_file()
    assert (output_dir / "archival_validation.tsv").is_file()
    assert (output_dir / "md5_checksums.tsv").is_file()
    checklist = (output_dir / "submission_checklist.md").read_text()
    assert "sample1" in checklist
    assert "annotation_gff" in checklist or "GenBank needs review" in checklist

def test_handles_no_samples_gracefully(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text("sample_id\tr1\tr2\n")
    sra_tsv = tmp_path / "sra.tsv"
    sra_tsv.write_text("sample_name\n")
    biosample_tsv = tmp_path / "biosample.tsv"
    biosample_tsv.write_text("sample_name\n")
    genbank_records_dir = tmp_path / "genbank_records"
    output_dir = tmp_path / "output"

    monkeypatch.setattr(sys, "argv", [
        "run_archival_validation.py", "--job-id", "job1", "--runtime-manifest", str(manifest),
        "--sra-metadata-tsv", str(sra_tsv), "--biosample-metadata-tsv", str(biosample_tsv),
        "--genbank-records-dir", str(genbank_records_dir), "--job-dir", str(tmp_path),
        "--output-dir", str(output_dir), "--log-file", str(tmp_path / "logs/s.log"),
        "--global-log-file", str(tmp_path / "logs/g.log"),
    ])

    assert run_archival_validation.main() == 0
