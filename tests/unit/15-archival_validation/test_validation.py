"""
test_validation.py

Unit tests for archival cross-record validation (#374/#450): a
sample present in one real table but not another is a reportable
inconsistency, and the per-sample readiness record correctly
reflects every real upstream artifact and table review_status.
"""

# Imports
from pathlib import Path
from mitopipeline.utils.archival.validation import (
    build_archival_validation,
    validate_archival_readiness,
    validate_cross_record_identifiers,
)

def test_cross_record_identifiers_flags_a_sample_missing_from_one_table() -> None:
    sra_rows = {"sample1": {}, "sample2": {}}
    biosample_rows = {"sample1": {}}
    genbank_rows = {"sample1": {}, "sample2": {}}

    messages = validate_cross_record_identifiers(sra_rows, biosample_rows, genbank_rows)

    assert len(messages) == 1
    assert "sample2" in messages[0]
    assert "BioSample" in messages[0]

def test_cross_record_identifiers_clean_when_all_tables_agree() -> None:
    rows = {"sample1": {}}
    assert validate_cross_record_identifiers(rows, rows, rows) == []

def test_readiness_record_ready_when_every_source_is_ready() -> None:
    record = validate_archival_readiness(
        sample_id = "sample1", organism = "Danio rerio",
        raw_fastq_ready = True, assembly_ready = True, annotation_ready = True,
        sra_row = {"review_status": "READY", "missing_fields": ""},
        biosample_row = {"review_status": "READY", "missing_fields": ""},
        genbank_row = {"review_status": "READY", "missing_files": "", "missing_metadata": ""},
    )
    assert record["ready_for_review"] == "True"
    assert record["submission_blockers"] == ""

def test_readiness_record_lists_real_blockers_when_not_ready() -> None:
    record = validate_archival_readiness(
        sample_id = "sample1", organism = "Danio rerio",
        raw_fastq_ready = False, assembly_ready = True, annotation_ready = True,
        sra_row = None,
        biosample_row = {"review_status": "MANUAL_REVIEW_REQUIRED", "missing_fields": "collection_date"},
        genbank_row = {"review_status": "READY", "missing_files": "", "missing_metadata": ""},
    )
    assert record["ready_for_review"] == "False"
    assert "raw FASTQ files missing" in record["submission_blockers"]
    assert "no SRA record" in record["submission_blockers"]
    assert "collection_date" in record["submission_blockers"]

def test_readiness_record_dedupes_genbank_blocker_across_files_and_metadata() -> None:
    """Regression: GenBank's own missing_files and missing_metadata can
    legitimately share a field NAME for two different real reasons (a
    missing FILE vs a missing required METADATA field that reads empty
    precisely because that file was never copied) -- the combined
    submission_blockers string must not show that name twice."""
    record = validate_archival_readiness(
        sample_id = "sample1", organism = "Danio rerio",
        raw_fastq_ready = True, assembly_ready = False, annotation_ready = True,
        sra_row = {"review_status": "READY", "missing_fields": ""},
        biosample_row = {"review_status": "READY", "missing_fields": ""},
        genbank_row = {"review_status": "MANUAL_REVIEW_REQUIRED",
                       "missing_files": "annotation_gff;assembly_fasta",
                       "missing_metadata": "assembly_fasta"},
    )
    blocker = next(b for b in record["submission_blockers"].split("; ") if b.startswith("GenBank"))
    fields = blocker.removeprefix("GenBank needs review: ").split(";")
    assert fields.count("assembly_fasta") == 1
    assert set(fields) == {"annotation_gff", "assembly_fasta"}

def test_readiness_record_missing_metadata_covers_sra_biosample_and_genbank() -> None:
    """missing_metadata is the single source render_submission_checklist()
    reads to attach a hint per field -- it must carry every real
    missing field name across all three record types, not just SRA/
    BioSample (GenBank's own fields used to be visible only inside the
    free-text submission_blockers string)."""
    record = validate_archival_readiness(
        sample_id = "sample1", organism = "Danio rerio",
        raw_fastq_ready = True, assembly_ready = True, annotation_ready = True,
        sra_row = {"review_status": "MANUAL_REVIEW_REQUIRED", "missing_fields": "instrument_model"},
        biosample_row = {"review_status": "MANUAL_REVIEW_REQUIRED", "missing_fields": "collection_date;geo_loc_name"},
        genbank_row = {"review_status": "MANUAL_REVIEW_REQUIRED", "missing_files": "", "missing_metadata": "topology"},
    )
    fields = set(record["missing_metadata"].split(";"))
    assert fields == {"instrument_model", "collection_date", "geo_loc_name", "topology"}

def test_build_archival_validation_from_real_tsv_files(tmp_path: Path) -> None:
    r1, r2 = tmp_path / "r1.fq", tmp_path / "r2.fq"
    r1.write_text("read1")
    r2.write_text("read2")
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text(f"sample_id\tr1\tr2\nsample1\t{r1}\t{r2}\n")

    sra_tsv = tmp_path / "sra.tsv"
    sra_tsv.write_text("sample_name\treview_status\tmissing_fields\nsample1\tREADY\t\n")

    biosample_tsv = tmp_path / "biosample.tsv"
    biosample_tsv.write_text("sample_name\torganism\treview_status\tmissing_fields\nsample1\tDanio rerio\tREADY\t\n")

    genbank_tsv = tmp_path / "genbank.tsv"
    genbank_tsv.write_text(
        "sample_id\torganism\tassembly_fasta\tannotation_gff\treview_status\tmissing_files\tmissing_metadata\n"
        "sample1\tDanio rerio\t/real/path.fasta\t/real/path.gff\tREADY\t\t\n"
    )

    records = build_archival_validation(manifest, sra_tsv, biosample_tsv, genbank_tsv, tmp_path)

    assert len(records) == 1
    assert records[0]["ready_for_review"] == "True"
