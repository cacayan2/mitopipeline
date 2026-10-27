"""
test_extract.py

Unit tests for per-sample statistics extraction (#353-356): assembly
statistics from contig_selection's own record, BLAST statistics from
the ranked top-hit TSV plus the validation gate's own verdict,
annotation statistics against the core animal mitochondrial gene set,
and merging -- including missing-stage and sample-isolation handling
(#357).
"""

# Imports
import json
from pathlib import Path
import pytest
from mitopipeline.utils.blast.blast_hits import write_top_blast_hits
from mitopipeline.utils.blast.models import BlastHit
from mitopipeline.utils.statistics.extract import (
    extract_annotation_statistics,
    extract_assembly_statistics,
    extract_blast_statistics,
    merge_sample_statistics,
)

def _write_json(path: Path, data: dict) -> Path:
    path.parent.mkdir(parents = True, exist_ok = True)
    path.write_text(json.dumps(data))
    return path

def _write_fasta(path: Path, records: list[tuple[str, str]]) -> Path:
    path.parent.mkdir(parents = True, exist_ok = True)
    path.write_text("".join(f">{name}\n{seq}\n" for name, seq in records))
    return path

# --- Assembly statistics (#353) ---

def test_extract_assembly_statistics_from_a_real_selection(tmp_path: Path) -> None:
    fasta = _write_fasta(tmp_path / "data.fasta", [("gargammel_circular", "ACGTACGTAC")])  # 10bp, 50% GC
    selection_json = _write_json(tmp_path / "selection.json", {
        "sample_id": "sample1", "selected": True, "status": "circular",
        "fasta_path": str(fasta), "source_attempt": 0, "length_bp": 10,
    })
    qc_json = _write_json(tmp_path / "qc.json", {"r1_read_count": 100, "r2_read_count": 100})

    stats = extract_assembly_statistics("sample1", selection_json, recruited_reads_qc_json = qc_json, runtime_seconds = 42.0)

    assert stats.sample_id == "sample1"
    assert stats.status == "circular"
    assert stats.source_attempt == 0
    assert stats.length_bp == 10
    assert stats.contig_count == 1
    assert stats.gc_content_percent == pytest.approx(50.0)
    assert stats.recruited_read_count == 200
    assert stats.runtime_seconds == 42.0

def test_extract_assembly_statistics_without_optional_fields(tmp_path: Path) -> None:
    fasta = _write_fasta(tmp_path / "data.fasta", [("s", "AAAA")])
    selection_json = _write_json(tmp_path / "selection.json", {
        "sample_id": "sample1", "selected": True, "status": "fragmented",
        "fasta_path": str(fasta), "source_attempt": 1, "length_bp": 4,
    })

    stats = extract_assembly_statistics("sample1", selection_json)
    assert stats.recruited_read_count is None
    assert stats.runtime_seconds is None

def test_extract_assembly_statistics_none_when_nothing_selected(tmp_path: Path) -> None:
    selection_json = _write_json(tmp_path / "selection.json", {
        "sample_id": "sample1", "selected": False, "attempts_considered": 3,
    })
    assert extract_assembly_statistics("sample1", selection_json) is None

def test_extract_assembly_statistics_none_when_selection_file_missing(tmp_path: Path) -> None:
    assert extract_assembly_statistics("sample1", tmp_path / "does_not_exist.json") is None

# --- BLAST statistics (#354) ---

def _real_blast_hit(**overrides) -> BlastHit:
    defaults = dict(
        qseqid = "sample1", sseqid = "NC_002333.2", pident = 98.5, length = 16000,
        mismatch = 10, gapopen = 0, qstart = 1, qend = 16000, sstart = 1, send = 16000,
        evalue = 0.0, bitscore = 30000.0, qcovs = 99.0, staxids = "7998",
        sscinames = "Danio rerio", stitle = "Danio rerio mitochondrion, complete genome", rank = 1,
    )
    defaults.update(overrides)
    return BlastHit(**defaults)

def test_extract_blast_statistics_from_the_real_top_hit(tmp_path: Path) -> None:
    top_hits_tsv = tmp_path / "top_hits.tsv"
    write_top_blast_hits([_real_blast_hit()], top_hits_tsv, sample_id = "sample1")
    validation_json = _write_json(tmp_path / "validation.json", {
        "sample_id": "sample1", "selected_attempt": 0,
        "results": [{"attempt": 0, "checks": [{"name": "blast_identity", "verdict": "PASS", "reason": "strong match"}]}],
    })

    stats = extract_blast_statistics("sample1", top_hits_tsv, biological_validation_json = validation_json)

    assert stats.top_accession == "NC_002333.2"
    assert stats.top_scientific_name == "Danio rerio"
    assert stats.top_identity_percent == pytest.approx(98.5)
    assert stats.top_query_coverage_percent == pytest.approx(99.0)
    assert stats.decision_status == "PASS"

def test_extract_blast_statistics_without_validation_json(tmp_path: Path) -> None:
    top_hits_tsv = tmp_path / "top_hits.tsv"
    write_top_blast_hits([_real_blast_hit()], top_hits_tsv, sample_id = "sample1")
    stats = extract_blast_statistics("sample1", top_hits_tsv)
    assert stats.decision_status is None

def test_extract_blast_statistics_none_when_no_hits(tmp_path: Path) -> None:
    top_hits_tsv = tmp_path / "top_hits.tsv"
    write_top_blast_hits([], top_hits_tsv, sample_id = "sample1")
    assert extract_blast_statistics("sample1", top_hits_tsv) is None

def test_extract_blast_statistics_none_when_tsv_missing(tmp_path: Path) -> None:
    assert extract_blast_statistics("sample1", tmp_path / "does_not_exist.tsv") is None

def test_extract_blast_statistics_uses_circular_coverage_correction(tmp_path: Path) -> None:
    top_hits_tsv = tmp_path / "top_hits.tsv"
    write_top_blast_hits([_real_blast_hit(qcovs = 53.0)], top_hits_tsv, sample_id = "sample1")
    circular_coverage_json = _write_json(tmp_path / "circular_coverage.json", {
        "success": True, "corrected_query_coverage": 100.0,
    })

    stats = extract_blast_statistics("sample1", top_hits_tsv, circular_coverage_json = circular_coverage_json)

    assert stats.top_query_coverage_percent == pytest.approx(100.0)

def test_extract_blast_statistics_ignores_a_failed_circular_coverage_correction(tmp_path: Path) -> None:
    top_hits_tsv = tmp_path / "top_hits.tsv"
    write_top_blast_hits([_real_blast_hit(qcovs = 53.0)], top_hits_tsv, sample_id = "sample1")
    circular_coverage_json = _write_json(tmp_path / "circular_coverage.json", {"success": False})

    stats = extract_blast_statistics("sample1", top_hits_tsv, circular_coverage_json = circular_coverage_json)

    assert stats.top_query_coverage_percent == pytest.approx(53.0)

# --- Annotation statistics (#355) ---

def test_extract_annotation_statistics_counts_core_genes(tmp_path: Path) -> None:
    result_json = _write_json(tmp_path / "annotation.json", {
        "sample_id": "sample1", "success": True, "feature_count": 15,
        "gene_names": ["cox1", "cox2", "cox3", "cob", "atp6", "atp8",
                       "nad1", "nad2", "nad3", "nad4", "nad5", "nad6", "trnL", "trnM"],
    })

    stats = extract_annotation_statistics("sample1", result_json)

    assert stats.gene_count == 15
    assert stats.core_gene_count == 12
    assert stats.missing_core_genes == ("nad4l",)

def test_extract_annotation_statistics_none_when_run_failed(tmp_path: Path) -> None:
    result_json = _write_json(tmp_path / "annotation.json", {"sample_id": "sample1", "success": False})
    assert extract_annotation_statistics("sample1", result_json) is None

def test_extract_annotation_statistics_none_when_missing(tmp_path: Path) -> None:
    assert extract_annotation_statistics("sample1", tmp_path / "does_not_exist.json") is None

# --- Merge + sample isolation (#356/#357) ---

def test_merge_sample_statistics_keeps_missing_stages_as_none() -> None:
    merged = merge_sample_statistics("sample1", assembly = None, blast = None, annotation = None)
    assert merged.sample_id == "sample1"
    assert merged.assembly is None and merged.blast is None and merged.annotation is None
    assert merged.to_dict() == {"sample_id": "sample1", "assembly": None, "blast": None, "annotation": None}

def test_one_samples_missing_stage_does_not_affect_another(tmp_path: Path) -> None:
    fasta = _write_fasta(tmp_path / "good" / "data.fasta", [("s", "ACGT")])
    good_selection = _write_json(tmp_path / "good" / "selection.json", {
        "sample_id": "good_sample", "selected": True, "status": "circular",
        "fasta_path": str(fasta), "source_attempt": 0, "length_bp": 4,
    })
    bad_selection_path = tmp_path / "bad" / "does_not_exist.json"

    good_stats = extract_assembly_statistics("good_sample", good_selection)
    bad_stats = extract_assembly_statistics("bad_sample", bad_selection_path)

    assert good_stats is not None
    assert bad_stats is None
