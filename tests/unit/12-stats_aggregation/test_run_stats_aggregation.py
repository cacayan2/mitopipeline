"""
test_run_stats_aggregation.py

Fixture tests for the statistics aggregation exec (#356/#357): a
full real set of inputs merges into one record, and a sample missing
every optional input still records cleanly with null fields -- never
a crash, never a fabricated value.
"""

# Imports
import importlib
import json
import sys
from pathlib import Path
import pytest
from mitopipeline.utils.blast.blast_hits import write_top_blast_hits
from mitopipeline.utils.blast.models import BlastHit

run_stats_aggregation = importlib.import_module("mitopipeline.exec.12-stats_aggregation.run_stats_aggregation")

def _write_json(path: Path, data: dict) -> Path:
    path.parent.mkdir(parents = True, exist_ok = True)
    path.write_text(json.dumps(data))
    return path

def _argv(tmp_path: Path, sample_id: str, contig_selection_json: Path, **optional: Path) -> list[str]:
    argv = [
        "run_stats_aggregation.py",
        "--sample-id", sample_id,
        "--contig-selection-json", str(contig_selection_json),
        "--output-dir", str(tmp_path / "output"),
        "--log-file", str(tmp_path / "logs" / "stats.log"),
        "--global-log-file", str(tmp_path / "logs" / "global.log"),
    ]
    flag_names = {
        "recruited_reads_qc_json": "--recruited-reads-qc-json",
        "top_hits_tsv": "--top-hits-tsv",
        "biological_validation_json": "--biological-validation-json",
        "annotation_result_json": "--annotation-result-json",
    }
    for key, value in optional.items():
        argv.extend([flag_names[key], str(value)])
    return argv

def test_full_statistics_merge(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    sample_id = "sample1"
    fasta = tmp_path / "data.fasta"
    fasta.write_text(">s\nACGTACGTAC\n")
    selection_json = _write_json(tmp_path / "selection.json", {
        "sample_id": sample_id, "selected": True, "status": "circular",
        "fasta_path": str(fasta), "source_attempt": 0, "length_bp": 10,
    })
    top_hits_tsv = tmp_path / "top_hits.tsv"
    write_top_blast_hits([BlastHit(
        qseqid = sample_id, sseqid = "NC_000001.1", pident = 99.0, length = 16000, mismatch = 0,
        gapopen = 0, qstart = 1, qend = 16000, sstart = 1, send = 16000, evalue = 0.0,
        bitscore = 30000.0, qcovs = 99.0, staxids = "1", sscinames = "Example species",
        stitle = "Example species mitochondrion", rank = 1,
    )], top_hits_tsv, sample_id = sample_id)
    annotation_json = _write_json(tmp_path / "annotation.json", {
        "sample_id": sample_id, "success": True, "feature_count": 13,
        "gene_names": ["atp6", "atp8", "cob", "cox1", "cox2", "cox3",
                       "nad1", "nad2", "nad3", "nad4", "nad4l", "nad5", "nad6"],
    })

    monkeypatch.setattr(sys, "argv", _argv(
        tmp_path, sample_id, selection_json,
        top_hits_tsv = top_hits_tsv, annotation_result_json = annotation_json,
    ))

    assert run_stats_aggregation.main() == 0
    result = json.loads((tmp_path / "output" / f"{sample_id}.statistics.json").read_text())
    assert result["assembly"]["status"] == "circular"
    assert result["blast"]["top_accession"] == "NC_000001.1"
    assert result["annotation"]["core_gene_count"] == 13
    assert result["annotation"]["missing_core_genes"] == []

def test_missing_optional_inputs_records_null_fields_not_a_crash(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    sample_id = "sample_no_hits"
    selection_json = _write_json(tmp_path / "selection.json", {
        "sample_id": sample_id, "selected": False, "attempts_considered": 1,
    })

    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, selection_json))

    assert run_stats_aggregation.main() == 0
    result = json.loads((tmp_path / "output" / f"{sample_id}.statistics.json").read_text())
    assert result == {"sample_id": sample_id, "assembly": None, "blast": None, "annotation": None}
