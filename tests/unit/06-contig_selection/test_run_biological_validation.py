"""
test_run_biological_validation.py

Fixture tests for the biological validation exec (#414): applies
identical checks to the primary and every rescue attempt, and wires
BLAST/annotation evidence only onto the attempt contig_selection
actually chose.
"""

# Imports
import importlib
import json
import sys
from pathlib import Path
import pytest

run_biological_validation = importlib.import_module("mitopipeline.exec.06-contig_selection.run_biological_validation")

def _write_attempt(attempt_dir: Path, sample_id: str, status: str, length_bp: int = 16000) -> None:
    attempt_dir.mkdir(parents = True, exist_ok = True)
    (attempt_dir / f"{sample_id}.assembly_status.json").write_text(json.dumps({
        "sample_id": sample_id, "status": status, "success": status != "missing",
    }))
    if status in ("circular", "fragmented"):
        (attempt_dir / "data.fasta").write_text(f">seq\n{'ACGT' * (length_bp // 4)}\n")

def _write_selection(path: Path, sample_id: str, selected: bool, source_attempt: int | None = None, reason: str = "") -> None:
    path.parent.mkdir(parents = True, exist_ok = True)
    path.write_text(json.dumps({
        "sample_id": sample_id, "selected": selected, "source_attempt": source_attempt, "reason": reason,
    }))

def _argv(tmp_path: Path, sample_id: str, attempt_dirs: list[Path], selection_json: Path, **extra) -> list[str]:
    argv = ["run_biological_validation.py", "--sample-id", sample_id, "--contig-selection-json", str(selection_json),
            "--min-length", str(extra.get("min_length", 14000)), "--max-length", str(extra.get("max_length", 22000)),
            "--output-dir", str(tmp_path / "validation"),
            "--log-file", str(tmp_path / "logs" / "validation.log"), "--global-log-file", str(tmp_path / "logs" / "global.log")]
    for d in attempt_dirs:
        argv += ["--attempt-dir", str(d)]
    if "blast_result_json" in extra:
        argv += ["--blast-result-json", str(extra["blast_result_json"])]
    if "annotation_result_json" in extra:
        argv += ["--annotation-result-json", str(extra["annotation_result_json"])]
    return argv

def test_validates_primary_only_when_no_rescue_attempts(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    sample_id = "SAMPLE_01"
    primary_dir = tmp_path / "assembly"
    _write_attempt(primary_dir, sample_id, "circular")
    selection_json = tmp_path / f"{sample_id}.contig_selection.json"
    _write_selection(selection_json, sample_id, selected = True, source_attempt = 0, reason = "circular, 16000bp")

    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, [primary_dir], selection_json))
    assert run_biological_validation.main() == 0

    result = json.loads((tmp_path / "validation" / f"{sample_id}.biological_validation.json").read_text())
    assert len(result["results"]) == 1
    assert result["results"][0]["source"] == "primary"
    assert result["results"][0]["overall_verdict"] == "PASS"

def test_validates_primary_and_rescue_identically(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    sample_id = "SAMPLE_02"
    primary_dir = tmp_path / "assembly"
    rescue_dir = tmp_path / "rescue_attempt_1"
    _write_attempt(primary_dir, sample_id, "fragmented")
    _write_attempt(rescue_dir, sample_id, "circular")
    selection_json = tmp_path / f"{sample_id}.contig_selection.json"
    _write_selection(selection_json, sample_id, selected = True, source_attempt = 1, reason = "circular beats fragmented")

    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, [primary_dir, rescue_dir], selection_json))
    assert run_biological_validation.main() == 0

    result = json.loads((tmp_path / "validation" / f"{sample_id}.biological_validation.json").read_text())
    assert len(result["results"]) == 2
    sources = {r["attempt"]: r["source"] for r in result["results"]}
    assert sources == {0: "primary", 1: "rescue"}
    assert result["selected_attempt"] == 1

def test_blast_and_annotation_evidence_only_applies_to_the_selected_attempt(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    sample_id = "SAMPLE_03"
    primary_dir = tmp_path / "assembly"
    rescue_dir = tmp_path / "rescue_attempt_1"
    _write_attempt(primary_dir, sample_id, "fragmented")
    _write_attempt(rescue_dir, sample_id, "circular")
    selection_json = tmp_path / f"{sample_id}.contig_selection.json"
    _write_selection(selection_json, sample_id, selected = True, source_attempt = 1, reason = "circular beats fragmented")

    blast_result_json = tmp_path / f"{sample_id}.blast_result.json"
    blast_result_json.write_text(json.dumps({"success": True, "hit_count": 1, "output_file": None}))

    # A real top-hits TSV, written through the real write_top_blast_hits()
    # function this exec script actually reads back -- exercises the
    # real file format, not an assumption about its columns.
    from mitopipeline.utils.blast.blast_hits import write_top_blast_hits
    from mitopipeline.utils.blast.models import BlastHit
    top_hits_tsv = tmp_path / "top_hits.tsv"
    write_top_blast_hits([BlastHit(qseqid = "q1", sseqid = "s1", pident = 99.0, length = 500, mismatch = 1,
                                   gapopen = 0, qstart = 1, qend = 500, sstart = 1, send = 500, evalue = 1e-50,
                                   bitscore = 200.0, qcovs = 95.0, staxids = "9606", sscinames = "Danio rerio",
                                   stitle = "title", rank = 1)], top_hits_tsv, sample_id = sample_id)

    annotation_result_json = tmp_path / f"{sample_id}.mitos2_annotation_result.json"
    annotation_result_json.write_text(json.dumps({"success": True, "gene_names": ["cox1"]}))

    argv = _argv(tmp_path, sample_id, [primary_dir, rescue_dir], selection_json,
                blast_result_json = blast_result_json, annotation_result_json = annotation_result_json)
    argv += ["--blast-top-hits-tsv", str(top_hits_tsv)]
    monkeypatch.setattr(sys, "argv", argv)
    assert run_biological_validation.main() == 0

    result = json.loads((tmp_path / "validation" / f"{sample_id}.biological_validation.json").read_text())
    by_attempt = {r["attempt"]: r for r in result["results"]}
    selected_checks = {c["name"] for c in by_attempt[1]["checks"]}
    unselected_checks = {c["name"] for c in by_attempt[0]["checks"]}
    assert "gene_complement" in selected_checks
    assert "blast_identity" in selected_checks
    assert "gene_complement" not in unselected_checks   # the primary never got MITOS2 run on it
    assert "blast_identity" not in unselected_checks

def test_missing_assembly_attempt_is_skipped_not_crashed(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    sample_id = "SAMPLE_04"
    primary_dir = tmp_path / "assembly"
    _write_attempt(primary_dir, sample_id, "missing")
    selection_json = tmp_path / f"{sample_id}.contig_selection.json"
    _write_selection(selection_json, sample_id, selected = False)

    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, [primary_dir], selection_json))
    assert run_biological_validation.main() == 1   # nothing usable to validate

    result = json.loads((tmp_path / "validation" / f"{sample_id}.biological_validation.json").read_text())
    assert result["results"] == []

def test_missing_contig_selection_json_is_recorded_not_a_crash(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    sample_id = "SAMPLE_05"
    primary_dir = tmp_path / "assembly"
    _write_attempt(primary_dir, sample_id, "circular")
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, [primary_dir], tmp_path / "does_not_exist.json"))
    assert run_biological_validation.main() == 1
