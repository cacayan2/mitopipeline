"""
test_vertical_slice_blast_mitos2.py

Real, end-to-end integration test for #434: runs a candidate assembly
(the #403 gargammel fixture's own source reference, a real zebrafish
mitogenome) through BLAST identification, top-hit selection, reference
retrieval, and MITOS2 annotation -- the actual installed blastn and
runmitos.py binaries, not stubs -- and verifies provenance carries
correctly from one stage's output into the next.

Marked slow: requires the 05-blast and mitos2 conda environments and
the real Actinopterygii BLAST/MITOS2 reference data on this machine.
Skipped automatically if either isn't available, rather than failing
the whole suite on a machine that doesn't have them installed.
"""

# Imports
import importlib
import json
import shutil
import sys
from pathlib import Path
import pytest

run_blast_config = importlib.import_module("mitopipeline.exec.05a-blast_config.run_blast_config")
run_blast = importlib.import_module("mitopipeline.exec.05b-blast.run_blast")
run_blast_reference_retrieval = importlib.import_module("mitopipeline.exec.05c-blast_hits.run_blast_reference_retrieval")
run_mitos2_config = importlib.import_module("mitopipeline.exec.08a-mitos2_config.run_mitos2_config")
run_mitos2_annotation = importlib.import_module("mitopipeline.exec.08b-annotation.run_mitos2_annotation")

from mitopipeline.utils.blast.blast_hits import parse_blast_tsv, select_top_blast_hits, write_top_blast_hits

FIXTURE_FASTA = Path(__file__).parent.parent / "fixtures" / "gargammel_e2e" / "source_reference_NC_002333.2.fasta"
BLAST_DATABASE = Path("/home/eacac/mitogenome-batch-processor/v0.2.0/resources/blast/actinopterygii_mt_refs")
MITOS_REFDIR = Path("/home/eacac/Projects/mitogenome-batch-processor/resources/mitos")

pytestmark = pytest.mark.slow

def _require(condition: bool, message: str) -> None:
    if not condition:
        pytest.skip(message)

@pytest.fixture(autouse = True)
def _require_real_environment():
    _require(FIXTURE_FASTA.is_file(), f"Fixture not found: {FIXTURE_FASTA}")
    _require(shutil.which("blastn") is not None, "blastn not on PATH -- run inside the 05-blast conda env")
    _require(shutil.which("runmitos.py") is not None, "runmitos.py not on PATH -- run inside the mitos2 conda env")
    _require(BLAST_DATABASE.with_suffix(".nhr").is_file() or BLAST_DATABASE.with_suffix(".nin").is_file(),
              f"Real BLAST database not found: {BLAST_DATABASE}")
    _require(MITOS_REFDIR.is_dir(), f"Real MITOS2 reference data not found: {MITOS_REFDIR}")

def test_blast_selection_and_annotation_provenance_chains_correctly(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    sample_id = "VERTICAL_SLICE"

    # --- Stage 1: BLAST config + exec, against the real database ---
    monkeypatch.setattr(sys, "argv", [
        "run_blast_config.py",
        "--sample-id", sample_id,
        "--query-fasta", str(FIXTURE_FASTA),
        "--database", str(BLAST_DATABASE),
        "--mode", "local",
        "--blast-output-dir", str(tmp_path / "blast"),
        "--output-dir", str(tmp_path / "blast_config"),
        "--log-file", str(tmp_path / "logs" / "blast_config.log"),
        "--global-log-file", str(tmp_path / "logs" / "global.log"),
    ])
    assert run_blast_config.main() == 0
    blast_config_path = tmp_path / "blast_config" / f"{sample_id}.blast_config.json"

    monkeypatch.setattr(sys, "argv", [
        "run_blast.py",
        "--sample-id", sample_id,
        "--config-json", str(blast_config_path),
        "--working-dir", str(tmp_path),
        "--threads", "2",
        "--log-file", str(tmp_path / "logs" / "blast.log"),
        "--global-log-file", str(tmp_path / "logs" / "global.log"),
    ])
    assert run_blast.main() == 0
    blast_result = json.loads((tmp_path / f"{sample_id}.blast_result.json").read_text())
    assert blast_result["success"] is True
    assert blast_result["hit_count"] >= 1

    # --- Stage 2: rank/select top hit, write the TSV #401's retrieval reads ---
    hits = parse_blast_tsv(Path(blast_result["output_file"]))
    top_hits = select_top_blast_hits(hits, maximum_hits = 1)
    assert top_hits[0].sseqid == "ref|NC_002333.2|" or "NC_002333.2" in top_hits[0].sseqid
    top_hits_tsv = tmp_path / "top_hits.tsv"
    write_top_blast_hits(top_hits, top_hits_tsv, sample_id = sample_id)

    # --- Stage 3: reference retrieval, from the real database ---
    monkeypatch.setattr(sys, "argv", [
        "run_blast_reference_retrieval.py",
        "--sample-id", sample_id,
        "--top-hits-tsv", str(top_hits_tsv),
        "--database", str(BLAST_DATABASE),
        "--output-dir", str(tmp_path / "blast_reference"),
        "--working-dir", str(tmp_path),
        "--log-file", str(tmp_path / "logs" / "retrieval.log"),
        "--global-log-file", str(tmp_path / "logs" / "global.log"),
    ])
    assert run_blast_reference_retrieval.main() == 0
    retrieval_result = json.loads((tmp_path / "blast_reference" / f"{sample_id}.blast_reference_retrieval.json").read_text())
    assert retrieval_result["retrieved"] is True
    # Provenance: the retrieved accession must be the SAME one BLAST
    # actually selected as the top hit, not re-derived independently.
    assert retrieval_result["accession"] in top_hits[0].sseqid

    # --- Stage 4: MITOS2 config + annotation, on the SAME candidate assembly ---
    monkeypatch.setattr(sys, "argv", [
        "run_mitos2_config.py",
        "--sample-id", sample_id,
        "--input-fasta", str(FIXTURE_FASTA),
        "--genetic-code", "2",
        "--refdir", str(MITOS_REFDIR),
        "--refseqver", "refseq89m",
        "--noplots",
        "--annotation-output-dir", str(tmp_path / "annotation"),
        "--output-dir", str(tmp_path / "mitos2_config"),
        "--log-file", str(tmp_path / "logs" / "mitos2_config.log"),
        "--global-log-file", str(tmp_path / "logs" / "global.log"),
    ])
    assert run_mitos2_config.main() == 0
    mitos2_config_path = tmp_path / "mitos2_config" / f"{sample_id}.mitos2_config.json"

    monkeypatch.setattr(sys, "argv", [
        "run_mitos2_annotation.py",
        "--sample-id", sample_id,
        "--config-json", str(mitos2_config_path),
        "--working-dir", str(tmp_path),
        "--log-file", str(tmp_path / "logs" / "mitos2.log"),
        "--global-log-file", str(tmp_path / "logs" / "global.log"),
    ])
    assert run_mitos2_annotation.main() == 0
    annotation_result = json.loads((tmp_path / "annotation" / f"{sample_id}.mitos2_annotation_result.json").read_text())
    assert annotation_result["success"] is True
    # A real zebrafish mitogenome has 13 protein-coding genes, 22
    # tRNAs, and 2 rRNAs -- a real annotation run should find a
    # meaningful number of real genes, not zero and not a fabricated
    # fixed count.
    assert annotation_result["feature_count"] >= 10
    assert "cox1" in annotation_result["gene_names"] or "COX1" in annotation_result["gene_names"]

    # --- Provenance check across the whole slice ---
    # Every stage's own record carries the SAME sample_id -- nothing
    # downstream silently operated on a different sample.
    assert blast_result["sample_id"] == sample_id
    assert retrieval_result["sample_id"] == sample_id
    assert annotation_result["sample_id"] == sample_id
