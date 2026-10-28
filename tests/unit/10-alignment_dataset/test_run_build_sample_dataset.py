"""
test_run_build_sample_dataset.py

Fixture tests for the per-sample phylogeny-dataset exec (restoring
v0.1.0's per-sample design): runs build_sample_phylogeny_dataset()
through the exec script's own argument parsing and result recording.
A missing assembly or missing BLAST hits is a normal, non-error
outcome this script must still complete for.
"""

# Imports
import importlib
import json
import sys
from pathlib import Path
import pytest
from mitopipeline.utils.blast.blast_hits import write_top_blast_hits
from mitopipeline.utils.blast.models import BlastHit

run_build_sample_dataset = importlib.import_module("mitopipeline.exec.10-alignment_dataset.run_build_sample_dataset")

def _hit(**overrides) -> BlastHit:
    defaults = dict(
        qseqid = "sample1", sseqid = "NC_000001.1", pident = 98.5, length = 16000,
        mismatch = 10, gapopen = 0, qstart = 1, qend = 16000, sstart = 1, send = 16000,
        evalue = 0.0, bitscore = 30000.0, qcovs = 99.0, staxids = "7998",
        sscinames = "Danio rerio", stitle = "Danio rerio mitochondrion, complete genome", rank = 1,
    )
    defaults.update(overrides)
    return BlastHit(**defaults)

def _argv(sample_id: str, assembly_fasta: Path, top_hits_tsv: Path, output_fasta: Path,
          output_dir: Path, tmp_path: Path) -> list[str]:
    return [
        "run_build_sample_dataset.py",
        "--sample-id", sample_id,
        "--assembly-fasta", str(assembly_fasta),
        "--top-hits-tsv", str(top_hits_tsv),
        "--output-fasta", str(output_fasta),
        "--entrez-email", "test@example.com",
        "--output-dir", str(output_dir),
        "--log-file", str(tmp_path / "logs" / "dataset.log"),
        "--global-log-file", str(tmp_path / "logs" / "global.log"),
    ]

def test_missing_assembly_is_recorded_not_a_crash(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    sample_id = "sample1"
    top_hits_tsv = tmp_path / "sample1.top_hits.tsv"
    write_top_blast_hits([_hit()], top_hits_tsv, sample_id = sample_id)
    output_dir = tmp_path / "output"

    monkeypatch.setattr(sys, "argv", _argv(
        sample_id, tmp_path / "does_not_exist.fasta", top_hits_tsv, tmp_path / "out.fasta", output_dir, tmp_path,
    ))

    assert run_build_sample_dataset.main() == 0
    result = json.loads((output_dir / f"{sample_id}.phylogeny_dataset_result.json").read_text())
    assert result["success"] is False
    assert result["reason"] == "no selected assembly"

def test_missing_top_hits_is_recorded_not_a_crash(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    sample_id = "sample1"
    assembly_fasta = tmp_path / "sample1.fasta"
    assembly_fasta.write_text(">contig_1\nATGCGT\n")
    output_dir = tmp_path / "output"

    monkeypatch.setattr(sys, "argv", _argv(
        sample_id, assembly_fasta, tmp_path / "does_not_exist.tsv", tmp_path / "out.fasta", output_dir, tmp_path,
    ))

    assert run_build_sample_dataset.main() == 0
    result = json.loads((output_dir / f"{sample_id}.phylogeny_dataset_result.json").read_text())
    assert result["success"] is False
    assert result["reason"] == "no BLAST hits available"

def test_successful_dataset_build_is_recorded(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    sample_id = "sample1"
    assembly_fasta = tmp_path / "sample1.fasta"
    assembly_fasta.write_text(">contig_1\nATGCGT\n")
    top_hits_tsv = tmp_path / "sample1.top_hits.tsv"
    write_top_blast_hits([_hit()], top_hits_tsv, sample_id = sample_id)
    output_fasta = tmp_path / "sample1.phylogeny_dataset.fasta"
    output_dir = tmp_path / "output"

    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from mitopipeline.utils.phylogeny import sample_dataset

    def fake_retrieve(accessions, entrez_email, entrez_api_key = None, logger = None):
        return [SeqRecord(Seq("AAAA"), id = "NC_000001.1")]

    monkeypatch.setattr(sample_dataset, "retrieve_reference_sequences", fake_retrieve)
    monkeypatch.setattr(sys, "argv", _argv(sample_id, assembly_fasta, top_hits_tsv, output_fasta, output_dir, tmp_path))

    assert run_build_sample_dataset.main() == 0
    result = json.loads((output_dir / f"{sample_id}.phylogeny_dataset_result.json").read_text())
    assert result["success"] is True
    assert result["reference_count"] == 1
    assert Path(result["output_fasta"]).is_file()
