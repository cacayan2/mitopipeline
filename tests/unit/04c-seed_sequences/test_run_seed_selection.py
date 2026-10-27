"""
test_run_seed_selection.py

Fixture tests for the seed-selection exec: wires select_seed() to a
CLI and records its outcome, including the fallback case.
"""

# Imports
import importlib
import json
import sys
from pathlib import Path
import pytest
from mitopipeline.utils.seed.models import SeedSelectionResult

run_seed_selection = importlib.import_module("mitopipeline.exec.04c-seed_sequences.run_seed_selection")

def _argv(tmp_path: Path, sample_id: str, scientific_name: str = "Danio rerio") -> list[str]:
    return [
        "run_seed_selection.py",
        "--sample-id", sample_id,
        "--scientific-name", scientific_name,
        "--database-dir", str(tmp_path / "mitofish_db"),
        "--output-dir", str(tmp_path / "seed_out"),
        "--log-file", str(tmp_path / "logs" / "seed.log"),
        "--global-log-file", str(tmp_path / "logs" / "global.log"),
    ]

def test_main_records_a_real_seed_selection(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    sample_id = "SAMPLE_01"
    seed_fasta = tmp_path / "seed_out" / f"{sample_id}.12S_rRNA.fasta"

    def _stub_select_seed(**kwargs):
        seed_fasta.parent.mkdir(parents = True, exist_ok = True)
        seed_fasta.write_text(">seed\nACGT\n")
        return SeedSelectionResult(
            sample_id = sample_id, scientific_name = "Danio rerio", gene = "12S_rRNA",
            accession = "NC_002333.2", source = "mitofish", fasta_path = seed_fasta,
            fallback_used = False, fallback_reason = None, match_level = "exact",
            checksum = "abc123", ncbi_taxon_id = 7954, sequence_length = 950,
        )

    monkeypatch.setattr(run_seed_selection, "select_seed", _stub_select_seed)
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id))

    assert run_seed_selection.main() == 0
    record = json.loads((tmp_path / "seed_out" / f"{sample_id}.seed_selection.json").read_text())
    assert record["fallback_used"] is False
    assert record["gene"] == "12S_rRNA"
    assert record["fasta_path"] == str(seed_fasta)

def test_main_records_a_fallback_result_as_a_valid_non_error_outcome(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    sample_id = "SAMPLE_02"

    def _stub_select_seed(**kwargs):
        return SeedSelectionResult(
            sample_id = sample_id, scientific_name = "Unknown species", gene = "none",
            accession = "GETORGANELLE_DEFAULT", source = "fallback_default_database",
            fasta_path = None, fallback_used = True, fallback_reason = "no_mitofish_candidate",
            match_level = None, checksum = None, ncbi_taxon_id = None, sequence_length = None,
        )

    monkeypatch.setattr(run_seed_selection, "select_seed", _stub_select_seed)
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, "Unknown species"))

    assert run_seed_selection.main() == 0   # fallback is valid, never an error
    record = json.loads((tmp_path / "seed_out" / f"{sample_id}.seed_selection.json").read_text())
    assert record["fallback_used"] is True
    assert record["fasta_path"] is None
