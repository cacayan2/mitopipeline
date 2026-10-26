"""
test_run_blast_reference_retrieval.py

Fixture tests for the BLAST reference-retrieval exec (#401): reads
the rank-1 hit's accession from a top-hit TSV and retrieves its FASTA.
"""

# Imports
import importlib
import json
import sys
from datetime import datetime
from pathlib import Path
import pytest
from mitopipeline.api.blastdbcmd import BlastDbCmdRunner
from mitopipeline.models.command_result.command_result import CommandResult

run_blast_reference_retrieval = importlib.import_module(
    "mitopipeline.exec.05c-blast_hits.run_blast_reference_retrieval")

def _fake_result(success: bool = True, return_code: int = 0) -> CommandResult:
    now = datetime.now()
    return CommandResult(
        command = [], return_code = return_code, stdout = "", stderr = "", runtime_seconds = 0.0,
        success = success, tool_name = "blastdbcmd", started_at = now, ended_at = now,
        peak_memory_mb = 0.0, user_cpu_seconds = 0.0, system_cpu_seconds = 0.0,
    )

def _write_top_hits_tsv(path: Path, rows: list[dict]) -> None:
    import csv
    path.parent.mkdir(parents = True, exist_ok = True)
    columns = ["sample_id", "rank", "qseqid", "sseqid", "pident", "length", "mismatch",
               "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore",
               "qcovs", "staxids", "sscinames", "stitle"]
    with path.open("w", encoding = "utf-8", newline = "") as handle:
        writer = csv.DictWriter(handle, fieldnames = columns, delimiter = "\t")
        writer.writeheader()
        for row in rows:
            full = {c: "" for c in columns}
            full.update(row)
            writer.writerow(full)

def _argv(tmp_path: Path, sample_id: str, top_hits_tsv: Path) -> list[str]:
    return [
        "run_blast_reference_retrieval.py",
        "--sample-id", sample_id,
        "--top-hits-tsv", str(top_hits_tsv),
        "--database", "actinopterygii_mt_refs",
        "--output-dir", str(tmp_path / "blast_reference"),
        "--working-dir", str(tmp_path),
        "--log-file", str(tmp_path / "logs" / "retrieval.log"),
        "--global-log-file", str(tmp_path / "logs" / "global.log"),
    ]

def test_retrieves_the_rank_1_accession(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    sample_id = "SAMPLE_01"
    tsv = tmp_path / "top_hits.tsv"
    _write_top_hits_tsv(tsv, [
        {"sample_id": sample_id, "rank": "2", "sseqid": "NC_999999.1"},
        {"sample_id": sample_id, "rank": "1", "sseqid": "NC_025767.1"},
    ])

    def _stub_run(self: BlastDbCmdRunner) -> CommandResult:
        self.output_dir.mkdir(parents = True, exist_ok = True)
        self.output_fasta.write_text(f">{self.accession}\nACGT\n")
        return _fake_result()

    monkeypatch.setattr(BlastDbCmdRunner, "run", _stub_run)
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, tsv))

    assert run_blast_reference_retrieval.main() == 0

    record = json.loads((tmp_path / "blast_reference" / f"{sample_id}.blast_reference_retrieval.json").read_text())
    assert record["retrieved"] is True
    assert record["accession"] == "NC_025767.1"   # rank 1, not rank 2
    assert record["database"] == "actinopterygii_mt_refs"

def test_missing_tsv_is_recorded_not_a_crash(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    sample_id = "SAMPLE_02"
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, tmp_path / "does_not_exist.tsv"))
    assert run_blast_reference_retrieval.main() == 1
    record = json.loads((tmp_path / "blast_reference" / f"{sample_id}.blast_reference_retrieval.json").read_text())
    assert record["retrieved"] is False

def test_no_rank_1_row_for_this_sample_is_recorded(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    sample_id = "SAMPLE_03"
    tsv = tmp_path / "top_hits.tsv"
    _write_top_hits_tsv(tsv, [{"sample_id": "OTHER_SAMPLE", "rank": "1", "sseqid": "NC_111111.1"}])
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, tsv))
    assert run_blast_reference_retrieval.main() == 1
    record = json.loads((tmp_path / "blast_reference" / f"{sample_id}.blast_reference_retrieval.json").read_text())
    assert record["retrieved"] is False
    assert record["accession"] is None

def test_unavailable_accession_is_recorded_not_a_crash(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    sample_id = "SAMPLE_04"
    tsv = tmp_path / "top_hits.tsv"
    _write_top_hits_tsv(tsv, [{"sample_id": sample_id, "rank": "1", "sseqid": "NC_NOT_IN_DB.1"}])

    def _failing_run(self):
        raise SystemExit(1)

    monkeypatch.setattr(BlastDbCmdRunner, "run", _failing_run)
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, tsv))

    assert run_blast_reference_retrieval.main() == 1
    record = json.loads((tmp_path / "blast_reference" / f"{sample_id}.blast_reference_retrieval.json").read_text())
    assert record["retrieved"] is False
    assert record["accession"] == "NC_NOT_IN_DB.1"
    assert record["fasta_path"] is None
