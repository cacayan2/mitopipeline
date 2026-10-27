"""
test_run_blast_hit_selection.py

Fixture tests for the BLAST hit ranking/selection exec: wires
rank_blast_hits/select_top_blast_hits/write_top_blast_hits to a CLI.
"""

# Imports
import importlib
import sys
from pathlib import Path
import pytest
from mitopipeline.utils.blast.blast_hits import read_top_blast_hits

run_blast_hit_selection = importlib.import_module("mitopipeline.exec.05c-blast_hits.run_blast_hit_selection")

def _row(sseqid: str, evalue: str = "1e-50") -> str:
    return "\t".join(["q1", sseqid, "98.5", "500", "2", "0", "1", "501", "1", "501",
                      evalue, "200", "95", "9606", "Homo sapiens", "title"])

def _argv(tmp_path: Path, sample_id: str, blast_tsv: Path, **extra) -> list[str]:
    argv = ["run_blast_hit_selection.py", "--sample-id", sample_id, "--blast-tsv", str(blast_tsv),
            "--output-dir", str(tmp_path / "out"),
            "--log-file", str(tmp_path / "logs" / "sel.log"), "--global-log-file", str(tmp_path / "logs" / "global.log")]
    if "maximum_hits" in extra:
        argv += ["--maximum-hits", str(extra["maximum_hits"])]
    return argv

def test_selects_and_writes_top_hits(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    sample_id = "SAMPLE_01"
    blast_tsv = tmp_path / "blast.tsv"
    blast_tsv.write_text(_row("strong", "1e-80") + "\n" + _row("weak", "1e-5") + "\n")
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, blast_tsv))

    assert run_blast_hit_selection.main() == 0
    hits = read_top_blast_hits(tmp_path / "out" / f"{sample_id}.top_hits.tsv")
    assert [h.sseqid for h in hits] == ["strong", "weak"]

def test_zero_hits_writes_an_empty_valid_tsv(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    sample_id = "SAMPLE_02"
    blast_tsv = tmp_path / "blast.tsv"
    blast_tsv.write_text("")
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, blast_tsv))

    assert run_blast_hit_selection.main() == 0
    assert read_top_blast_hits(tmp_path / "out" / f"{sample_id}.top_hits.tsv") == []

def test_missing_blast_tsv_is_recorded_not_a_crash(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, "SAMPLE_03", tmp_path / "does_not_exist.tsv"))
    assert run_blast_hit_selection.main() == 1

def test_maximum_hits_caps_selection(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    sample_id = "SAMPLE_04"
    blast_tsv = tmp_path / "blast.tsv"
    blast_tsv.write_text("\n".join(_row(f"s{i}") for i in range(10)) + "\n")
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, blast_tsv, maximum_hits = 3))

    assert run_blast_hit_selection.main() == 0
    assert len(read_top_blast_hits(tmp_path / "out" / f"{sample_id}.top_hits.tsv")) == 3
