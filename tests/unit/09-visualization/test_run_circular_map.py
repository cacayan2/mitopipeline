"""
test_run_circular_map.py

Fixture tests for the circular genome map rendering exec (#439): a
successful render records the PNG and feature count; a malformed
GFF or missing FASTA is recorded as this sample's own failed result,
never raised uncaught, so no other sample's run is affected (#440).
"""

# Imports
import importlib
import json
import sys
from pathlib import Path
import pytest

run_circular_map = importlib.import_module("mitopipeline.exec.09-visualization.run_circular_map")

REAL_GFF = (
    "##gff-version 3\n"
    "test\tmitos\tgene\t1\t600\t.\t+\t.\tID=gene_cox2;Name=cox2\n"
    "test\tmitos\tgene\t650\t1200\t.\t-\t.\tID=gene_nad1;Name=nad1\n"
)

def _argv(tmp_path: Path, sample_id: str, gff_path: Path, fasta_path: Path, output_dir: Path) -> list[str]:
    return [
        "run_circular_map.py",
        "--sample-id", sample_id,
        "--gff-path", str(gff_path),
        "--fasta-path", str(fasta_path),
        "--output-dir", str(output_dir),
        "--log-file", str(tmp_path / "logs" / "circular_map.log"),
        "--global-log-file", str(tmp_path / "logs" / "global.log"),
    ]

def test_successful_render_records_png_and_feature_count(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    sample_id = "SAMPLE_01"
    gff_path = tmp_path / "result.gff"
    gff_path.write_text(REAL_GFF)
    fasta_path = tmp_path / "data.fasta"
    fasta_path.write_text(">test\n" + "ACGT" * 400 + "\n")
    output_dir = tmp_path / "visualization" / sample_id

    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, gff_path, fasta_path, output_dir))

    assert run_circular_map.main() == 0
    result = json.loads((output_dir / f"{sample_id}.circular_map_result.json").read_text())
    assert result["success"] is True
    assert result["feature_count"] == 2
    output_png = Path(result["output_png"])
    assert output_png.is_file()
    assert output_png.stat().st_size > 0

def test_missing_fasta_is_recorded_not_raised(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    sample_id = "SAMPLE_02"
    gff_path = tmp_path / "result.gff"
    gff_path.write_text(REAL_GFF)
    missing_fasta = tmp_path / "does_not_exist.fasta"
    output_dir = tmp_path / "visualization" / sample_id

    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, gff_path, missing_fasta, output_dir))

    assert run_circular_map.main() == 1
    result = json.loads((output_dir / f"{sample_id}.circular_map_result.json").read_text())
    assert result["success"] is False
    assert result["output_png"] is None
    assert result["feature_count"] is None

def test_malformed_gff_is_recorded_not_raised(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    sample_id = "SAMPLE_03"
    gff_path = tmp_path / "result.gff"
    gff_path.write_text("test\tmitos\tgene\tnotanumber\t600\t.\t+\t.\tID=gene_cox2;Name=cox2\n")
    fasta_path = tmp_path / "data.fasta"
    fasta_path.write_text(">test\n" + "ACGT" * 400 + "\n")
    output_dir = tmp_path / "visualization" / sample_id

    monkeypatch.setattr(sys, "argv", _argv(tmp_path, sample_id, gff_path, fasta_path, output_dir))

    assert run_circular_map.main() == 1
    result = json.loads((output_dir / f"{sample_id}.circular_map_result.json").read_text())
    assert result["success"] is False

def test_one_samples_failure_does_not_affect_another(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    failing_sample, passing_sample = "FAIL_SAMPLE", "PASS_SAMPLE"

    gff_path = tmp_path / "result.gff"
    gff_path.write_text(REAL_GFF)
    fasta_path = tmp_path / "data.fasta"
    fasta_path.write_text(">test\n" + "ACGT" * 400 + "\n")

    failing_output_dir = tmp_path / "visualization" / failing_sample
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, failing_sample, gff_path, tmp_path / "missing.fasta", failing_output_dir))
    assert run_circular_map.main() == 1

    passing_output_dir = tmp_path / "visualization" / passing_sample
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, passing_sample, gff_path, fasta_path, passing_output_dir))
    assert run_circular_map.main() == 0

    passing_result = json.loads((passing_output_dir / f"{passing_sample}.circular_map_result.json").read_text())
    assert passing_result["success"] is True
