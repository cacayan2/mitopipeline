"""
test_run_data_aggregation.py

Fixture tests for the data-aggregation exec (Emil's own direct ask,
2026-10-28): every sample's own real assembly/annotation/phylogeny
outputs are consolidated into aggregated_data/, and a sample with no
assembly is skipped gracefully, not a crash.
"""

# Imports
import importlib
import json
import sys
from pathlib import Path
import pytest

run_data_aggregation = importlib.import_module("mitopipeline.exec.12b-data_aggregation.run_data_aggregation")

def _write_json(path: Path, data: dict) -> Path:
    path.parent.mkdir(parents = True, exist_ok = True)
    path.write_text(json.dumps(data))
    return path

def test_aggregates_a_real_sample_with_full_outputs(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    job_dir = tmp_path / "job"
    output_dir = job_dir / "aggregated_data"

    assembly_fasta = job_dir / "06-contig_selection" / "sample1" / "primary" / "data.fasta"
    assembly_fasta.parent.mkdir(parents = True)
    assembly_fasta.write_text(">sample1\nACGT\n")
    _write_json(job_dir / "06-contig_selection" / "sample1" / "sample1.contig_selection.json", {
        "sample_id": "sample1", "selected": True, "fasta_path": str(assembly_fasta),
    })

    map_png = job_dir / "09-visualization" / "sample1" / "sample1.circular_map.png"
    map_png.parent.mkdir(parents = True)
    map_png.write_bytes(b"fake png")
    _write_json(job_dir / "09-visualization" / "sample1" / "sample1.circular_map_result.json", {
        "success": True, "output_png": str(map_png),
    })

    gff = job_dir / "08b-annotation" / "sample1" / "sample1.annotation.gff"
    gff.parent.mkdir(parents = True)
    gff.write_text("##gff-version 3\n")
    _write_json(job_dir / "08b-annotation" / "sample1" / "sample1.mitos2_annotation_result.json", {
        "success": True, "gff_path": str(gff),
    })

    newick = job_dir / "11-phylogeny" / "sample1" / "sample1.treefile.nwk"
    newick.parent.mkdir(parents = True)
    newick.write_text("(a,b,c);\n")
    tree_png = job_dir / "11-phylogeny" / "sample1" / "tree_figure.png"
    tree_png.write_bytes(b"fake tree png")

    monkeypatch.setattr(sys, "argv", [
        "run_data_aggregation.py", "--job-id", "job1", "--sample-id", "sample1",
        "--job-dir", str(job_dir), "--output-dir", str(output_dir),
        "--log-file", str(tmp_path / "logs" / "agg.log"),
        "--global-log-file", str(tmp_path / "logs" / "global.log"),
    ])

    assert run_data_aggregation.main() == 0
    assert (output_dir / "assemblies" / "sample1.assembly.fasta").is_file()
    assert (output_dir / "annotations" / "sample1.annotated_map.png").is_file()
    assert (output_dir / "annotations" / "gff" / "sample1.annotation.gff").is_file()
    assert (output_dir / "phylogenies" / "sample1.treefile.nwk").is_file()
    assert (output_dir / "phylogenies" / "sample1.tree_figure.png").is_file()

    result = json.loads((output_dir / "job1.data_aggregation_result.json").read_text())
    assert result["samples"][0]["sample_id"] == "sample1"
    assert result["samples"][0]["assembly_path"] is not None

def test_sample_with_no_assembly_is_skipped_not_a_crash(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    job_dir = tmp_path / "job"
    output_dir = job_dir / "aggregated_data"

    monkeypatch.setattr(sys, "argv", [
        "run_data_aggregation.py", "--job-id", "job1", "--sample-id", "sample_bare",
        "--job-dir", str(job_dir), "--output-dir", str(output_dir),
        "--log-file", str(tmp_path / "logs" / "agg.log"),
        "--global-log-file", str(tmp_path / "logs" / "global.log"),
    ])

    assert run_data_aggregation.main() == 0
    result = json.loads((output_dir / "job1.data_aggregation_result.json").read_text())
    assert result["samples"][0]["assembly_path"] is None
    assert not (output_dir / "assemblies").exists()

def test_aggregates_real_archival_materials(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    job_dir = tmp_path / "job"
    output_dir = job_dir / "aggregated_data"

    sra = job_dir / "14a-sra_archival" / "sra_metadata.tsv"
    sra.parent.mkdir(parents = True)
    sra.write_text("sample_name\treview_status\nsample1\tREADY\n")

    biosample = job_dir / "14b-biosample_archival" / "biosample_metadata.tsv"
    biosample.parent.mkdir(parents = True)
    biosample.write_text("sample_name\treview_status\nsample1\tREADY\n")

    genbank_submission = job_dir / "14c-mitogenome_archival" / "sample1_Danio_rerio"
    genbank_submission.mkdir(parents = True)
    (genbank_submission / "sample1.fasta").write_text(">sample1\nACGT\n")
    (genbank_submission / "sample1.gff").write_text("##gff-version 3\n")
    _write_json(job_dir / "14c-mitogenome_archival" / "sample1.genbank_record.json", {
        "sample_id": "sample1", "submission_folder": str(genbank_submission),
    })

    validation = job_dir / "15-archival_validation" / "archival_validation.tsv"
    validation.parent.mkdir(parents = True)
    validation.write_text("sample_id\tready_for_review\nsample1\tTrue\n")
    checklist = job_dir / "15-archival_validation" / "submission_checklist.md"
    checklist.write_text("# checklist\n")
    checksums = job_dir / "15-archival_validation" / "md5_checksums.tsv"
    checksums.write_text("file\tmd5\na.fasta\tabc123\n")
    genbank_metadata = job_dir / "15-archival_validation" / "genbank_metadata.tsv"
    genbank_metadata.write_text("sample_name\treview_status\nsample1\tREADY\n")

    monkeypatch.setattr(sys, "argv", [
        "run_data_aggregation.py", "--job-id", "job1", "--sample-id", "sample1",
        "--job-dir", str(job_dir), "--output-dir", str(output_dir),
        "--sra-metadata-tsv", str(sra), "--biosample-metadata-tsv", str(biosample),
        "--genbank-metadata-tsv", str(genbank_metadata), "--archival-validation-tsv", str(validation),
        "--submission-checklist-md", str(checklist), "--checksums-tsv", str(checksums),
        "--log-file", str(tmp_path / "logs" / "agg.log"),
        "--global-log-file", str(tmp_path / "logs" / "global.log"),
    ])

    assert run_data_aggregation.main() == 0
    assert (output_dir / "archival" / "sra_metadata.tsv").is_file()
    assert (output_dir / "archival" / "biosample_metadata.tsv").is_file()
    assert (output_dir / "archival" / "genbank_metadata.tsv").is_file()
    assert (output_dir / "archival" / "archival_validation.tsv").is_file()
    assert (output_dir / "archival" / "submission_checklist.md").is_file()
    assert (output_dir / "archival" / "md5_checksums.tsv").is_file()
    assert (output_dir / "archival" / "genbank" / "sample1" / "sample1.fasta").read_text() == ">sample1\nACGT\n"
    assert (output_dir / "archival" / "genbank" / "sample1" / "sample1.gff").is_file()

    result = json.loads((output_dir / "job1.data_aggregation_result.json").read_text())
    assert result["samples"][0]["genbank_submission_dir_path"] is not None
    assert result["archival"]["sra_metadata_path"] is not None
    assert result["archival"]["submission_checklist_path"] is not None

def test_multiple_samples_are_isolated(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    job_dir = tmp_path / "job"
    output_dir = job_dir / "aggregated_data"

    for sample_id, content in (("sample_a", "AAAA"), ("sample_b", "BBBB")):
        fasta = job_dir / "06-contig_selection" / sample_id / "primary" / "data.fasta"
        fasta.parent.mkdir(parents = True)
        fasta.write_text(f">{sample_id}\n{content}\n")
        _write_json(job_dir / "06-contig_selection" / sample_id / f"{sample_id}.contig_selection.json", {
            "sample_id": sample_id, "selected": True, "fasta_path": str(fasta),
        })

    monkeypatch.setattr(sys, "argv", [
        "run_data_aggregation.py", "--job-id", "job1",
        "--sample-id", "sample_a", "--sample-id", "sample_b",
        "--job-dir", str(job_dir), "--output-dir", str(output_dir),
        "--log-file", str(tmp_path / "logs" / "agg.log"),
        "--global-log-file", str(tmp_path / "logs" / "global.log"),
    ])

    assert run_data_aggregation.main() == 0
    assert (output_dir / "assemblies" / "sample_a.assembly.fasta").read_text() == ">sample_a\nAAAA\n"
    assert (output_dir / "assemblies" / "sample_b.assembly.fasta").read_text() == ">sample_b\nBBBB\n"
