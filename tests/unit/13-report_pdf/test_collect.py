"""
test_collect.py

Unit tests for assembling report data (#358-361) from real result
JSON: sample metadata, fastp stats, merged statistics, biological
validation verdicts, visualization results, and resolved tool
configs -- including missing-data handling (#359).
"""

# Imports
import json
from pathlib import Path
import pytest
from mitopipeline.utils.reporting.collect import (
    collect_phylogeny_summary,
    collect_run_report_data,
    collect_sample_report_data,
    parse_fastp_stats,
    read_sample_metadata,
)

def _write_json(path: Path, data: dict) -> Path:
    path.parent.mkdir(parents = True, exist_ok = True)
    path.write_text(json.dumps(data))
    return path

def test_read_sample_metadata_finds_the_right_row(tmp_path: Path) -> None:
    manifest = tmp_path / "runtime_manifest.tsv"
    manifest.write_text("sample_id\tr1\tr2\tgenus\tspecies\nsample1\ta\tb\tDanio\trerio\nsample2\tc\td\tCyprinus\tcarpio\n")
    metadata = read_sample_metadata(manifest, "sample2")
    assert metadata["genus"] == "Cyprinus"
    assert metadata["species"] == "carpio"

def test_read_sample_metadata_missing_manifest_returns_empty(tmp_path: Path) -> None:
    assert read_sample_metadata(tmp_path / "missing.tsv", "sample1") == {}

def test_parse_fastp_stats_from_a_real_shaped_report(tmp_path: Path) -> None:
    fastp_json = _write_json(tmp_path / "fastp.json", {
        "summary": {
            "before_filtering": {"total_reads": 1000, "q20_rate": 0.98, "q30_rate": 0.95, "gc_content": 0.45},
            "after_filtering": {"total_reads": 900, "q20_rate": 0.99, "q30_rate": 0.97, "gc_content": 0.46},
        },
        "filtering_result": {"low_quality_reads": 50},
    })
    stats = parse_fastp_stats(fastp_json)
    assert stats["reads_in"] == 1000
    assert stats["reads_out"] == 900
    assert stats["read_retention_rate"] == pytest.approx(0.9)
    assert stats["low_quality_reads"] == 50

def test_parse_fastp_stats_none_when_missing(tmp_path: Path) -> None:
    assert parse_fastp_stats(tmp_path / "missing.json") is None
    assert parse_fastp_stats(None) is None

def test_collect_sample_report_data_full(tmp_path: Path) -> None:
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text("sample_id\tr1\tr2\tgenus\tspecies\nsample1\ta\tb\tDanio\trerio\n")

    statistics_json = _write_json(tmp_path / "statistics.json", {
        "sample_id": "sample1",
        "assembly": {"sample_id": "sample1", "status": "circular", "source_attempt": 0, "length_bp": 16596,
                    "contig_count": 1, "gc_content_percent": 45.0, "recruited_read_count": 5000, "runtime_seconds": 12.5},
        "blast": {"sample_id": "sample1", "top_accession": "NC_1", "top_scientific_name": "Danio rerio",
                 "top_identity_percent": 99.0, "top_query_coverage_percent": 98.0, "top_evalue": 0.0, "decision_status": "PASS"},
        "annotation": {"sample_id": "sample1", "gene_count": 13, "core_gene_count": 13, "missing_core_genes": []},
    })
    validation_json = _write_json(tmp_path / "validation.json", {
        "sample_id": "sample1", "selected_attempt": 0,
        "results": [{"attempt": 0, "overall_verdict": "PASS",
                    "checks": [{"name": "length_bounds", "verdict": "PASS", "reason": "within bounds"}]}],
    })
    viz_json = _write_json(tmp_path / "viz.json", {"success": True, "output_png": str(tmp_path / "map.png")})
    (tmp_path / "map.png").write_bytes(b"fake")
    go_config = _write_json(tmp_path / "go_config.json", {"organelle_type": "animal_mt", "seeded": True, "seed_file": "12S_rRNA.fasta"})
    blast_config = _write_json(tmp_path / "blast_config.json", {"database": "actinopterygii_mt_refs", "mode": "blastn"})
    mitos_config = _write_json(tmp_path / "mitos_config.json", {"genetic_code": 2, "refseqver": "refseq89m"})

    data = collect_sample_report_data(
        "sample1", manifest, statistics_json = statistics_json, biological_validation_json = validation_json,
        visualization_result_json = viz_json, getorganelle_config_json = go_config,
        blast_config_json = blast_config, mitos2_config_json = mitos_config,
    )

    assert data.sample_id == "sample1"
    assert data.metadata["genus"] == "Danio"
    assert data.assembly.status == "circular"
    assert data.blast.top_accession == "NC_1"
    assert data.annotation.core_gene_count == 13
    assert data.overall_verdict == "PASS"
    assert len(data.validation_checks) == 1
    assert data.visualization_png == Path(tmp_path / "map.png")
    assert data.assembly_config.seeded is True
    assert data.blast_config.database == "actinopterygii_mt_refs"
    assert data.annotation_config.refseqver == "refseq89m"
    assert data.rescued is False

def test_collect_sample_report_data_handles_missing_stages(tmp_path: Path) -> None:
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text("sample_id\tr1\tr2\nsample1\ta\tb\n")

    data = collect_sample_report_data("sample1", manifest)

    assert data.assembly is None
    assert data.blast is None
    assert data.annotation is None
    assert data.overall_verdict is None
    assert data.visualization_png is None
    assert data.rescued is False

def test_rescued_true_when_source_attempt_is_nonzero(tmp_path: Path) -> None:
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text("sample_id\tr1\tr2\nsample1\ta\tb\n")
    statistics_json = _write_json(tmp_path / "statistics.json", {
        "sample_id": "sample1",
        "assembly": {"sample_id": "sample1", "status": "fragmented", "source_attempt": 1, "length_bp": 16000,
                    "contig_count": 1, "gc_content_percent": 44.0},
        "blast": None, "annotation": None,
    })
    data = collect_sample_report_data("sample1", manifest, statistics_json = statistics_json)
    assert data.rescued is True

def test_collect_phylogeny_summary(tmp_path: Path) -> None:
    result_json = _write_json(tmp_path / "iqtree_result.json", {
        "success": True, "tip_count": 5, "newick_path": str(tmp_path / "tree.nwk"), "reason": "completed",
    })
    summary = collect_phylogeny_summary(result_json, model = "GTR+G")
    assert summary.success is True
    assert summary.tip_count == 5
    assert summary.model == "GTR+G"

def test_collect_phylogeny_summary_none_when_missing(tmp_path: Path) -> None:
    assert collect_phylogeny_summary(tmp_path / "missing.json") is None

def test_collect_phylogeny_summary_includes_the_tree_figure_when_it_exists(tmp_path: Path) -> None:
    result_json = _write_json(tmp_path / "iqtree_result.json", {
        "success": True, "tip_count": 5, "newick_path": str(tmp_path / "tree.nwk"), "reason": "completed",
    })
    figure_png = tmp_path / "tree_figure.png"
    figure_png.write_bytes(b"fake png bytes")

    summary = collect_phylogeny_summary(result_json, tree_figure_png = figure_png)
    assert summary.tree_figure_png == figure_png

def test_collect_phylogeny_summary_omits_a_tree_figure_that_was_never_rendered(tmp_path: Path) -> None:
    result_json = _write_json(tmp_path / "iqtree_result.json", {
        "success": True, "tip_count": 5, "newick_path": str(tmp_path / "tree.nwk"), "reason": "completed",
    })
    summary = collect_phylogeny_summary(result_json, tree_figure_png = tmp_path / "never_rendered.png")
    assert summary.tree_figure_png is None

def test_collect_run_report_data_joins_samples(tmp_path: Path) -> None:
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text("sample_id\tr1\tr2\nsample1\ta\tb\n")
    sample = collect_sample_report_data("sample1", manifest)
    run_data = collect_run_report_data("job1", "2026-10-27T00:00:00Z", [sample])
    assert run_data.job_id == "job1"
    assert len(run_data.samples) == 1
