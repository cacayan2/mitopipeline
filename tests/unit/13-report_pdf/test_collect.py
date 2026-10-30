"""
test_collect.py

Unit tests for assembling report data (#358-361) from real result
JSON: sample metadata, fastp stats, merged statistics, biological
validation verdicts, visualization results, and resolved tool
configs -- including missing-data handling (#359).
"""

# Imports
import json
import zipfile
from pathlib import Path
import pytest
from mitopipeline.utils.reporting.collect import (
    collect_job_metadata,
    collect_phylogeny_summary,
    collect_read_qc,
    collect_rescue_summary,
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

def test_parse_fastp_stats_extracts_real_settings_from_the_command_string(tmp_path: Path) -> None:
    """fastp's own JSON report records the exact real command it was
    invoked with (confirmed directly) -- settings are parsed from
    that, never re-derived from a separately-maintained config copy
    that could drift from what actually ran."""
    fastp_json = _write_json(tmp_path / "fastp.json", {
        "summary": {
            "before_filtering": {"total_reads": 1000, "q20_rate": 0.98, "q30_rate": 0.95, "gc_content": 0.45},
            "after_filtering": {"total_reads": 900, "q20_rate": 0.99, "q30_rate": 0.97, "gc_content": 0.46},
        },
        "filtering_result": {"low_quality_reads": 50},
        "command": (
            "fastp --in1 a.fastq.gz --in2 b.fastq.gz --qualified_quality_phred 20 "
            "--length_required 50 --cut_window_size 4 --cut_mean_quality 20 "
            "--n_base_limit 5 --unqualified_percent_limit 40 --detect_adapter_for_pe"
        ),
    })
    stats = parse_fastp_stats(fastp_json)
    assert stats["qualified_quality_phred"] == 20
    assert stats["length_required"] == 50
    assert stats["cut_window_size"] == 4
    assert stats["cut_mean_quality"] == 20
    assert stats["n_base_limit"] == 5
    assert stats["unqualified_percent_limit"] == 40
    assert stats["adapter_detection"] is True

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

def _write_fastqc_zip(zip_path: Path, total_sequences: int = 1000, percent_gc: float = 40.0,
                      sequence_length: str = "100", verdicts: str = "PASS\tBasic Statistics\tx\n") -> Path:
    zip_path.parent.mkdir(parents = True, exist_ok = True)
    data_text = (
        "##FastQC\t0.13.0\n"
        ">>Basic Statistics\tpass\n"
        "#Measure\tValue\n"
        f"Total Sequences\t{total_sequences}\n"
        f"Sequence length\t{sequence_length}\n"
        f"%GC\t{percent_gc}\n"
        ">>END_MODULE\n"
    )
    with zipfile.ZipFile(zip_path, "w") as archive:
        archive.writestr("x_fastqc/fastqc_data.txt", data_text)
        archive.writestr("x_fastqc/summary.txt", verdicts)
        archive.writestr("x_fastqc/Images/per_base_quality.png", b"\x89PNG\r\n\x1a\nfake")
    return zip_path

def test_collect_read_qc_combines_r1_and_r2(tmp_path: Path) -> None:
    r1 = _write_fastqc_zip(tmp_path / "r1.zip", total_sequences = 1000, percent_gc = 40.0)
    r2 = _write_fastqc_zip(tmp_path / "r2.zip", total_sequences = 900, percent_gc = 42.0)

    summary = collect_read_qc(r1, r2, image_output_path = tmp_path / "out.png")

    assert summary is not None
    assert summary.total_sequences == 1900
    assert summary.percent_gc == pytest.approx(41.0)
    assert summary.sequence_length == "100"
    assert summary.overall_verdict == "PASS"
    assert summary.representative_image_png == tmp_path / "out.png"
    assert (tmp_path / "out.png").is_file()

def test_collect_read_qc_overall_verdict_is_worst_of_r1_and_r2(tmp_path: Path) -> None:
    r1 = _write_fastqc_zip(tmp_path / "r1.zip", verdicts = "PASS\tBasic Statistics\tx\n")
    r2 = _write_fastqc_zip(tmp_path / "r2.zip", verdicts = "WARN\tPer sequence GC content\tx\n")

    summary = collect_read_qc(r1, r2)

    assert summary.overall_verdict == "WARN"
    assert summary.failed_modules == ()
    assert summary.warned_modules == ("Per sequence GC content",)

def test_collect_read_qc_reports_failed_and_warned_module_names(tmp_path: Path) -> None:
    r1 = _write_fastqc_zip(
        tmp_path / "r1.zip",
        verdicts = "PASS\tBasic Statistics\tx\nFAIL\tPer base sequence quality\tx\nWARN\tAdapter Content\tx\n",
    )
    r2 = _write_fastqc_zip(
        tmp_path / "r2.zip",
        verdicts = "WARN\tPer base sequence quality\tx\nFAIL\tSequence Duplication Levels\tx\n",
    )

    summary = collect_read_qc(r1, r2)

    assert summary.overall_verdict == "FAIL"
    # "Per base sequence quality" FAILs on R1 and only WARNs on R2 --
    # the worse of the two (FAIL) must win, not whichever read is
    # processed last.
    assert summary.failed_modules == ("Per base sequence quality", "Sequence Duplication Levels")
    assert summary.warned_modules == ("Adapter Content",)

def test_collect_read_qc_none_when_both_zips_missing(tmp_path: Path) -> None:
    assert collect_read_qc(tmp_path / "missing_r1.zip", tmp_path / "missing_r2.zip") is None

def test_collect_read_qc_handles_one_missing_zip(tmp_path: Path) -> None:
    r1 = _write_fastqc_zip(tmp_path / "r1.zip", total_sequences = 500)
    summary = collect_read_qc(r1, tmp_path / "missing_r2.zip")
    assert summary is not None
    assert summary.total_sequences == 500

def test_collect_sample_report_data_enriches_assembly_config_with_seed_selection(tmp_path: Path) -> None:
    """getorganelle_config.json only ever records the seed FILE path --
    the real gene/accession/source comes from the separate real
    seed_selection.json (run_seed_selection.py's own output)."""
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text("sample_id\tr1\tr2\nsample1\ta\tb\n")
    go_config = _write_json(tmp_path / "go_config.json", {"organelle_type": "animal_mt", "seeded": True, "seed_file": "sample1.seed.fasta"})
    seed_selection = _write_json(tmp_path / "seed_selection.json", {
        "sample_id": "sample1", "scientific_name": "Danio rerio", "gene": "COXI",
        "accession": "NC_002333.2", "source": "mitofish", "fallback_used": False,
    })

    data = collect_sample_report_data(
        "sample1", manifest, getorganelle_config_json = go_config, seed_selection_json = seed_selection,
    )

    assert data.assembly_config.seed_gene == "COXI"
    assert data.assembly_config.seed_accession == "NC_002333.2"
    assert data.assembly_config.seed_source == "mitofish"
    assert data.assembly_config.seed_fallback_reason is None

def test_collect_sample_report_data_records_seed_fallback_reason_when_unseeded(tmp_path: Path) -> None:
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text("sample_id\tr1\tr2\nsample1\ta\tb\n")
    go_config = _write_json(tmp_path / "go_config.json", {"organelle_type": "animal_mt", "seeded": False, "seed_file": None})
    seed_selection = _write_json(tmp_path / "seed_selection.json", {
        "sample_id": "sample1", "scientific_name": "Unknownus fishicus", "gene": "COXI",
        "accession": "GETORGANELLE_DEFAULT", "source": "fallback_default_database",
        "fallback_used": True, "fallback_reason": "no_fetchable_mitofish_match",
    })

    data = collect_sample_report_data(
        "sample1", manifest, getorganelle_config_json = go_config, seed_selection_json = seed_selection,
    )

    assert data.assembly_config.seed_accession is None
    assert data.assembly_config.seed_fallback_reason == "no_fetchable_mitofish_match"

def test_collect_rescue_summary_reads_every_real_attempt(tmp_path: Path) -> None:
    """Regression (Emil's direct ask, 2026-10-28): the methods section
    needs every real rescue attempt, not only the one that succeeded
    -- a real failed/skipped attempt (e.g. gargammel_fragmented's own
    attempt 1 producing "missing" instead of a better result, and
    attempt 2 never running because of it) must still be reported."""
    rescue_dir = tmp_path / "07-assembly_rescue" / "sample1"
    attempt_1 = rescue_dir / "attempt_1"
    attempt_1.mkdir(parents = True)
    _write_json(attempt_1 / "sample1.rescue_mapping.json", {
        "sample_id": "sample1", "attempted": True,
        "reason": "status 'fragmented' triggers rescue; attempt 1 of 2",
    })
    _write_json(attempt_1 / "reassembly" / "sample1.assembly_status.json", {
        "sample_id": "sample1", "status": "missing", "success": False, "profile_name": "reduced_coverage",
    })
    attempt_2 = rescue_dir / "attempt_2"
    attempt_2.mkdir(parents = True)
    _write_json(attempt_2 / "sample1.rescue_mapping.json", {
        "sample_id": "sample1", "attempted": False,
        "reason": "status 'missing' is not in trigger_statuses ['fragmented', 'partial']",
    })

    attempts = collect_rescue_summary(rescue_dir)

    assert len(attempts) == 2
    assert attempts[0].attempt == 1
    assert attempts[0].attempted is True
    assert attempts[0].profile_name == "reduced_coverage"
    assert attempts[0].status == "missing"
    assert attempts[1].attempt == 2
    assert attempts[1].attempted is False
    assert attempts[1].profile_name is None

def test_collect_rescue_summary_empty_when_no_rescue_dir(tmp_path: Path) -> None:
    assert collect_rescue_summary(None) == ()
    assert collect_rescue_summary(tmp_path / "does_not_exist") == ()

def test_collect_run_report_data_joins_samples(tmp_path: Path) -> None:
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text("sample_id\tr1\tr2\nsample1\ta\tb\n")
    sample = collect_sample_report_data("sample1", manifest)
    run_data = collect_run_report_data("job1", "2026-10-27T00:00:00Z", [sample])
    assert run_data.job_id == "job1"
    assert len(run_data.samples) == 1

def test_collect_job_metadata_merges_current_and_historical(tmp_path: Path) -> None:
    current = _write_json(tmp_path / "current_metadata.json", {
        "status": "completed", "start_time": "2026-10-28T09:28:39", "runtime_seconds": 217.6, "sample_count": 3,
    })
    historical = tmp_path / "historical_metadata.json"
    historical.write_text(json.dumps([
        {"status": "created", "runtime_seconds": 0.0},
        {"status": "completed", "runtime_seconds": 217.6},
        {"status": "created", "runtime_seconds": 0.0},
        {"status": "failed", "runtime_seconds": 5.2},
    ]))

    meta = collect_job_metadata(current, historical, "job1")

    assert meta.status == "completed"
    assert meta.runtime_seconds == 217.6
    assert meta.sample_count == 3
    assert meta.attempt_count == 2
    assert round(meta.cumulative_runtime_seconds, 1) == 222.8
    assert meta.incomplete_attempts == 0

def test_collect_job_metadata_none_when_nothing_available(tmp_path: Path) -> None:
    assert collect_job_metadata(None, None, "job1") is None
    assert collect_job_metadata(tmp_path / "missing.json", tmp_path / "also_missing.json", "job1") is None

def test_collect_job_metadata_handles_current_only(tmp_path: Path) -> None:
    current = _write_json(tmp_path / "current_metadata.json", {"status": "completed", "runtime_seconds": 10.0})
    meta = collect_job_metadata(current, None, "job1")
    assert meta.status == "completed"
    assert meta.attempt_count == 0

def test_collect_job_metadata_substitutes_real_elapsed_time_for_a_still_running_job(tmp_path: Path) -> None:
    """Regression: the report-generation rule runs as part of the same
    Snakemake DAG it's reporting on, so current_metadata.json's own
    status is always still "created" (with the original
    runtime_seconds=0.0 placeholder) at the moment this is read --
    the terminal "completed" entry is written only after the whole
    run, including this report, finishes. A real elapsed-so-far
    duration must be substituted, not a stale zero."""
    from datetime import datetime, timedelta
    start_time = (datetime.now() - timedelta(seconds = 30)).isoformat()
    current = _write_json(tmp_path / "current_metadata.json", {
        "status": "created", "start_time": start_time, "runtime_seconds": 0.0,
    })
    meta = collect_job_metadata(current, None, "job1")
    assert meta.status == "created"
    assert meta.runtime_seconds >= 29.0
