"""
test_circular_map.py

Unit tests for circular mitochondrial genome map normalization and
rendering (#349-352): origin-crossing splitting, reverse-strand and
missing-annotation handling, and PNG render/error handling.
"""

# Imports
from pathlib import Path
import pytest
from mitopipeline.utils.visualization.circular_map import (
    DrawableFeature,
    build_circular_map_data,
    render_circular_genome_map,
    render_circular_map_data,
    split_origin_crossing_feature,
)
from mitopipeline.utils.annotation.models import GeneFeature

GENOME_LENGTH = 20
FASTA_TEXT = f">sample1\n{'ACGT' * (GENOME_LENGTH // 4)}\n"   # 50% GC, deterministic

# A real MITOS2 result.gff layout: one protein-coding gene (type
# "gene"), one tRNA (its own "tRNA" transcript row, not the
# "ncRNA_gene" row -- confirmed against MITOS2's own real gffstr()),
# one rRNA, and one gene wrapping the circular origin (end > genome
# length, per MITOS2's own real convention).
REAL_GFF = (
    "##gff-version 3\n"
    "#!gff-spec-version 1.21\n"
    "sample1\tmitos\tregion\t1\t20\t.\t+\t.\tID=sample1:1..20;Is_circular=True;Name=sample1\n"
    "sample1\tmitos\tgene\t2\t6\t.\t+\t.\tID=gene_cox1;Name=cox1;gene_id=cox1\n"
    "sample1\tmitos\ttRNA\t8\t10\t.\t-\t.\tID=transcript_trnL;Name=trnL;gene_id=trnL\n"
    "sample1\tmitos\trRNA\t12\t14\t.\t+\t.\tID=transcript_rrnS;Name=rrnS;gene_id=rrnS\n"
    "sample1\tmitos\tgene\t18\t25\t.\t+\t.\tID=gene_nad1;Name=nad1;gene_id=nad1\n"
)

def _write(tmp_path: Path, name: str, text: str) -> Path:
    path = tmp_path / name
    path.write_text(text)
    return path

def test_build_circular_map_data_parses_all_three_drawable_types(tmp_path: Path) -> None:
    gff = _write(tmp_path, "result.gff", REAL_GFF)
    fasta = _write(tmp_path, "genome.fasta", FASTA_TEXT)

    map_data = build_circular_map_data("sample1", gff, fasta)

    assert map_data.sample_id == "sample1"
    assert map_data.sequence_id == "sample1"
    assert map_data.genome_length == GENOME_LENGTH
    assert map_data.gc_content_percent == pytest.approx(50.0)
    gene_ids = sorted({d.feature.gene_id for d in map_data.features})
    assert gene_ids == ["cox1", "nad1", "rrnS", "trnL"]

def test_origin_crossing_feature_is_split_into_two_segments(tmp_path: Path) -> None:
    gff = _write(tmp_path, "result.gff", REAL_GFF)
    fasta = _write(tmp_path, "genome.fasta", FASTA_TEXT)

    map_data = build_circular_map_data("sample1", gff, fasta)

    nad1_segments = [d for d in map_data.features if d.feature.gene_id == "nad1"]
    assert len(nad1_segments) == 2
    starts_ends = sorted((d.feature.start, d.feature.end) for d in nad1_segments)
    assert starts_ends == [(1, 5), (18, 20)]

def test_non_crossing_feature_is_not_split() -> None:
    feature = GeneFeature(seqid = "s", gene_id = "cox1", start = 2, end = 6, strand = "+", score = None)
    segments = split_origin_crossing_feature(DrawableFeature(feature, "gene"), GENOME_LENGTH)
    assert segments == [DrawableFeature(feature, "gene")]

def test_feature_starting_beyond_genome_length_raises() -> None:
    feature = GeneFeature(seqid = "s", gene_id = "bad", start = 25, end = 30, strand = "+", score = None)
    with pytest.raises(ValueError):
        split_origin_crossing_feature(DrawableFeature(feature, "gene"), GENOME_LENGTH)

def test_reverse_strand_feature_is_preserved_through_normalization(tmp_path: Path) -> None:
    gff = _write(tmp_path, "result.gff", REAL_GFF)
    fasta = _write(tmp_path, "genome.fasta", FASTA_TEXT)

    map_data = build_circular_map_data("sample1", gff, fasta)

    trnl = next(d for d in map_data.features if d.feature.gene_id == "trnL")
    assert trnl.feature.strand == "-"
    assert trnl.feature_type == "tRNA"

def test_missing_annotations_is_an_empty_feature_set_not_an_error(tmp_path: Path) -> None:
    empty_gff = "##gff-version 3\nsample1\tmitos\tregion\t1\t20\t.\t+\t.\tID=sample1:1..20\n"
    gff = _write(tmp_path, "result.gff", empty_gff)
    fasta = _write(tmp_path, "genome.fasta", FASTA_TEXT)

    map_data = build_circular_map_data("sample1", gff, fasta)
    assert map_data.features == ()

def test_missing_gff_raises(tmp_path: Path) -> None:
    fasta = _write(tmp_path, "genome.fasta", FASTA_TEXT)
    with pytest.raises(FileNotFoundError):
        build_circular_map_data("sample1", tmp_path / "missing.gff", fasta)

def test_missing_fasta_raises(tmp_path: Path) -> None:
    gff = _write(tmp_path, "result.gff", REAL_GFF)
    with pytest.raises(FileNotFoundError):
        build_circular_map_data("sample1", gff, tmp_path / "missing.fasta")

def test_render_writes_a_real_png(tmp_path: Path) -> None:
    gff = _write(tmp_path, "result.gff", REAL_GFF)
    fasta = _write(tmp_path, "genome.fasta", FASTA_TEXT)
    output_png = tmp_path / "map.png"

    map_data = render_circular_genome_map("sample1", gff, fasta, output_png)

    assert output_png.is_file()
    assert output_png.stat().st_size > 0
    assert map_data.sample_id == "sample1"

def test_render_with_empty_features_still_produces_a_png(tmp_path: Path) -> None:
    empty_gff = "##gff-version 3\nsample1\tmitos\tregion\t1\t20\t.\t+\t.\tID=sample1:1..20\n"
    gff = _write(tmp_path, "result.gff", empty_gff)
    fasta = _write(tmp_path, "genome.fasta", FASTA_TEXT)
    output_png = tmp_path / "map.png"

    render_circular_genome_map("sample1", gff, fasta, output_png)
    assert output_png.is_file()
    assert output_png.stat().st_size > 0

def test_render_rejects_non_positive_dpi(tmp_path: Path) -> None:
    gff = _write(tmp_path, "result.gff", REAL_GFF)
    fasta = _write(tmp_path, "genome.fasta", FASTA_TEXT)
    map_data = build_circular_map_data("sample1", gff, fasta)

    with pytest.raises(ValueError):
        render_circular_map_data(map_data, tmp_path / "map.png", dpi = 0)

def test_render_rejects_invalid_genome_length(tmp_path: Path) -> None:
    gff = _write(tmp_path, "result.gff", REAL_GFF)
    fasta = _write(tmp_path, "genome.fasta", FASTA_TEXT)
    map_data = build_circular_map_data("sample1", gff, fasta)
    invalid_map_data = type(map_data)(
        sample_id = map_data.sample_id, genome_length = 0, sequence_id = map_data.sequence_id,
        gc_content_percent = map_data.gc_content_percent, features = map_data.features,
    )

    with pytest.raises(ValueError):
        render_circular_map_data(invalid_map_data, tmp_path / "map.png")
