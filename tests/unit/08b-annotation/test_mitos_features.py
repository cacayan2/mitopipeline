"""
test_mitos_features.py

Unit tests for parsing MITOS2's own GFF3 gene-feature output (#410):
coordinate/strand/qualifier normalization, and missing/malformed
feature handling.
"""

# Imports
from pathlib import Path
import pytest
from mitopipeline.utils.annotation.mitos_features import parse_mitos_gff
from mitopipeline.utils.annotation.models import GeneFeature

# Captured verbatim from a real runmitos.py run against the installed
# MITOS2 binary (#409's real end-to-end verification) -- a genuine
# result.gff, not a hand-built fixture.
REAL_GFF = (
    "##gff-version 3\n"
    "#!gff-spec-version 1.21\n"
    "test_mitogenome\tmitos\tregion\t1\t808\t.\t+\t.\tID=test_mitogenome:1..808;Is_circular=True;"
    "Name=test_mitogenome;genome=mitochondrion;mol_type=genomic DNA\n"
    "test_mitogenome\tmitos\tgene\t1\t600\t.\t+\t.\tID=gene_cox2;Name=cox2;gene_id=cox2\n"
    "test_mitogenome\tmitos\texon\t1\t600\t415800021.6\t+\t0\tParent=transcript_cox2;Name=cox2\n"
)

def test_parses_the_real_gene_row_and_skips_region_and_exon(tmp_path: Path) -> None:
    path = tmp_path / "result.gff"
    path.write_text(REAL_GFF)
    features = parse_mitos_gff(path)
    assert len(features) == 1
    feature = features[0]
    assert feature == GeneFeature(
        seqid = "test_mitogenome", gene_id = "cox2", start = 1, end = 600,
        strand = "+", score = None,
        qualifiers = {"ID": "gene_cox2", "Name": "cox2", "gene_id": "cox2"},
    )

def test_skips_header_and_blank_lines(tmp_path: Path) -> None:
    path = tmp_path / "result.gff"
    path.write_text("##gff-version 3\n\n" + REAL_GFF.splitlines()[2] + "\n")
    features = parse_mitos_gff(path)
    assert features == []   # the region row, not a gene

def test_empty_file_is_zero_genes_not_an_error(tmp_path: Path) -> None:
    path = tmp_path / "result.gff"
    path.write_text("")
    assert parse_mitos_gff(path) == []

def test_no_gene_rows_is_zero_genes_not_an_error(tmp_path: Path) -> None:
    path = tmp_path / "result.gff"
    path.write_text("##gff-version 3\n" + REAL_GFF.splitlines()[2] + "\n" + REAL_GFF.splitlines()[4] + "\n")
    assert parse_mitos_gff(path) == []

def test_missing_file_raises(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        parse_mitos_gff(tmp_path / "does_not_exist.gff")

def test_directory_instead_of_file_raises(tmp_path: Path) -> None:
    with pytest.raises(ValueError):
        parse_mitos_gff(tmp_path)

def test_wrong_column_count_raises_with_line_number(tmp_path: Path) -> None:
    path = tmp_path / "result.gff"
    path.write_text("test_mitogenome\tmitos\tgene\t1\t600\t.\t+\t.\n")   # missing attributes column
    with pytest.raises(ValueError, match = "row 1"):
        parse_mitos_gff(path)

def test_non_integer_coordinates_raise(tmp_path: Path) -> None:
    path = tmp_path / "result.gff"
    path.write_text("test_mitogenome\tmitos\tgene\tNOT_A_NUMBER\t600\t.\t+\t.\tID=gene_cox2;Name=cox2\n")
    with pytest.raises(ValueError):
        parse_mitos_gff(path)

def test_gene_row_missing_name_qualifier_raises(tmp_path: Path) -> None:
    path = tmp_path / "result.gff"
    path.write_text("test_mitogenome\tmitos\tgene\t1\t600\t.\t+\t.\tID=gene_cox2\n")
    with pytest.raises(ValueError, match = "missing Name/gene_id"):
        parse_mitos_gff(path)

def test_gene_id_falls_back_to_gene_id_qualifier_when_name_absent(tmp_path: Path) -> None:
    path = tmp_path / "result.gff"
    path.write_text("test_mitogenome\tmitos\tgene\t1\t600\t.\t+\t.\tID=gene_cox2;gene_id=cox2\n")
    features = parse_mitos_gff(path)
    assert features[0].gene_id == "cox2"

def test_unknown_strand_and_score_normalize_to_none(tmp_path: Path) -> None:
    path = tmp_path / "result.gff"
    path.write_text("test_mitogenome\tmitos\tgene\t1\t600\t.\t.\t.\tID=gene_cox2;Name=cox2\n")
    feature = parse_mitos_gff(path)[0]
    assert feature.strand is None
    assert feature.score is None

def test_real_score_value_parses_as_float(tmp_path: Path) -> None:
    path = tmp_path / "result.gff"
    path.write_text("test_mitogenome\tmitos\tgene\t1\t600\t123.45\t-\t.\tID=gene_cox2;Name=cox2\n")
    feature = parse_mitos_gff(path)[0]
    assert feature.score == 123.45
    assert feature.strand == "-"

def test_custom_feature_type_filters_differently(tmp_path: Path) -> None:
    path = tmp_path / "result.gff"
    path.write_text(REAL_GFF)
    features = parse_mitos_gff(path, feature_type = "region")
    assert len(features) == 1
    assert features[0].gene_id == "test_mitogenome"

def test_multiple_genes_preserve_file_order(tmp_path: Path) -> None:
    path = tmp_path / "result.gff"
    path.write_text(
        "test_mitogenome\tmitos\tgene\t1\t600\t.\t+\t.\tID=gene_cox2;Name=cox2\n"
        "test_mitogenome\tmitos\tgene\t650\t1200\t.\t-\t.\tID=gene_nad1;Name=nad1\n"
    )
    features = parse_mitos_gff(path)
    assert [f.gene_id for f in features] == ["cox2", "nad1"]
