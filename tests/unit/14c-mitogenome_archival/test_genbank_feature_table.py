"""
test_genbank_feature_table.py

Unit tests for building a real NCBI GenBank 5-column feature table
(.tbl) from MITOS2's own parsed gene/tRNA/rRNA features: location-
line formatting (including origin-crossing and minus-strand), real
gene/product name resolution, unrecognized-feature handling, and
vertebrate-mitochondrial-genetic-code CDS translation checks.
"""

# Imports
from pathlib import Path
import pytest
from mitopipeline.utils.annotation.models import GeneFeature
from mitopipeline.utils.archival.genbank_feature_table import (
    build_feature_table,
    build_sample_feature_table,
    check_cds_translations,
    trna_product,
)


def _gene(gene_id: str, start: int, end: int, strand: str = "+") -> GeneFeature:
    return GeneFeature(seqid = "seq1", gene_id = gene_id, start = start, end = end,
                        strand = strand, score = None, qualifiers = {})


def test_trna_product_resolves_known_codes() -> None:
    assert trna_product("trnG") == ("trnG", "tRNA-Gly")
    assert trna_product("trnL1") == ("trnL1", "tRNA-Leu")
    assert trna_product("trnS2") == ("trnS2", "tRNA-Ser")


def test_trna_product_strips_the_real_anticodon_suffix() -> None:
    # MITOS2's own "tRNA" row embeds the real anticodon in parens
    # (confirmed against a real GFF output) -- useful provenance,
    # but not a valid NCBI /gene symbol.
    assert trna_product("trnG(tcc)") == ("trnG", "tRNA-Gly")


def test_trna_product_rejects_unrecognized_code() -> None:
    with pytest.raises(ValueError):
        trna_product("trnZ")


def test_build_feature_table_writes_source_gene_and_cds_lines() -> None:
    genes = [_gene("cox2", 1, 600)]
    text, warnings = build_feature_table("sample1", "seq1", genes, [], [], genome_length = 1000,
                                           organism = "Danio rerio")
    assert warnings == []
    lines = text.splitlines()
    assert lines[0] == ">Feature seq1"
    assert "1\t1000\tsource" in lines
    assert "\t\t\torganism\tDanio rerio" in lines
    assert "\t\t\tmol_type\tgenomic DNA" in lines
    assert "1\t600\tgene" in lines
    assert "\t\t\tgene\tCOX2" in lines
    assert "1\t600\tCDS" in lines
    assert "\t\t\tproduct\tcytochrome c oxidase subunit II" in lines
    assert "\t\t\tcodon_start\t1" in lines
    assert "\t\t\ttransl_table\t2" in lines


def test_build_feature_table_writes_trna_and_rrna_lines() -> None:
    text, warnings = build_feature_table(
        "sample1", "seq1", [], [_gene("trnG", 100, 169)], [_gene("rrnS", 200, 1149)],
        genome_length = 2000, organism = "Danio rerio",
    )
    assert warnings == []
    assert "100\t169\ttRNA" in text.splitlines()
    assert "\t\t\tproduct\ttRNA-Gly" in text.splitlines()
    assert "200\t1149\trRNA" in text.splitlines()
    assert "\t\t\tproduct\t12S ribosomal RNA" in text.splitlines()


def test_minus_strand_feature_writes_reversed_location_columns() -> None:
    # GFF coordinates are always given in ascending order (100 < 500)
    # regardless of strand -- the strand column alone carries the
    # biological direction, matching real MITOS2 output.
    text, _ = build_feature_table("sample1", "seq1", [_gene("nad6", 100, 500, strand = "-")], [], [],
                                    genome_length = 1000, organism = "Danio rerio")
    # Minus strand is encoded by writing the larger coordinate first --
    # never a separate qualifier.
    assert "500\t100\tgene" in text.splitlines()


def test_origin_crossing_feature_splits_into_two_location_lines() -> None:
    # A feature MITOS2 reports as crossing the circular origin: end
    # (16620) exceeds genome_length (16596), MITOS2's own real
    # convention (confirmed against circular_map.py's identical
    # handling), never start > end.
    text, warnings = build_feature_table("sample1", "seq1", [], [_gene("trnQ", 16560, 16620)], [],
                                           genome_length = 16596, organism = "Danio rerio")
    assert warnings == []
    lines = text.splitlines()
    assert "16560\t16596\ttRNA" in lines
    assert "1\t24\t" in lines  # continuation line: no feature key, 16620 - 16596 = 24


def test_unrecognized_gene_skipped_with_warning() -> None:
    text, warnings = build_feature_table("sample1", "seq1", [_gene("mystery_gene", 1, 100)], [], [],
                                           genome_length = 1000, organism = "Danio rerio")
    assert len(warnings) == 1
    assert "mystery_gene" in warnings[0]
    assert "1\t100\tgene" not in text.splitlines()


def test_check_cds_translations_passes_a_clean_orf(tmp_path: Path) -> None:
    fasta = tmp_path / "genome.fasta"
    fasta.write_text(">seq1\nATGTAA\n")
    warnings = check_cds_translations([_gene("cox2", 1, 6)], fasta, genome_length = 6)
    assert warnings == []


def test_check_cds_translations_flags_a_real_premature_stop(tmp_path: Path) -> None:
    fasta = tmp_path / "genome.fasta"
    # ATG (Met) / TAA (stop, vertebrate mt table 2) / ATG (Met) -- a
    # real internal stop before the CDS's own last codon.
    fasta.write_text(">seq1\nATGTAAATG\n")
    warnings = check_cds_translations([_gene("cox2", 1, 9)], fasta, genome_length = 9)
    assert len(warnings) == 1
    assert "cox2" in warnings[0]
    assert "premature internal stop" in warnings[0]


def test_build_sample_feature_table_writes_a_real_file(tmp_path: Path) -> None:
    gff = tmp_path / "result.gff"
    gff.write_text(
        "##gff-version 3\n"
        "seq1\tmitos\tgene\t1\t600\t.\t+\t.\tID=gene_cox2;Name=cox2;gene_id=cox2\n"
        "seq1\tmitfi\ttRNA\t601\t670\t.\t+\t.\tID=transcript_trnG;Name=trnG;gene_id=trnG\n"
    )
    fasta = tmp_path / "genome.fasta"
    # 600bp CDS (1-600) + 70bp tRNA (601-670) = 670bp total, matching the GFF above exactly.
    fasta.write_text(">seq1\n" + "ATG" + "C" * 594 + "TAA" + "A" * 70 + "\n")
    output = tmp_path / "sample1.tbl"

    warnings = build_sample_feature_table("sample1", gff, fasta, "Danio rerio", output)

    assert output.is_file()
    text = output.read_text()
    assert ">Feature seq1" in text
    assert "1\t600\tCDS" in text.splitlines()
    assert "601\t670\ttRNA" in text.splitlines()
    assert warnings == []
