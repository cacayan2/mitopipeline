"""
test_data_aggregation.py

Unit tests for consolidating each sample's own final assembly/
annotation/phylogeny outputs into aggregated_data/ (Emil's own direct
ask, 2026-10-28).
"""

# Imports
from pathlib import Path
from mitopipeline.utils.aggregation.data_aggregation import aggregate_archival_data, aggregate_sample_data

def test_aggregate_sample_data_copies_every_real_file(tmp_path: Path) -> None:
    assembly_fasta = tmp_path / "source" / "sample1.data.fasta"
    assembly_fasta.parent.mkdir(parents = True)
    assembly_fasta.write_text(">sample1\nACGT\n")

    annotation_image = tmp_path / "source" / "sample1.circular_map.png"
    annotation_image.write_bytes(b"fake png bytes")

    annotation_gff = tmp_path / "source" / "result.gff"
    annotation_gff.write_text("##gff-version 3\n")

    newick = tmp_path / "source" / "sample1.treefile.nwk"
    newick.write_text("(a,b,c);\n")

    tree_image = tmp_path / "source" / "tree_figure.png"
    tree_image.write_bytes(b"fake tree png bytes")

    output_root = tmp_path / "aggregated_data"
    result = aggregate_sample_data(
        "sample1", output_root,
        assembly_fasta = assembly_fasta, annotation_image_png = annotation_image,
        annotation_gff = annotation_gff, phylogeny_newick = newick, phylogeny_image_png = tree_image,
    )

    assert result.assembly_path == output_root / "assemblies" / "sample1.assembly.fasta"
    assert result.assembly_path.is_file()
    assert result.assembly_path.read_text() == ">sample1\nACGT\n"

    assert result.annotation_image_path == output_root / "annotations" / "sample1.annotated_map.png"
    assert result.annotation_image_path.is_file()

    assert result.annotation_gff_path == output_root / "annotations" / "gff" / "sample1.annotation.gff"
    assert result.annotation_gff_path.is_file()

    assert result.phylogeny_newick_path == output_root / "phylogenies" / "sample1.treefile.nwk"
    assert result.phylogeny_newick_path.is_file()

    assert result.phylogeny_image_path == output_root / "phylogenies" / "sample1.tree_figure.png"
    assert result.phylogeny_image_path.is_file()

def test_aggregate_sample_data_skips_missing_inputs_gracefully(tmp_path: Path) -> None:
    """A sample with no assembly (and therefore no real annotation/
    phylogeny either) must not crash -- every field is simply None."""
    result = aggregate_sample_data("sample_bare", tmp_path / "aggregated_data")

    assert result.assembly_path is None
    assert result.annotation_image_path is None
    assert result.annotation_gff_path is None
    assert result.phylogeny_newick_path is None
    assert result.phylogeny_image_path is None

def test_aggregate_sample_data_skips_a_source_that_does_not_exist(tmp_path: Path) -> None:
    result = aggregate_sample_data(
        "sample1", tmp_path / "aggregated_data",
        assembly_fasta = tmp_path / "does_not_exist.fasta",
    )
    assert result.assembly_path is None

def test_aggregate_sample_data_skips_an_empty_source_file(tmp_path: Path) -> None:
    empty_fasta = tmp_path / "empty.fasta"
    empty_fasta.write_text("")
    result = aggregate_sample_data("sample1", tmp_path / "aggregated_data", assembly_fasta = empty_fasta)
    assert result.assembly_path is None

def test_aggregate_sample_data_copies_the_real_genbank_submission_folder(tmp_path: Path) -> None:
    submission_folder = tmp_path / "source" / "sample1_Danio_rerio"
    submission_folder.mkdir(parents = True)
    (submission_folder / "sample1.fasta").write_text(">sample1\nACGT\n")
    (submission_folder / "sample1.gff").write_text("##gff-version 3\n")

    output_root = tmp_path / "aggregated_data"
    result = aggregate_sample_data("sample1", output_root, genbank_submission_dir = submission_folder)

    assert result.genbank_submission_dir_path == output_root / "archival" / "sample1" / "genbank"
    assert (result.genbank_submission_dir_path / "sample1.fasta").read_text() == ">sample1\nACGT\n"
    assert (result.genbank_submission_dir_path / "sample1.gff").is_file()

def test_aggregate_sample_data_skips_an_empty_genbank_submission_folder(tmp_path: Path) -> None:
    empty_folder = tmp_path / "source" / "sample1_empty"
    empty_folder.mkdir(parents = True)
    result = aggregate_sample_data("sample1", tmp_path / "aggregated_data", genbank_submission_dir = empty_folder)
    assert result.genbank_submission_dir_path is None

def test_aggregate_sample_data_copies_table2asn_files_into_the_same_genbank_folder(tmp_path: Path) -> None:
    sqn = tmp_path / "source" / "sample1.sqn"
    sqn.parent.mkdir(parents = True)
    sqn.write_text("fake sqn")
    val = tmp_path / "source" / "sample1.val"
    val.write_text("fake val")
    stats = tmp_path / "source" / "sample1.stats"
    stats.write_text("fake stats")
    dr = tmp_path / "source" / "sample1.dr"
    dr.write_text("fake dr")
    gbf = tmp_path / "source" / "sample1.gbf"
    gbf.write_text("fake gbf")

    output_root = tmp_path / "aggregated_data"
    result = aggregate_sample_data(
        "sample1", output_root,
        table2asn_sqn = sqn, table2asn_validation_file = val, table2asn_stats = stats,
        table2asn_discrepancy = dr, table2asn_genbank_flatfile = gbf,
    )

    genbank_dir = output_root / "archival" / "sample1" / "genbank"
    assert result.table2asn_sqn_path == genbank_dir / "sample1.sqn"
    assert result.table2asn_validation_path == genbank_dir / "sample1.val"
    assert result.table2asn_stats_path == genbank_dir / "sample1.stats"
    assert result.table2asn_discrepancy_path == genbank_dir / "sample1.dr"
    assert result.table2asn_genbank_flatfile_path == genbank_dir / "sample1.gbf"
    assert result.table2asn_sqn_path.read_text() == "fake sqn"

def test_aggregate_sample_data_lands_table2asn_files_alongside_the_genbank_submission_folder(tmp_path: Path) -> None:
    """The real intended use: table2asn files end up in the exact same directory as the GenBank submission folder's own copied fasta/gff/tbl, not a second, separate location."""
    submission_folder = tmp_path / "source" / "sample1_Danio_rerio"
    submission_folder.mkdir(parents = True)
    (submission_folder / "sample1.fasta").write_text(">sample1\nACGT\n")
    (submission_folder / "sample1.tbl").write_text(">Feature sample1\n")
    sqn = tmp_path / "source" / "sample1.sqn"
    sqn.write_text("fake sqn")

    output_root = tmp_path / "aggregated_data"
    result = aggregate_sample_data(
        "sample1", output_root, genbank_submission_dir = submission_folder, table2asn_sqn = sqn,
    )

    assert result.genbank_submission_dir_path == result.table2asn_sqn_path.parent
    assert (result.genbank_submission_dir_path / "sample1.fasta").is_file()
    assert (result.genbank_submission_dir_path / "sample1.tbl").is_file()
    assert (result.genbank_submission_dir_path / "sample1.sqn").is_file()

def test_aggregate_sample_data_skips_missing_table2asn_files_gracefully(tmp_path: Path) -> None:
    result = aggregate_sample_data("sample1", tmp_path / "aggregated_data")
    assert result.table2asn_sqn_path is None
    assert result.table2asn_validation_path is None
    assert result.table2asn_stats_path is None
    assert result.table2asn_discrepancy_path is None
    assert result.table2asn_genbank_flatfile_path is None

def test_aggregate_sample_data_isolates_different_samples(tmp_path: Path) -> None:
    fasta_a = tmp_path / "a.fasta"
    fasta_a.write_text(">a\nACGT\n")
    fasta_b = tmp_path / "b.fasta"
    fasta_b.write_text(">b\nTTTT\n")
    output_root = tmp_path / "aggregated_data"

    result_a = aggregate_sample_data("sample_a", output_root, assembly_fasta = fasta_a)
    result_b = aggregate_sample_data("sample_b", output_root, assembly_fasta = fasta_b)

    assert result_a.assembly_path.read_text() == ">a\nACGT\n"
    assert result_b.assembly_path.read_text() == ">b\nTTTT\n"
    assert result_a.assembly_path != result_b.assembly_path

def test_aggregate_archival_data_copies_every_real_job_level_file(tmp_path: Path) -> None:
    sra = tmp_path / "source" / "sra_metadata.tsv"
    sra.parent.mkdir(parents = True)
    sra.write_text("sample_name\treview_status\nsample1\tREADY\n")
    biosample = tmp_path / "source" / "biosample_metadata.tsv"
    biosample.write_text("sample_name\treview_status\nsample1\tREADY\n")
    bioproject = tmp_path / "source" / "bioproject_metadata.json"
    bioproject.write_text('{"title": "test"}')
    genbank_metadata = tmp_path / "source" / "genbank_metadata.tsv"
    genbank_metadata.write_text("sample_name\treview_status\nsample1\tREADY\n")
    validation = tmp_path / "source" / "archival_validation.tsv"
    validation.write_text("sample_id\tready_for_review\nsample1\tTrue\n")
    checklist = tmp_path / "source" / "submission_checklist.md"
    checklist.write_text("# checklist\n")
    checksums = tmp_path / "source" / "md5_checksums.tsv"
    checksums.write_text("file\tmd5\na.fasta\tabc123\n")

    output_root = tmp_path / "aggregated_data"
    result = aggregate_archival_data(
        output_root, sra_metadata_tsv = sra, biosample_metadata_tsv = biosample,
        bioproject_metadata_json = bioproject,
        genbank_metadata_tsv = genbank_metadata, archival_validation_tsv = validation,
        submission_checklist_md = checklist, checksums_tsv = checksums,
    )

    assert result.sra_metadata_path == output_root / "archival" / "sra" / "sra_metadata.tsv"
    assert result.sra_metadata_path.is_file()
    assert result.biosample_metadata_path == output_root / "archival" / "biosample" / "biosample_metadata.tsv"
    assert result.bioproject_metadata_path == output_root / "archival" / "bioproject" / "bioproject_metadata.json"
    assert result.bioproject_metadata_path.is_file()
    assert result.genbank_metadata_path == output_root / "archival" / "genbank" / "genbank_metadata.tsv"
    assert result.archival_validation_path == output_root / "archival" / "archival_validation.tsv"
    assert result.submission_checklist_path == output_root / "archival" / "submission_checklist.md"
    assert result.checksums_path == output_root / "archival" / "md5_checksums.tsv"

def test_aggregate_archival_data_skips_missing_inputs_gracefully(tmp_path: Path) -> None:
    """A run with archival disabled must not crash -- every field is simply None."""
    result = aggregate_archival_data(tmp_path / "aggregated_data")
    assert result.sra_metadata_path is None
    assert result.biosample_metadata_path is None
    assert result.bioproject_metadata_path is None
    assert result.genbank_metadata_path is None
    assert result.archival_validation_path is None
    assert result.submission_checklist_path is None
    assert result.checksums_path is None
