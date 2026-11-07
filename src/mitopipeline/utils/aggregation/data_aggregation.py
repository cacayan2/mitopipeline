"""
data_aggregation.py

Consolidates each sample's own final assembly/annotation/phylogeny
outputs -- scattered across their own per-stage directories -- into
one flat, easy-to-browse aggregated_data/ tree (Emil's own direct
ask, 2026-10-28): assemblies/, annotations/ (+ annotations/gff/), and
phylogenies/. Copies real files from their own real, already-produced
locations; never regenerates, re-derives, or fabricates anything.

Also consolidates the run's real NCBI archival materials (Emil's own
direct follow-up ask) into an archival/ subdirectory -- organized per
organism first, then by which record type it belongs to (Emil's own
direct ask, 2026-11-06, replacing an earlier database-type-first
layout): archival/<sample_id>/genbank/, /sra/, /biosample/,
/bioproject/, so a human reviewer doing the actual manual NCBI upload
finds one self-contained folder per specimen rather than hunting the
same sample's records across four separate type folders. SRA and
BioSample are genuinely built as one combined table per run (not one
file per sample), so each organism's copy is a real extraction of
that sample's own row, not a derived value; BioProject is genuinely
one record for the whole submission, not per organism, so every
sample's copy is identical by design -- the canonical copy still
lives at the run level (archival/bioproject/), alongside the other
run-level, genuinely cross-cutting materials (the combined GenBank
metadata table, the submission checklist/readiness report/checksums)
that don't belong to any one organism.
"""

# Imports
import csv
import json
import shutil
from dataclasses import dataclass
from logging import Logger
from pathlib import Path

@dataclass(frozen = True, slots = True)
class SampleAggregationResult:
    """
    What got copied into aggregated_data/ for one sample.

    Attributes:
        sample_id (str): The sample id.
        assembly_path (Path | None): Where the assembly FASTA was copied to, or None if no assembly was available.
        annotation_image_path (Path | None): Where the annotated genome map PNG was copied to, or None if unavailable.
        annotation_gff_path (Path | None): Where the annotation GFF was copied to, or None if unavailable.
        phylogeny_newick_path (Path | None): Where the Newick tree was copied to, or None if unavailable.
        phylogeny_image_path (Path | None): Where the tree figure PNG was copied to, or None if unavailable.
        genbank_submission_dir_path (Path | None): Where this sample's GenBank submission folder was copied to, or None if unavailable.
        table2asn_sqn_path (Path | None): Where the real table2asn .sqn submission file was copied to, or None if unavailable.
        table2asn_validation_path (Path | None): Where the real table2asn .val validation report was copied to, or None if unavailable.
        table2asn_stats_path (Path | None): Where the real table2asn .stats file was copied to, or None if unavailable.
        table2asn_discrepancy_path (Path | None): Where the real table2asn .dr discrepancy report was copied to, or None if unavailable.
        table2asn_genbank_flatfile_path (Path | None): Where the real table2asn .gbf GenBank flatfile was copied to, or None if unavailable.
        sra_metadata_path (Path | None): Where this sample's own single-row SRA metadata extract was written, or None if unavailable.
        biosample_metadata_path (Path | None): Where this sample's own single-row BioSample metadata extract was written, or None if unavailable.
        bioproject_metadata_path (Path | None): Where this sample's own reference copy of the run's one real BioProject record was copied to, or None if unavailable -- BioProject is genuinely one record per whole submission, not per organism, so every sample's copy is identical by design.
    """
    sample_id: str
    assembly_path: Path | None = None
    annotation_image_path: Path | None = None
    annotation_gff_path: Path | None = None
    phylogeny_newick_path: Path | None = None
    phylogeny_image_path: Path | None = None
    genbank_submission_dir_path: Path | None = None
    table2asn_sqn_path: Path | None = None
    table2asn_validation_path: Path | None = None
    table2asn_stats_path: Path | None = None
    table2asn_discrepancy_path: Path | None = None
    table2asn_genbank_flatfile_path: Path | None = None
    sra_metadata_path: Path | None = None
    biosample_metadata_path: Path | None = None
    bioproject_metadata_path: Path | None = None

@dataclass(frozen = True, slots = True)
class ArchivalAggregationResult:
    """
    What got copied into aggregated_data/archival/ for the run's own
    job-level NCBI archival materials.

    Attributes:
        sra_metadata_path (Path | None): Where the SRA metadata table was copied to, or None if unavailable.
        biosample_metadata_path (Path | None): Where the BioSample metadata table was copied to, or None if unavailable.
        bioproject_metadata_path (Path | None): Where the BioProject metadata record was copied to, or None if unavailable.
        genbank_metadata_path (Path | None): Where the combined GenBank metadata table was copied to, or None if unavailable.
        archival_validation_path (Path | None): Where the per-sample readiness report was copied to, or None if unavailable.
        submission_checklist_path (Path | None): Where the human-readable submission checklist was copied to, or None if unavailable.
        checksums_path (Path | None): Where the MD5 checksum manifest was copied to, or None if unavailable.
    """
    sra_metadata_path: Path | None = None
    biosample_metadata_path: Path | None = None
    bioproject_metadata_path: Path | None = None
    genbank_metadata_path: Path | None = None
    archival_validation_path: Path | None = None
    submission_checklist_path: Path | None = None
    checksums_path: Path | None = None

def _copy_if_present(source: str | Path | None, destination: Path, logger: Logger | None = None) -> Path | None:
    """Copies source to destination if source is a real, existing, non-empty file. Returns destination, or None if source was unavailable."""
    if source is None:
        return None
    source = Path(source)
    if not source.is_file() or source.stat().st_size == 0:
        return None
    destination.parent.mkdir(parents = True, exist_ok = True)
    shutil.copy2(source, destination)
    if logger is not None:
        logger.info(f"Copied {source} to {destination}.")
    return destination

def _copytree_if_present(source: str | Path | None, destination: Path, logger: Logger | None = None) -> Path | None:
    """Copies source (a whole directory) to destination if source is a real, existing, non-empty directory. Returns destination, or None if source was unavailable."""
    if source is None:
        return None
    source = Path(source)
    if not source.is_dir() or not any(source.iterdir()):
        return None
    if destination.exists():
        shutil.rmtree(destination)
    destination.parent.mkdir(parents = True, exist_ok = True)
    shutil.copytree(source, destination)
    if logger is not None:
        logger.info(f"Copied {source} to {destination}.")
    return destination

def _copy_into_dir_if_present(source: str | Path | None, destination_dir: Path, logger: Logger | None = None) -> Path | None:
    """Copies source into destination_dir under its own real basename, if source is a real, existing, non-empty file. Used to land a sample's table2asn submission files alongside its GenBank submission folder, in the same destination directory. Returns the full destination path, or None if source was unavailable."""
    if source is None:
        return None
    source = Path(source)
    return _copy_if_present(source, destination_dir / source.name, logger)

def _extract_tsv_row_if_present(source: str | Path | None, sample_id: str, destination: Path,
                                 logger: Logger | None = None) -> Path | None:
    """
    Writes just one sample's own real row out of a combined, run-level
    SRA/BioSample metadata table (matched on that table's own
    "sample_name" column, which both writers set to the real sample_id)
    -- SRA and BioSample are built as one job-level table covering
    every sample, with no separate per-sample file to copy, so this is
    a real extraction, not a derived/fabricated value. Returns None
    when the source table, or this sample's own row in it, doesn't
    exist.
    """
    if source is None:
        return None
    source = Path(source)
    if not source.is_file() or source.stat().st_size == 0:
        return None
    with source.open(encoding = "utf-8", newline = "") as handle:
        reader = csv.DictReader(handle, delimiter = "\t")
        row = next((r for r in reader if r.get("sample_name") == sample_id), None)
    if row is None:
        return None
    destination.parent.mkdir(parents = True, exist_ok = True)
    with destination.open("w", encoding = "utf-8", newline = "") as handle:
        writer = csv.DictWriter(handle, fieldnames = list(row.keys()), delimiter = "\t")
        writer.writeheader()
        writer.writerow(row)
    if logger is not None:
        logger.info(f"Extracted {sample_id}'s own row from {source} to {destination}.")
    return destination

def aggregate_sample_data(sample_id: str, output_root: str | Path,
                           assembly_fasta: str | Path | None = None,
                           annotation_image_png: str | Path | None = None,
                           annotation_gff: str | Path | None = None,
                           phylogeny_newick: str | Path | None = None,
                           phylogeny_image_png: str | Path | None = None,
                           genbank_submission_dir: str | Path | None = None,
                           table2asn_sqn: str | Path | None = None,
                           table2asn_validation_file: str | Path | None = None,
                           table2asn_stats: str | Path | None = None,
                           table2asn_discrepancy: str | Path | None = None,
                           table2asn_genbank_flatfile: str | Path | None = None,
                           sra_metadata_tsv: str | Path | None = None,
                           biosample_metadata_tsv: str | Path | None = None,
                           bioproject_metadata_json: str | Path | None = None,
                           logger: Logger | None = None) -> SampleAggregationResult:
    """
    Copies one sample's own real final outputs into the consolidated
    aggregated_data/ tree (assemblies/, annotations/ [+
    annotations/gff/], phylogenies/, archival/<sample>/) under
    output_root. archival/<sample>/ is organized per organism first,
    then by which NCBI record type it belongs to (Emil's own direct
    ask, 2026-11-06, replacing the earlier database-type-first
    layout): genbank/ (that sample's real GenBank submission folder,
    plus its real table2asn submission-validation files copied in
    alongside it, so the full submission-ready package --
    .fasta/.gff/.tbl/.sqn/.val/.stats/.dr/.gbf -- lives in one place),
    sra/ and biosample/ (this one sample's own real row, extracted
    from the run's combined job-level table -- SRA/BioSample are built
    as one table per run, not one file per sample, so there's no
    separate per-sample source to copy; this is a genuine extraction
    of that sample's own data, not a derived value), and bioproject/
    (a reference copy of the run's one real BioProject record --
    BioProject is genuinely one record per whole submission, not per
    organism, so every sample's copy is identical by design; the
    canonical copy still lives at archival/bioproject/ via
    aggregate_archival_data() below). Any input that's None, missing,
    or empty on disk is skipped quietly -- a sample with no assembly
    (and therefore no real annotation/phylogeny/GenBank folder either)
    is a normal outcome, never a crash.

    Args:
        sample_id (str): The sample id.
        output_root (str | Path): The aggregated_data/ directory itself.
        assembly_fasta (str | Path | None): The selected attempt's own real assembly FASTA.
        annotation_image_png (str | Path | None): The sample's own real annotated genome map PNG.
        annotation_gff (str | Path | None): The sample's own real MITOS2 annotation GFF.
        phylogeny_newick (str | Path | None): The sample's own real Newick tree.
        phylogeny_image_png (str | Path | None): The sample's own real tree figure PNG.
        genbank_submission_dir (str | Path | None): The sample's own real GenBank submission folder, copied whole so it stays a self-contained, submission-ready package.
        table2asn_sqn (str | Path | None): The sample's own real table2asn .sqn submission file.
        table2asn_validation_file (str | Path | None): The sample's own real table2asn .val validation report.
        table2asn_stats (str | Path | None): The sample's own real table2asn .stats file.
        table2asn_discrepancy (str | Path | None): The sample's own real table2asn .dr discrepancy report.
        table2asn_genbank_flatfile (str | Path | None): The sample's own real table2asn .gbf GenBank flatfile.
        sra_metadata_tsv (str | Path | None): The run's own real combined SRA metadata table, to extract this one sample's row from.
        biosample_metadata_tsv (str | Path | None): The run's own real combined BioSample metadata table, to extract this one sample's row from.
        bioproject_metadata_json (str | Path | None): The run's own real BioProject metadata record, copied as-is (not per-sample data).
        logger (Logger | None): Logger to use. Defaults to None.

    Returns:
        SampleAggregationResult: What was actually copied for this sample.
    """
    output_root = Path(output_root)
    sample_archival_dir = output_root / "archival" / sample_id
    genbank_dir = sample_archival_dir / "genbank"

    assembly_path = _copy_if_present(
        assembly_fasta, output_root / "assemblies" / f"{sample_id}.assembly.fasta", logger,
    )
    annotation_image_path = _copy_if_present(
        annotation_image_png, output_root / "annotations" / f"{sample_id}.annotated_map.png", logger,
    )
    annotation_gff_path = _copy_if_present(
        annotation_gff, output_root / "annotations" / "gff" / f"{sample_id}.annotation.gff", logger,
    )
    phylogeny_newick_path = _copy_if_present(
        phylogeny_newick, output_root / "phylogenies" / f"{sample_id}.treefile.nwk", logger,
    )
    phylogeny_image_path = _copy_if_present(
        phylogeny_image_png, output_root / "phylogenies" / f"{sample_id}.tree_figure.png", logger,
    )
    genbank_submission_dir_path = _copytree_if_present(genbank_submission_dir, genbank_dir, logger)
    table2asn_sqn_path = _copy_into_dir_if_present(table2asn_sqn, genbank_dir, logger)
    table2asn_validation_path = _copy_into_dir_if_present(table2asn_validation_file, genbank_dir, logger)
    table2asn_stats_path = _copy_into_dir_if_present(table2asn_stats, genbank_dir, logger)
    table2asn_discrepancy_path = _copy_into_dir_if_present(table2asn_discrepancy, genbank_dir, logger)
    table2asn_genbank_flatfile_path = _copy_into_dir_if_present(table2asn_genbank_flatfile, genbank_dir, logger)
    sra_metadata_path = _extract_tsv_row_if_present(
        sra_metadata_tsv, sample_id, sample_archival_dir / "sra" / "sra_metadata.tsv", logger,
    )
    biosample_metadata_path = _extract_tsv_row_if_present(
        biosample_metadata_tsv, sample_id, sample_archival_dir / "biosample" / "biosample_metadata.tsv", logger,
    )
    bioproject_metadata_path = _copy_if_present(
        bioproject_metadata_json, sample_archival_dir / "bioproject" / "bioproject_metadata.json", logger,
    )

    return SampleAggregationResult(
        sample_id = sample_id, assembly_path = assembly_path,
        annotation_image_path = annotation_image_path, annotation_gff_path = annotation_gff_path,
        phylogeny_newick_path = phylogeny_newick_path, phylogeny_image_path = phylogeny_image_path,
        genbank_submission_dir_path = genbank_submission_dir_path,
        table2asn_sqn_path = table2asn_sqn_path, table2asn_validation_path = table2asn_validation_path,
        table2asn_stats_path = table2asn_stats_path, table2asn_discrepancy_path = table2asn_discrepancy_path,
        table2asn_genbank_flatfile_path = table2asn_genbank_flatfile_path,
        sra_metadata_path = sra_metadata_path, biosample_metadata_path = biosample_metadata_path,
        bioproject_metadata_path = bioproject_metadata_path,
    )

def aggregate_archival_data(output_root: str | Path,
                             sra_metadata_tsv: str | Path | None = None,
                             biosample_metadata_tsv: str | Path | None = None,
                             bioproject_metadata_json: str | Path | None = None,
                             genbank_metadata_tsv: str | Path | None = None,
                             archival_validation_tsv: str | Path | None = None,
                             submission_checklist_md: str | Path | None = None,
                             checksums_tsv: str | Path | None = None,
                             logger: Logger | None = None) -> ArchivalAggregationResult:
    """
    Copies the run's own real job-level, genuinely cross-cutting NCBI
    archival materials into aggregated_data/archival/ under
    output_root -- the combined SRA/BioSample/GenBank metadata tables,
    the canonical BioProject record, and the run-level readiness
    report/checklist/checksums -- at archival/'s own top level,
    alongside (not instead of) the per-organism archival/<sample_id>/
    folders aggregate_sample_data() builds for each specimen's own
    extracted records. Any input that's None, missing, or empty on disk is skipped quietly
    -- a run with archival disabled is a normal outcome.

    Args:
        output_root (str | Path): The aggregated_data/ directory itself.
        sra_metadata_tsv (str | Path | None): The run's own real SRA metadata table.
        biosample_metadata_tsv (str | Path | None): The run's own real BioSample metadata table.
        bioproject_metadata_json (str | Path | None): The run's own real BioProject metadata record.
        genbank_metadata_tsv (str | Path | None): The run's own real combined GenBank metadata table.
        archival_validation_tsv (str | Path | None): The run's own real per-sample submission-readiness report.
        submission_checklist_md (str | Path | None): The run's own real human-readable submission checklist.
        checksums_tsv (str | Path | None): The run's own real MD5 checksum manifest.
        logger (Logger | None): Logger to use. Defaults to None.

    Returns:
        ArchivalAggregationResult: What was actually copied.
    """
    archival_root = Path(output_root) / "archival"

    return ArchivalAggregationResult(
        sra_metadata_path = _copy_if_present(sra_metadata_tsv, archival_root / "sra" / "sra_metadata.tsv", logger),
        biosample_metadata_path = _copy_if_present(biosample_metadata_tsv, archival_root / "biosample" / "biosample_metadata.tsv", logger),
        bioproject_metadata_path = _copy_if_present(bioproject_metadata_json, archival_root / "bioproject" / "bioproject_metadata.json", logger),
        genbank_metadata_path = _copy_if_present(genbank_metadata_tsv, archival_root / "genbank" / "genbank_metadata.tsv", logger),
        archival_validation_path = _copy_if_present(archival_validation_tsv, archival_root / "archival_validation.tsv", logger),
        submission_checklist_path = _copy_if_present(submission_checklist_md, archival_root / "submission_checklist.md", logger),
        checksums_path = _copy_if_present(checksums_tsv, archival_root / "md5_checksums.tsv", logger),
    )
