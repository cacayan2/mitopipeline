"""
data_aggregation.py

Consolidates each sample's own final assembly/annotation/phylogeny
outputs -- scattered across their own per-stage directories -- into
one flat, easy-to-browse aggregated_data/ tree (Emil's own direct
ask, 2026-10-28): assemblies/, annotations/ (+ annotations/gff/), and
phylogenies/. Copies real files from their own real, already-produced
locations; never regenerates, re-derives, or fabricates anything.

Also consolidates the run's real NCBI archival materials (Emil's own
direct follow-up ask) into an archival/ subdirectory, one real
subfolder per metadata type -- archival/sra/, archival/biosample/,
archival/bioproject/, and archival/genbank/ (each sample's own real
GenBank submission folder, plus that sample's real table2asn
submission-validation files copied in alongside it) -- so a human
reviewer doing the actual manual NCBI upload finds each record type
already organized the way NCBI itself separates them, never a flat
pile of same-looking .tsv files. The run-level submission checklist/
readiness report/checksums span all four types at once, so they stay
at archival/'s own top level as a cross-cutting index into the
type-specific subfolders.
"""

# Imports
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
                           logger: Logger | None = None) -> SampleAggregationResult:
    """
    Copies one sample's own real final outputs into the consolidated
    aggregated_data/ tree (assemblies/, annotations/ [+
    annotations/gff/], phylogenies/, archival/genbank/<sample>/) under
    output_root. The sample's real table2asn submission-validation
    files are copied into that same archival/genbank/<sample>/
    directory, alongside its GenBank FASTA/GFF/.tbl, so the real
    submission-ready package (.fasta/.gff/.tbl/.sqn/.val/.stats/.dr/
    .gbf) lives in one place rather than split across two stage
    directories. Any input that's None, missing, or empty on disk is
    skipped quietly -- a sample with no assembly (and therefore no
    real annotation/phylogeny/GenBank folder either) is a normal
    outcome, never a crash.

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
        logger (Logger | None): Logger to use. Defaults to None.

    Returns:
        SampleAggregationResult: What was actually copied for this sample.
    """
    output_root = Path(output_root)
    genbank_dir = output_root / "archival" / "genbank" / sample_id

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

    return SampleAggregationResult(
        sample_id = sample_id, assembly_path = assembly_path,
        annotation_image_path = annotation_image_path, annotation_gff_path = annotation_gff_path,
        phylogeny_newick_path = phylogeny_newick_path, phylogeny_image_path = phylogeny_image_path,
        genbank_submission_dir_path = genbank_submission_dir_path,
        table2asn_sqn_path = table2asn_sqn_path, table2asn_validation_path = table2asn_validation_path,
        table2asn_stats_path = table2asn_stats_path, table2asn_discrepancy_path = table2asn_discrepancy_path,
        table2asn_genbank_flatfile_path = table2asn_genbank_flatfile_path,
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
    Copies the run's own real job-level NCBI archival materials into
    aggregated_data/archival/ under output_root, one real subfolder
    per metadata type -- archival/sra/, archival/biosample/,
    archival/bioproject/, archival/genbank/ (the combined metadata
    table here, alongside each sample's own GenBank submission folder
    copied per-sample by aggregate_sample_data()) -- so a human
    reviewer preparing the real manual NCBI upload finds each record
    type already separated the way NCBI itself separates them, never
    a flat pile of same-looking .tsv files. The run-level readiness
    report/checklist/checksums span all four types at once, so they
    stay at archival/'s own top level as a cross-cutting index. Any
    input that's None, missing, or empty on disk is skipped quietly
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
