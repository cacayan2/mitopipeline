"""
genbank_metadata.py

Collects one accepted sample's real selected assembly and MITOS2
annotation into an accession-independent GenBank submission folder
(#370/#373) -- reading contig_selection's own real selection record
and MITOS2's own real GFF rather than assuming a fixed path, since
v0.2.0's selected attempt may be the primary or any rescue attempt.
"""

# Imports
import csv
import json
import shutil
from logging import Logger
from pathlib import Path
from mitopipeline.utils.archival.models import GENBANK_COLUMNS, GENBANK_REQUIRED_FIELDS, MANUAL_REVIEW_REQUIRED, READY

def _safe_folder_name(value: str) -> str:
    """Builds a filesystem-safe submission folder name from a sample id/organism string."""
    cleaned = "".join(c if c.isalnum() or c in {"-", "_", "."} else "_" for c in value.strip())
    return cleaned.strip("_") or "unknown_sample"

def _read_json(path: Path) -> dict | None:
    """Reads a JSON file, returning None when it's absent or unreadable."""
    if not path.is_file():
        return None
    try:
        return json.loads(path.read_text(encoding = "utf-8"))
    except (OSError, json.JSONDecodeError):
        return None

def build_genbank_submission(sample_id: str, organism: str, contig_selection_json: str | Path,
                              annotation_result_json: str | Path, output_directory: str | Path,
                              manifest_row: dict[str, str] | None = None,
                              defaults: dict[str, object] | None = None, mode: str = "copy",
                              logger: Logger | None = None) -> dict[str, str]:
    """
    Collects one sample's real selected assembly FASTA and MITOS2
    annotation GFF into its own accession-independent submission
    folder.

    Args:
        sample_id (str): The sample id.
        organism (str): The sample's organism binomial, for the folder name and record.
        contig_selection_json (str | Path): The selection JSON written by run_contig_selection.py.
        annotation_result_json (str | Path): The result JSON written by run_mitos2_annotation.py.
        output_directory (str | Path): Where to create this sample's own submission folder.
        manifest_row (dict[str, str] | None): This sample's own runtime manifest row -- molecule_type/genome_location/assembly_method/annotation_method, when present, win over defaults (topology is never read from here; it always reflects this sample's own real assembly status).
        defaults (dict[str, object] | None): Configured GenBank defaults (stage_options.archival.biosample -- molecule_type, genome_location, topology, assembly_method, annotation_method), used only when the manifest itself has no value.
        mode (str): "copy" or "symlink" for collecting the real files. Defaults to "copy".
        logger (Logger | None): Logger to use. Defaults to None.

    Returns:
        dict[str, str]: The GenBank submission record, in GENBANK_COLUMNS order.

    Raises:
        ValueError: If mode is not "copy" or "symlink".
    """
    if mode not in ("copy", "symlink"):
        raise ValueError(f"mode must be 'copy' or 'symlink', got {mode!r}.")

    manifest_row = manifest_row or {}
    defaults = defaults or {}
    output_directory = Path(output_directory)
    folder_name = _safe_folder_name(f"{sample_id}_{organism}" if organism else sample_id)
    submission_folder = output_directory / folder_name
    submission_folder.mkdir(parents = True, exist_ok = True)

    selection = _read_json(Path(contig_selection_json))
    annotation = _read_json(Path(annotation_result_json))

    missing_files: list[str] = []
    assembly_fasta, sequence_length_bp, topology = "", "", str(defaults.get("topology", "circular"))

    if selection is not None and selection.get("selected"):
        source_fasta = Path(selection["fasta_path"])
        if source_fasta.is_file():
            destination = submission_folder / f"{sample_id}.fasta"
            if mode == "copy":
                shutil.copy2(source_fasta, destination)
            else:
                destination.unlink(missing_ok = True)
                destination.symlink_to(source_fasta.resolve())
            assembly_fasta = str(destination)
            sequence_length_bp = str(selection.get("length_bp", ""))
            topology = "circular" if selection.get("status") == "circular" else "linear"
        else:
            missing_files.append("assembly_fasta")
    else:
        missing_files.append("assembly_fasta")

    annotation_gff = ""
    if annotation is not None and annotation.get("success") and annotation.get("gff_path"):
        source_gff = Path(annotation["gff_path"])
        if source_gff.is_file():
            destination = submission_folder / f"{sample_id}.gff"
            if mode == "copy":
                shutil.copy2(source_gff, destination)
            else:
                destination.unlink(missing_ok = True)
                destination.symlink_to(source_gff.resolve())
            annotation_gff = str(destination)
        else:
            missing_files.append("annotation_gff")
    else:
        missing_files.append("annotation_gff")

    record = {
        "sample_id": sample_id,
        "organism": organism,
        "submission_folder": str(submission_folder),
        "assembly_fasta": assembly_fasta,
        "annotation_gff": annotation_gff,
        "sequence_length_bp": sequence_length_bp,
        "molecule_type": manifest_row.get("molecule_type", "").strip() or str(defaults.get("molecule_type", "DNA")),
        "genome_location": manifest_row.get("genome_location", "").strip() or str(defaults.get("genome_location", "mitochondrion")),
        # Never manifest/config-overridable -- always the real selected
        # attempt's own assembly status, never a human-asserted claim.
        "topology": topology,
        "assembly_method": manifest_row.get("assembly_method", "").strip() or str(defaults.get("assembly_method", "GetOrganelle")),
        "annotation_method": manifest_row.get("annotation_method", "").strip() or str(defaults.get("annotation_method", "MITOS2")),
    }

    missing_metadata = [field for field in GENBANK_REQUIRED_FIELDS if not record.get(field)]
    record["missing_files"] = ";".join(sorted(set(missing_files)))
    record["missing_metadata"] = ";".join(sorted(set(missing_metadata)))
    record["review_status"] = READY if not (missing_files or missing_metadata) else MANUAL_REVIEW_REQUIRED

    if logger is not None:
        logger.info(f"{{{sample_id}}} GenBank submission folder: {record['review_status']}"
                   f"{' (' + record['missing_files'] + record['missing_metadata'] + ')' if record['review_status'] == MANUAL_REVIEW_REQUIRED else ''}.")

    return {column: record[column] for column in GENBANK_COLUMNS}

def write_genbank_metadata(records: list[dict[str, str]], output_path: str | Path, logger: Logger | None = None) -> None:
    """
    Writes a GenBank metadata TSV summarizing every sample's own
    submission folder.

    Args:
        records (list[dict[str, str]]): The records built by build_genbank_submission(), one per sample.
        output_path (str | Path): Where to write the metadata TSV.
        logger (Logger | None): Logger to use. Defaults to None.

    Returns:
        None
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents = True, exist_ok = True)
    with output_path.open("w", encoding = "utf-8", newline = "") as handle:
        writer = csv.DictWriter(handle, fieldnames = GENBANK_COLUMNS, delimiter = "\t")
        writer.writeheader()
        writer.writerows(records)

    if logger is not None: logger.info(f"Wrote GenBank metadata table to {output_path}.")
