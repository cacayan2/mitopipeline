"""
validation.py

Cross-checks the three real archival tables (#374) -- SRA, BioSample,
GenBank -- for sample-identifier consistency, and produces one
submission-readiness record per sample (#450's own packaging
validation reads this next) -- a genuine "Submission Blockers"
report, matching the project's own NCBI archival requirements
research, not a rubber-stamped pass.
"""

# Imports
import csv
from logging import Logger
from pathlib import Path
from mitopipeline.utils.archival.models import READY, VALIDATION_COLUMNS

def _read_tsv_by_sample(path: Path, key_column: str) -> dict[str, dict[str, str]]:
    """Reads a TSV, indexed by its own key_column, or {} when absent."""
    if not path.is_file():
        return {}
    with path.open("r", encoding = "utf-8", newline = "") as handle:
        return {row[key_column]: row for row in csv.DictReader(handle, delimiter = "\t") if row.get(key_column)}

def validate_cross_record_identifiers(sra_rows: dict[str, dict], biosample_rows: dict[str, dict],
                                       genbank_rows: dict[str, dict]) -> list[str]:
    """
    Checks that every sample present in one archival table is present
    in the others too -- a sample dropped from one table but not
    another is a real, reportable inconsistency.

    Args:
        sra_rows (dict[str, dict]): SRA rows, keyed by sample_name.
        biosample_rows (dict[str, dict]): BioSample rows, keyed by sample_name.
        genbank_rows (dict[str, dict]): GenBank rows, keyed by sample_id.

    Returns:
        list[str]: Human-readable inconsistency messages, empty if every table agrees.
    """
    sample_sets = {"SRA": set(sra_rows), "BioSample": set(biosample_rows), "GenBank": set(genbank_rows)}
    all_ids = set().union(*sample_sets.values())

    messages: list[str] = []
    for sample_id in sorted(all_ids):
        missing_from = [table for table, ids in sample_sets.items() if sample_id not in ids]
        if missing_from:
            messages.append(f"Sample {sample_id} is missing from: {', '.join(missing_from)}.")
    return messages

def validate_archival_readiness(sample_id: str, organism: str, raw_fastq_ready: bool, assembly_ready: bool,
                                 annotation_ready: bool, sra_row: dict[str, str] | None,
                                 biosample_row: dict[str, str] | None, genbank_row: dict[str, str] | None) -> dict[str, str]:
    """
    Builds one sample's submission-readiness record from every real
    archival table's own review_status, plus the real upstream
    pipeline artifacts (raw reads, assembly, annotation).

    Args:
        sample_id (str): The sample id.
        organism (str): The sample's organism binomial.
        raw_fastq_ready (bool): Whether both raw FASTQ files are present.
        assembly_ready (bool): Whether a selected assembly FASTA exists.
        annotation_ready (bool): Whether a MITOS2 annotation GFF exists.
        sra_row (dict[str, str] | None): This sample's own SRA record, or None if absent.
        biosample_row (dict[str, str] | None): This sample's own BioSample record, or None if absent.
        genbank_row (dict[str, str] | None): This sample's own GenBank record, or None if absent.

    Returns:
        dict[str, str]: The readiness record, in VALIDATION_COLUMNS order.
    """
    sra_ready = sra_row is not None and sra_row.get("review_status") == READY
    biosample_ready = biosample_row is not None and biosample_row.get("review_status") == READY
    genbank_ready = genbank_row is not None and genbank_row.get("review_status") == READY

    # A missing GenBank assembly/annotation file and its corresponding
    # required-metadata field (both sourced from genbank_metadata.py's
    # own record) name the same real gap -- dedupe so it doesn't read
    # twice (e.g. "assembly_fasta;assembly_fasta").
    genbank_missing_fields = dict.fromkeys(
        f for part in (genbank_row.get("missing_files", ""), genbank_row.get("missing_metadata", "")) if genbank_row and part
        for f in part.split(";")
    )

    blockers: list[str] = []
    if not raw_fastq_ready: blockers.append("raw FASTQ files missing")
    if not assembly_ready: blockers.append("assembly not available")
    if not annotation_ready: blockers.append("annotation not available")
    if sra_row is None: blockers.append("no SRA record")
    elif not sra_ready: blockers.append(f"SRA fields need review: {sra_row.get('missing_fields', '')}")
    if biosample_row is None: blockers.append("no BioSample record")
    elif not biosample_ready: blockers.append(f"BioSample fields need review: {biosample_row.get('missing_fields', '')}")
    if genbank_row is None: blockers.append("no GenBank record")
    elif not genbank_ready: blockers.append(f"GenBank needs review: {';'.join(genbank_missing_fields)}")

    ready_for_review = not blockers

    record = {
        "sample_id": sample_id, "organism": organism,
        "raw_fastq_ready": str(raw_fastq_ready), "assembly_ready": str(assembly_ready),
        "annotation_ready": str(annotation_ready), "sra_ready": str(sra_ready),
        "biosample_ready": str(biosample_ready), "genbank_ready": str(genbank_ready),
        "ready_for_review": str(ready_for_review),
        "missing_files": ";".join(f for f in [
            "" if raw_fastq_ready else "raw_fastq", "" if assembly_ready else "assembly",
            "" if annotation_ready else "annotation",
        ] if f),
        # Every real missing SRA/BioSample/GenBank field name for this
        # sample, in one place -- the single source render_submission_checklist()
        # reads to attach a field_hints.py hint per name.
        "missing_metadata": ";".join(dict.fromkeys(f for f in [
            *((sra_row.get("missing_fields", "") if sra_row else "").split(";")),
            *((biosample_row.get("missing_fields", "") if biosample_row else "").split(";")),
            *genbank_missing_fields,
        ] if f)),
        "submission_blockers": "; ".join(blockers),
    }
    return {column: record[column] for column in VALIDATION_COLUMNS}

def build_archival_validation(runtime_manifest_path: str | Path, sra_tsv: str | Path, biosample_tsv: str | Path,
                               genbank_tsv: str | Path, job_directory: str | Path,
                               logger: Logger | None = None) -> list[dict[str, str]]:
    """
    Builds the full cross-sample archival readiness report.

    Args:
        runtime_manifest_path (str | Path): The runtime manifest TSV (the sample list).
        sra_tsv (str | Path): The SRA metadata table written by run_sra_metadata.py.
        biosample_tsv (str | Path): The BioSample metadata table written by run_biosample_metadata.py.
        genbank_tsv (str | Path): The GenBank metadata table written by run_mitogenome_metadata.py.
        job_directory (str | Path): The pipeline job's own output directory, used to check raw FASTQ/assembly/annotation presence directly.
        logger (Logger | None): Logger to use. Defaults to None.

    Returns:
        list[dict[str, str]]: One readiness record per sample.
    """
    runtime_manifest_path = Path(runtime_manifest_path)
    sra_rows = _read_tsv_by_sample(Path(sra_tsv), "sample_name")
    biosample_rows = _read_tsv_by_sample(Path(biosample_tsv), "sample_name")
    genbank_rows = _read_tsv_by_sample(Path(genbank_tsv), "sample_id")

    inconsistencies = validate_cross_record_identifiers(sra_rows, biosample_rows, genbank_rows)
    if logger is not None:
        for message in inconsistencies:
            logger.warning(message)

    records: list[dict[str, str]] = []
    with runtime_manifest_path.open("r", encoding = "utf-8", newline = "") as handle:
        for row in csv.DictReader(handle, delimiter = "\t"):
            sample_id = row.get("sample_id", "").strip()
            if not sample_id:
                continue

            sra_row, biosample_row, genbank_row = sra_rows.get(sample_id), biosample_rows.get(sample_id), genbank_rows.get(sample_id)
            organism = (genbank_row or {}).get("organism") or (biosample_row or {}).get("organism") or (sra_row or {}).get("organism", "")

            record = validate_archival_readiness(
                sample_id = sample_id, organism = organism,
                raw_fastq_ready = bool(row.get("r1")) and Path(row.get("r1", "")).is_file()
                                 and bool(row.get("r2")) and Path(row.get("r2", "")).is_file(),
                assembly_ready = bool(genbank_row) and bool(genbank_row.get("assembly_fasta")),
                annotation_ready = bool(genbank_row) and bool(genbank_row.get("annotation_gff")),
                sra_row = sra_row, biosample_row = biosample_row, genbank_row = genbank_row,
            )
            records.append(record)

    if logger is not None:
        ready_count = sum(1 for r in records if r["ready_for_review"] == "True")
        logger.info(f"Archival validation: {ready_count}/{len(records)} sample(s) ready for submission review.")

    return records

def write_archival_validation(records: list[dict[str, str]], output_path: str | Path, logger: Logger | None = None) -> None:
    """
    Writes the archival validation report as a TSV.

    Args:
        records (list[dict[str, str]]): The records built by build_archival_validation().
        output_path (str | Path): Where to write the report.
        logger (Logger | None): Logger to use. Defaults to None.

    Returns:
        None
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents = True, exist_ok = True)
    with output_path.open("w", encoding = "utf-8", newline = "") as handle:
        writer = csv.DictWriter(handle, fieldnames = VALIDATION_COLUMNS, delimiter = "\t")
        writer.writeheader()
        writer.writerows(records)

    if logger is not None: logger.info(f"Wrote archival validation report to {output_path}.")
