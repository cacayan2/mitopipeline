"""
sra_metadata.py

Builds a reviewable NCBI SRA submission table (#368/#371) from the
runtime manifest's own raw FASTQ paths plus the configured sequencing
defaults -- manifest values always win over defaults, and a field
neither source supplies is recorded in missing_fields rather than
fabricated.
"""

# Imports
import csv
from logging import Logger
from pathlib import Path
from mitopipeline.utils.archival.models import SRA_COLUMNS, SRA_REQUIRED_FIELDS, MANUAL_REVIEW_REQUIRED, READY

def _read_manifest_rows(runtime_manifest_path: Path) -> list[dict[str, str]]:
    """Reads the runtime manifest TSV into a list of row dicts."""
    with runtime_manifest_path.open("r", encoding = "utf-8", newline = "") as handle:
        return list(csv.DictReader(handle, delimiter = "\t"))

def _organism_name(row: dict[str, str]) -> str:
    """Builds an organism binomial from organism/scientific_name, or genus+species."""
    for column in ("organism", "scientific_name"):
        if row.get(column, "").strip():
            return row[column].strip()
    genus, species = row.get("genus", "").strip(), row.get("species", "").strip()
    return " ".join(part for part in (genus, species) if part)

def build_sra_metadata(runtime_manifest_path: str | Path, defaults: dict[str, object] | None = None,
                        logger: Logger | None = None) -> list[dict[str, str]]:
    """
    Builds one reviewable SRA metadata row per raw sequencing sample.

    Args:
        runtime_manifest_path (str | Path): The runtime manifest TSV.
        defaults (dict[str, object] | None): Configured sequencing defaults (stage_options.archival.sra), used only when the manifest itself has no value.
        logger (Logger | None): Logger to use. Defaults to None.

    Returns:
        list[dict[str, str]]: One row per sample, in SRA_COLUMNS order.

    Raises:
        FileNotFoundError: If the runtime manifest doesn't exist.
        ValueError: If the manifest is missing sample_id, r1, or r2.
    """
    runtime_manifest_path = Path(runtime_manifest_path)
    if not runtime_manifest_path.is_file():
        raise FileNotFoundError(f"Runtime manifest not found: {runtime_manifest_path}.")

    defaults = defaults or {}
    rows = _read_manifest_rows(runtime_manifest_path)
    if rows and not {"sample_id", "r1", "r2"} <= set(rows[0].keys()):
        raise ValueError("Runtime manifest is missing required columns: sample_id, r1, r2.")

    records: list[dict[str, str]] = []
    for row in rows:
        sample_id = row.get("sample_id", "").strip()
        if not sample_id:
            raise ValueError("Runtime manifest contains a blank sample_id.")

        library_layout = (row.get("library_layout", "").strip() or str(defaults.get("library_layout", "PAIRED"))).upper()
        title_template = str(defaults.get("title_template", "{sample_id} raw sequencing reads"))

        record = {
            "sample_name": sample_id,
            "library_ID": row.get("library_id", "").strip() or sample_id,
            "title": row.get("title", "").strip() or title_template.format(sample_id = sample_id),
            "library_strategy": row.get("library_strategy", "").strip() or str(defaults.get("library_strategy", "")),
            "library_source": row.get("library_source", "").strip() or str(defaults.get("library_source", "")),
            "library_selection": row.get("library_selection", "").strip() or str(defaults.get("library_selection", "")),
            "library_layout": library_layout,
            "platform": row.get("platform", "").strip() or str(defaults.get("platform", "")),
            "instrument_model": row.get("instrument_model", "").strip() or str(defaults.get("instrument_model") or ""),
            "design_description": row.get("design_description", "").strip() or str(defaults.get("design_description") or ""),
            "filetype": row.get("filetype", "").strip() or str(defaults.get("filetype", "fastq")),
            "filename": str(Path(row["r1"]).resolve()) if row.get("r1", "").strip() else "",
            "filename2": str(Path(row["r2"]).resolve()) if row.get("r2", "").strip() else "",
            "organism": _organism_name(row),
        }

        missing = [field for field in SRA_REQUIRED_FIELDS if not record[field]]
        if library_layout == "PAIRED" and not record["filename2"]:
            missing.append("filename2")

        record["missing_fields"] = ";".join(sorted(set(missing)))
        record["review_status"] = READY if not missing else MANUAL_REVIEW_REQUIRED
        records.append({column: record[column] for column in SRA_COLUMNS})

    if logger is not None:
        review_count = sum(1 for r in records if r["review_status"] == MANUAL_REVIEW_REQUIRED)
        logger.info(f"Built SRA metadata for {len(records)} sample(s); {review_count} require manual review.")

    return records

def write_sra_metadata(runtime_manifest_path: str | Path, output_path: str | Path,
                        defaults: dict[str, object] | None = None, logger: Logger | None = None) -> list[dict[str, str]]:
    """
    Builds and writes the SRA metadata table as a TSV.

    Args:
        runtime_manifest_path (str | Path): The runtime manifest TSV.
        output_path (str | Path): Where to write the SRA metadata TSV.
        defaults (dict[str, object] | None): Configured sequencing defaults.
        logger (Logger | None): Logger to use. Defaults to None.

    Returns:
        list[dict[str, str]]: The rows written.
    """
    records = build_sra_metadata(runtime_manifest_path, defaults = defaults, logger = logger)

    output_path = Path(output_path)
    output_path.parent.mkdir(parents = True, exist_ok = True)
    with output_path.open("w", encoding = "utf-8", newline = "") as handle:
        writer = csv.DictWriter(handle, fieldnames = SRA_COLUMNS, delimiter = "\t")
        writer.writeheader()
        writer.writerows(records)

    if logger is not None: logger.info(f"Wrote SRA metadata table to {output_path}.")
    return records
