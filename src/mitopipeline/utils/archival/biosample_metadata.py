"""
biosample_metadata.py

Builds a reviewable NCBI BioSample submission table (#369/#372) from
the runtime manifest's own sample metadata plus configured
defaults -- real investigator-supplied biological fields (collection
date, locality) are never fabricated; a missing one is recorded in
missing_fields instead.
"""

# Imports
import csv
from logging import Logger
from pathlib import Path
from mitopipeline.utils.archival.models import BIOSAMPLE_COLUMNS, BIOSAMPLE_REQUIRED_FIELDS, MANUAL_REVIEW_REQUIRED, READY

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

def build_biosample_metadata(runtime_manifest_path: str | Path, defaults: dict[str, object] | None = None,
                              logger: Logger | None = None) -> list[dict[str, str]]:
    """
    Builds one reviewable BioSample metadata row per biological sample.

    Args:
        runtime_manifest_path (str | Path): The runtime manifest TSV.
        defaults (dict[str, object] | None): Configured biological defaults (stage_options.archival.biosample), used only when the manifest itself has no value.
        logger (Logger | None): Logger to use. Defaults to None.

    Returns:
        list[dict[str, str]]: One row per sample, in BIOSAMPLE_COLUMNS order.

    Raises:
        FileNotFoundError: If the runtime manifest doesn't exist.
        ValueError: If the manifest is missing sample_id, or a row has a blank sample_id.
    """
    runtime_manifest_path = Path(runtime_manifest_path)
    if not runtime_manifest_path.is_file():
        raise FileNotFoundError(f"Runtime manifest not found: {runtime_manifest_path}.")

    defaults = defaults or {}
    rows = _read_manifest_rows(runtime_manifest_path)
    if rows and "sample_id" not in rows[0]:
        raise ValueError("Runtime manifest is missing required column: sample_id.")

    records: list[dict[str, str]] = []
    for row in rows:
        sample_id = row.get("sample_id", "").strip()
        if not sample_id:
            raise ValueError("Runtime manifest contains a blank sample_id.")

        description_template = str(defaults.get("description_template", "Biological sample for {sample_id}."))

        record = {
            "sample_name": sample_id,
            "bioproject_accession": row.get("bioproject_accession", "").strip() or str(defaults.get("bioproject_accession") or ""),
            "organism": _organism_name(row),
            "isolate": row.get("isolate", "").strip() or sample_id,
            "host": row.get("host", "").strip() or str(defaults.get("host") or ""),
            "isolation_source": row.get("isolation_source", "").strip() or str(defaults.get("isolation_source") or ""),
            "collection_date": row.get("collection_date", "").strip() or str(defaults.get("collection_date") or ""),
            "geo_loc_name": row.get("geo_loc_name", "").strip() or str(defaults.get("geo_loc_name") or ""),
            "lat_lon": row.get("lat_lon", "").strip() or str(defaults.get("lat_lon") or ""),
            "sample_type": row.get("sample_type", "").strip() or str(defaults.get("sample_type", "")),
            "tissue": row.get("tissue", "").strip() or str(defaults.get("tissue") or ""),
            "sex": row.get("sex", "").strip() or str(defaults.get("sex") or ""),
            "dev_stage": row.get("dev_stage", "").strip() or str(defaults.get("dev_stage") or ""),
            "description": row.get("description", "").strip() or description_template.format(sample_id = sample_id),
        }

        missing = [field for field in BIOSAMPLE_REQUIRED_FIELDS if not record[field]]
        record["missing_fields"] = ";".join(sorted(set(missing)))
        record["review_status"] = READY if not missing else MANUAL_REVIEW_REQUIRED
        records.append({column: record[column] for column in BIOSAMPLE_COLUMNS})

    if logger is not None:
        review_count = sum(1 for r in records if r["review_status"] == MANUAL_REVIEW_REQUIRED)
        logger.info(f"Built BioSample metadata for {len(records)} sample(s); {review_count} require manual review.")

    return records

def write_biosample_metadata(runtime_manifest_path: str | Path, output_path: str | Path,
                              defaults: dict[str, object] | None = None, logger: Logger | None = None) -> list[dict[str, str]]:
    """
    Builds and writes the BioSample metadata table as a TSV.

    Args:
        runtime_manifest_path (str | Path): The runtime manifest TSV.
        output_path (str | Path): Where to write the BioSample metadata TSV.
        defaults (dict[str, object] | None): Configured biological defaults.
        logger (Logger | None): Logger to use. Defaults to None.

    Returns:
        list[dict[str, str]]: The rows written.
    """
    records = build_biosample_metadata(runtime_manifest_path, defaults = defaults, logger = logger)

    output_path = Path(output_path)
    output_path.parent.mkdir(parents = True, exist_ok = True)
    with output_path.open("w", encoding = "utf-8", newline = "") as handle:
        writer = csv.DictWriter(handle, fieldnames = BIOSAMPLE_COLUMNS, delimiter = "\t")
        writer.writeheader()
        writer.writerows(records)

    if logger is not None: logger.info(f"Wrote BioSample metadata table to {output_path}.")
    return records
