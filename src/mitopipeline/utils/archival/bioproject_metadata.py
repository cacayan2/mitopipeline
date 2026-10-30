"""
bioproject_metadata.py

Prepares the real project-level BioProject metadata describing this
pipeline run's own research effort -- one record for the whole run,
never per-sample (a BioProject represents the research project
itself; BioSample is where the per-specimen detail belongs). Never
invents an accession, a release date, or a funding/publication
record this project genuinely doesn't have on file; any such field
is left explicitly blank and flagged, matching this project's own
BioSample/GenBank convention of tracking missing information rather
than guessing it.
"""

# Imports
from __future__ import annotations
from dataclasses import asdict, dataclass
import json
from logging import Logger
from pathlib import Path

# The real NCBI BioProject "Project Data Type" controlled-vocabulary
# value this project's own output actually is (genome assembly from
# raw sequencing reads) -- never a placeholder.
DEFAULT_DATA_TYPE = "Genome sequencing"


def _resolve_sample_scope(organisms: set[str], sample_count: int) -> str:
    """
    Resolves the real NCBI "Sample Scope" controlled-vocabulary value
    from what this run actually contains -- never a fixed default
    asserted regardless of the real data. "Monoisolate" (exactly one
    real sample, one organism), "Multiisolate" (more than one real
    sample, all the same organism), or "Multispecies" (more than one
    real organism present).

    Args:
        organisms (set[str]): The real, distinct organism names observed across this run's samples.
        sample_count (int): How many real samples this run actually produced data for.

    Returns:
        str: "Monoisolate", "Multiisolate", or "Multispecies".
    """
    if len(organisms) > 1:
        return "Multispecies"
    return "Monoisolate" if sample_count <= 1 else "Multiisolate"


@dataclass(frozen = True, slots = True)
class BioProjectMetadata:
    """
    One run's real project-level BioProject metadata.

    Attributes:
        title (str): The real project title.
        description (str): The real project description.
        data_type (str): The real NCBI "Project Data Type" value.
        sample_scope (str): The real NCBI "Sample Scope" value.
        organism (str): The real organism or taxonomic branch this project covers.
        sample_ids (str): A semicolon-joined list of every real sample id this run actually produced data for.
        release_status (str): "HUP" (hold until publication) by default -- never a fabricated public release date.
        publication (str): A real citation, or "" when none exists yet (never invented).
        funding (str): Real funding/grant information, or "" when none is on file (never invented).
        missing_metadata (str): A semicolon-joined list of fields genuinely missing, for honest tracking (never silently blank).
    """
    title: str
    description: str
    data_type: str
    sample_scope: str
    organism: str
    sample_ids: str
    release_status: str
    publication: str
    funding: str
    missing_metadata: str

    def to_dict(self) -> dict:
        """Serializes this metadata to a JSON-safe dictionary."""
        return asdict(self)


def build_bioproject_metadata(sample_organisms: dict[str, str], title: str | None = None,
                               description: str | None = None, data_type: str = DEFAULT_DATA_TYPE,
                               publication: str = "", funding: str = "",
                               logger: Logger | None = None) -> BioProjectMetadata:
    """
    Builds this run's real project-level BioProject metadata from the
    real organism each real sample actually resolved to -- never an
    assumed or generic taxonomic claim.

    Args:
        sample_organisms (dict[str, str]): {sample_id: organism}, the real organism each real sample this run produced data for actually resolved to.
        title (str | None): The real project title. None uses a real, non-fabricated default built from the real organism(s) present.
        description (str | None): The real project description. None uses this project's own real, documented purpose.
        data_type (str): The real NCBI "Project Data Type" value. Defaults to "Genome sequencing".
        publication (str): A real citation, if one exists. Defaults to "" (none on file).
        funding (str): Real funding/grant information, if any is on file. Defaults to "" (none on file).
        logger (Logger | None): Logger to use. Defaults to None.

    Returns:
        BioProjectMetadata: The real, typed project-level record.
    """
    sample_ids = sorted(sample_organisms)
    organisms = {organism for organism in sample_organisms.values() if organism}
    organism_label = "; ".join(sorted(organisms)) or "unknown organism (not yet resolved)"
    sample_scope = _resolve_sample_scope(organisms, len(sample_ids))

    resolved_title = title or f"Mitochondrial genome assembly and annotation of {organism_label}"
    resolved_description = description or (
        f"Automated assembly, annotation, and comparative analysis of complete mitochondrial genomes "
        f"from {organism_label}, generated by the MitoPipeline pipeline from raw paired-end sequencing reads."
    )

    missing: list[str] = []
    if not publication:
        missing.append("publication")
    if not funding:
        missing.append("funding")
    if not sample_ids:
        missing.append("sample_ids")
    if not organisms:
        missing.append("organism")

    if logger is not None:
        logger.info(f"BioProject metadata prepared for {len(sample_ids)} sample(s), scope={sample_scope}"
                   f"{' (missing: ' + ';'.join(missing) + ')' if missing else ''}.")

    return BioProjectMetadata(
        title = resolved_title, description = resolved_description, data_type = data_type,
        sample_scope = sample_scope, organism = organism_label, sample_ids = ";".join(sample_ids),
        # Never a fabricated public release date -- "HUP" (hold until
        # publication) is NCBI's own real default for exactly this
        # situation: data prepared ahead of a release decision only
        # the investigator can actually make.
        release_status = "HUP", publication = publication, funding = funding,
        missing_metadata = ";".join(missing),
    )


def write_bioproject_metadata(metadata: BioProjectMetadata, output_path: str | Path,
                               logger: Logger | None = None) -> None:
    """
    Writes this run's real BioProject metadata as JSON.

    Args:
        metadata (BioProjectMetadata): The record built by build_bioproject_metadata().
        output_path (str | Path): Where to write the metadata JSON.
        logger (Logger | None): Logger to use. Defaults to None.

    Returns:
        None
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents = True, exist_ok = True)
    output_path.write_text(json.dumps(metadata.to_dict(), indent = 4), encoding = "utf-8")
    if logger is not None:
        logger.info(f"Wrote BioProject metadata to {output_path}.")
