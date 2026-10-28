"""
field_hints.py

Human-facing hints for every NCBI SRA/BioSample/GenBank field this
pipeline can report as missing -- each hint names the field's real
expected format/example and exactly where a human supplies it (a
runtime-manifest column, which always takes priority per sample, or
a stage_options.archival config default applied to the whole run).
"""

FIELD_HINTS: dict[str, str] = {
    # SRA (stage_options.archival.sra / manifest columns)
    "sample_name": "Taken from the manifest's own `sample_id` column.",
    "library_ID": "Defaults to the sample id. Add a `library_id` manifest column for a different value.",
    "title": "A short description of the sequencing run. Add a `title` manifest column, or rely on the default `stage_options.archival.sra.title_template`.",
    "library_strategy": "e.g. `WGS`. Add a `library_strategy` manifest column, or set `stage_options.archival.sra.library_strategy`.",
    "library_source": "e.g. `GENOMIC`. Add a `library_source` manifest column, or set `stage_options.archival.sra.library_source`.",
    "library_selection": "e.g. `RANDOM`. Add a `library_selection` manifest column, or set `stage_options.archival.sra.library_selection`.",
    "library_layout": "`PAIRED` or `SINGLE`. Add a `library_layout` manifest column, or set `stage_options.archival.sra.library_layout`.",
    "platform": "e.g. `ILLUMINA`, `OXFORD_NANOPORE`. Add a `platform` manifest column, or set `stage_options.archival.sra.platform`.",
    "instrument_model": "The real sequencer model used, e.g. `Illumina NovaSeq 6000`. Add an `instrument_model` manifest column (can vary per sample), or set `stage_options.archival.sra.instrument_model` if every sample used the same instrument.",
    "filetype": "e.g. `fastq`. Add a `filetype` manifest column, or set `stage_options.archival.sra.filetype`.",
    "filename": "The raw R1 FASTQ path. Check the manifest's own `r1` column points to a real, existing file.",
    "filename2": "The raw R2 FASTQ path, required when `library_layout` is `PAIRED`. Check the manifest's own `r2` column points to a real, existing file.",
    # BioSample (stage_options.archival.biosample / manifest columns)
    "organism": "The scientific name (genus + species). Add `genus`/`species` (or a combined `organism`) manifest column.",
    "isolate": "Defaults to the sample id. Add an `isolate` manifest column for a different specimen identifier.",
    "host": "The host organism, for a parasite/symbiont sample. Add a `host` manifest column, or leave blank for a free-living organism.",
    "isolation_source": "The real environment the specimen was collected from, e.g. `freshwater river sediment`. Add an `isolation_source` manifest column (per sample), or set `stage_options.archival.biosample.isolation_source` if the whole run shares one source.",
    "collection_date": "The real specimen collection date, ISO 8601 (`YYYY-MM-DD`) or NCBI's own partial-date/`missing: not collected` placeholders. Add a `collection_date` manifest column (per sample), or set `stage_options.archival.biosample.collection_date`.",
    "geo_loc_name": "NCBI's controlled format: `Country: Region, Locality` (e.g. `USA: Illinois, Chicago`). Add a `geo_loc_name` manifest column (per sample), or set `stage_options.archival.biosample.geo_loc_name`.",
    "lat_lon": "Decimal degrees, e.g. `41.8781 N 87.6298 W`. Add a `lat_lon` manifest column, or set `stage_options.archival.biosample.lat_lon`.",
    "bioproject_accession": "An existing NCBI BioProject accession (e.g. `PRJNA123456`), registered separately through the NCBI Submission Portal before these records can reference it. Add a `bioproject_accession` manifest column, or set `stage_options.archival.biosample.bioproject_accession`.",
    "sample_type": "e.g. `environmental sample`. Add a `sample_type` manifest column, or set `stage_options.archival.biosample.sample_type`.",
    "tissue": "The tissue sampled, if applicable. Add a `tissue` manifest column, or set `stage_options.archival.biosample.tissue`.",
    "sex": "The specimen's sex, if known/applicable. Add a `sex` manifest column, or set `stage_options.archival.biosample.sex`.",
    "dev_stage": "The specimen's developmental stage, if known. Add a `dev_stage` manifest column, or set `stage_options.archival.biosample.dev_stage`.",
    "description": "A free-text sample description. Add a `description` manifest column, or rely on the default `stage_options.archival.biosample.description_template`.",
    # GenBank (molecule/genome fields reuse stage_options.archival.biosample; the rest are pipeline outputs)
    "sample_id": "Taken from the manifest's own `sample_id` column.",
    "assembly_fasta": "This sample's own selected mitochondrial assembly. Missing means the assembly stage did not produce a usable genome for this sample -- see its own Assembly Results section, not a config value to set.",
    "annotation_gff": "This sample's own MITOS2 annotation. Missing means annotation did not run or did not succeed for this sample -- see its own Annotation Results section, not a config value to set.",
    "molecule_type": "e.g. `DNA`. Add a `molecule_type` manifest column (per sample), or set `stage_options.archival.biosample.molecule_type`.",
    "genome_location": "e.g. `mitochondrion`. Add a `genome_location` manifest column (per sample), or set `stage_options.archival.biosample.genome_location`.",
    "topology": "`circular` or `linear` -- always derived automatically from this sample's own real assembly status, never a manifest column or config value to set.",
    "assembly_method": "The real assembly tool/version used. Add an `assembly_method` manifest column (per sample), or set `stage_options.archival.biosample.assembly_method`.",
    "annotation_method": "The real annotation tool/version used. Add an `annotation_method` manifest column (per sample), or set `stage_options.archival.biosample.annotation_method`.",
}

RERUN_NOTE = (
    "Add these values to the runtime manifest (preferred -- one row per sample) or the "
    "matching `stage_options.archival.*` config default, then re-run the pipeline. Only the "
    "affected archival and report steps regenerate; already-completed assembly, annotation, "
    "and other expensive stages are not re-run."
)

def field_hint(field_name: str) -> str:
    """
    Returns a field's own real hint describing its expected format and
    exactly where to supply it, or a generic fallback for any field
    name not yet catalogued here.

    Args:
        field_name (str): The real SRA/BioSample/GenBank column name.

    Returns:
        str: The hint text.
    """
    return FIELD_HINTS.get(
        field_name,
        "Check the runtime manifest or the relevant `stage_options.archival.*` config block for this field.",
    )
