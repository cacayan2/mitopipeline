"""
models.py

Column schemas for NCBI archival preparation, matching
the real NCBI Submission Portal template column names for SRA,
BioSample, and GenBank/mitogenome source-modifier tables -- confirmed
against the project's own NCBI archival requirements research
(zmisc/docs/ncbi-archival-requirements.md, v0.1.0).

The pipeline never submits to NCBI and never fabricates
investigator-supplied biological metadata (collection date, locality,
etc.) -- every record instead carries its own missing_fields and
review_status so a human reviewer sees exactly what still needs
filling in before upload.
"""

# Imports
READY = "READY"
MANUAL_REVIEW_REQUIRED = "MANUAL_REVIEW_REQUIRED"

# SRA Submission Portal's own real column names.
SRA_COLUMNS = [
    "sample_name", "library_ID", "title", "library_strategy", "library_source",
    "library_selection", "library_layout", "platform", "instrument_model",
    "design_description", "filetype", "filename", "filename2", "organism",
    "missing_fields", "review_status",
]
SRA_REQUIRED_FIELDS = [
    "sample_name", "library_ID", "title", "library_strategy", "library_source",
    "library_selection", "library_layout", "platform", "instrument_model",
    "filetype", "filename",
]

# BioSample Submission Portal's own real column names.
BIOSAMPLE_COLUMNS = [
    "sample_name", "bioproject_accession", "organism", "isolate", "host",
    "isolation_source", "collection_date", "geo_loc_name", "lat_lon",
    "sample_type", "tissue", "sex", "dev_stage", "description",
    "missing_fields", "review_status",
]
BIOSAMPLE_REQUIRED_FIELDS = [
    "sample_name", "organism", "collection_date", "geo_loc_name", "isolation_source",
]

# GenBank/mitogenome submission folder contract -- an
# accession-independent FASTA + annotation set per accepted assembly,
# ready for table2asn or the Submission Portal's own genome upload.
GENBANK_COLUMNS = [
    "sample_id", "organism", "submission_folder", "assembly_fasta",
    "annotation_gff", "annotation_tbl", "sequence_length_bp", "molecule_type", "genome_location",
    "topology", "assembly_method", "annotation_method",
    "missing_files", "missing_metadata", "feature_table_warnings", "review_status",
]
GENBANK_REQUIRED_FIELDS = [
    "sample_id", "organism", "assembly_fasta", "annotation_tbl", "molecule_type",
    "genome_location", "topology", "assembly_method", "annotation_method",
]

# Cross-record consistency + submission-readiness report.
VALIDATION_COLUMNS = [
    "sample_id", "organism", "raw_fastq_ready", "assembly_ready", "annotation_ready",
    "sra_ready", "biosample_ready", "genbank_ready", "table2asn_ready", "ready_for_review",
    "missing_files", "missing_metadata", "submission_blockers",
]
