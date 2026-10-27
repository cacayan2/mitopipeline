"""
mitos_features.py

Parse MITOS2's own GFF3 gene-feature output (result.gff) into typed,
normalized GeneFeature records.
"""

# Imports
from logging import Logger
from pathlib import Path
from mitopipeline.utils.annotation.models import GeneFeature

# The real column layout of MITOS2's own result.gff (confirmed
# against a real installed-runmitos.py run): standard GFF3 columns.
GFF_COLUMNS = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]

def _parse_attributes(raw: str) -> dict[str, str]:
    """
    Parses a GFF3 attributes column ("key=value;key=value;...") into
    a dict. A malformed pair with no "=" is skipped rather than
    raising -- MITOS2's own attributes column is otherwise
    well-formed, and this is incidental metadata, not a correctness
    gate on the feature itself.

    Args:
        raw (str): The raw attributes column value.

    Returns:
        dict[str, str]: The parsed qualifiers, in column order.
    """
    qualifiers: dict[str, str] = {}
    for pair in raw.split(";"):
        pair = pair.strip()
        if not pair or "=" not in pair:
            continue
        key, _, value = pair.partition("=")
        qualifiers[key.strip()] = value.strip()
    return qualifiers

def parse_mitos_gff(gff_path: str | Path, logger: Logger | None = None, feature_type: str = "gene") -> list[GeneFeature]:
    """
    Parses MITOS2's own result.gff into typed GeneFeature records,
    normalizing coordinates, strand, and qualifiers. Only rows whose
    GFF "type" column matches feature_type are kept ("gene" by
    default -- MITOS2 also writes "region" and "exon" rows in the
    same file, which this intentionally excludes).

    Blank lines and "#"-prefixed header/comment lines (MITOS2 writes
    "##gff-version 3" and "#!gff-spec-version ...") are skipped. A
    row with the wrong column count, non-integer coordinates, or a
    gene row missing its Name/gene_id qualifier is a malformed file
    -- raised immediately with the offending line number, rather than
    silently skipped, since a corrupt annotation should never be
    mistaken for a clean one with fewer genes.

    Args:
        gff_path (str | Path): Path to MITOS2's result.gff.
        logger (Logger | None): Logger to use. Defaults to None.
        feature_type (str): The GFF "type" column value to keep. Defaults to "gene".

    Returns:
        list[GeneFeature]: The parsed features, in file order. An empty list is a valid, non-error outcome (zero genes found).

    Raises:
        FileNotFoundError: If gff_path does not exist.
        ValueError: If gff_path is not a file, or a kept row is malformed.
    """
    gff_path = Path(gff_path)

    if not gff_path.exists():
        if logger is not None: logger.error(f"File {gff_path} not found.")
        raise FileNotFoundError(f"File {gff_path} not found.")
    if not gff_path.is_file():
        if logger is not None: logger.error(f"File {gff_path} is not a file.")
        raise ValueError(f"File {gff_path} is not a file.")

    features: list[GeneFeature] = []
    with gff_path.open("r", encoding = "utf-8", errors = "replace") as handle:
        for line_number, raw_line in enumerate(handle, start = 1):
            line = raw_line.rstrip("\n")

            # Skipping blank lines and MITOS2's own header/comment lines.
            if not line.strip() or line.startswith("#"):
                continue

            columns = line.split("\t")
            if len(columns) != len(GFF_COLUMNS):
                message = f"MITOS2 GFF row {line_number} contains {len(columns)} columns, expected {len(GFF_COLUMNS)}."
                if logger is not None: logger.error(message)
                raise ValueError(message)

            fields = dict(zip(GFF_COLUMNS, columns, strict = True))

            # Filtering to the requested feature type before any
            # further normalization -- a malformed "region"/"exon"
            # row (MITOS2 always writes well-formed ones, but a
            # hand-edited file might not) should never fail parsing
            # of the "gene" rows this is actually meant to validate.
            if fields["type"] != feature_type:
                continue

            try:
                start = int(fields["start"])
                end = int(fields["end"])
            except ValueError as error:
                message = f"MITOS2 GFF row {line_number}: non-integer coordinates ({fields['start']!r}, {fields['end']!r})."
                if logger is not None: logger.error(message)
                raise ValueError(message) from error

            # GFF coordinates are 1-based and positive; a feature
            # spanning a circular genome's origin is real and valid
            # (MITOS2's own gffstr() writes it with end > the genome's
            # length rather than start > end, never rejected here) --
            # only a non-positive coordinate is actually malformed.
            if start <= 0 or end <= 0:
                message = f"MITOS2 GFF row {line_number}: non-positive coordinates ({start}, {end})."
                if logger is not None: logger.error(message)
                raise ValueError(message)

            strand = fields["strand"] if fields["strand"] in ("+", "-") else None
            score = None if fields["score"] == "." else float(fields["score"])
            qualifiers = _parse_attributes(fields["attributes"])

            gene_id = qualifiers.get("Name") or qualifiers.get("gene_id")
            if not gene_id:
                message = f"MITOS2 GFF row {line_number}: {feature_type} feature missing Name/gene_id qualifier."
                if logger is not None: logger.error(message)
                raise ValueError(message)

            features.append(GeneFeature(
                seqid = fields["seqid"],
                gene_id = gene_id,
                start = start,
                end = end,
                strand = strand,
                score = score,
                qualifiers = qualifiers,
            ))

    if logger is not None: logger.info(f"Parsed {len(features)} '{feature_type}' features from {gff_path}.")
    return features
