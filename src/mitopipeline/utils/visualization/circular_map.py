"""
circular_map.py

Normalizes MITOS2 gene features and renders a circular mitochondrial
genome map: splits a feature crossing the circular origin
-- represented by MITOS2 with an end coordinate beyond the linear
FASTA length, confirmed against MITOS2's own real gffstr() -- into
two drawable segments, then draws every protein-coding gene, tRNA,
and rRNA as an arc around a circular backbone.
"""

# Imports
from __future__ import annotations
import math
from dataclasses import dataclass, replace
from logging import Logger
from pathlib import Path
from Bio import SeqIO
import matplotlib.patches as patches
import matplotlib.pyplot as plt
from mitopipeline.utils.annotation.models import GeneFeature
from mitopipeline.utils.annotation.mitos_features import parse_mitos_gff

# The GFF "type" column values MITOS2 writes for a gene-level row of
# each biological category (confirmed against MITOS2's own real
# gffstr(): protein-coding genes are type "gene"; tRNA/rRNA genes are
# written twice -- once as "ncRNA_gene", once as their own "tRNA"/
# "rRNA" transcript row at the identical coordinates -- so reading
# only the latter two avoids double-drawing the same gene).
DRAWABLE_FEATURE_TYPES: tuple[str, ...] = ("gene", "tRNA", "rRNA")

FEATURE_COLORS: dict[str, str] = {
    "gene": "#4C78A8",
    "tRNA": "#54A24B",
    "rRNA": "#B279A2",
}

@dataclass(frozen = True, slots = True)
class DrawableFeature:
    """
    One gene feature tagged with the drawable category it was parsed
    under, since GeneFeature itself carries no such
    distinction -- it is purely a drawing-time concern.

    Attributes:
        feature (GeneFeature): The underlying parsed feature.
        feature_type (str): One of DRAWABLE_FEATURE_TYPES.
    """
    feature: GeneFeature
    feature_type: str

@dataclass(frozen = True, slots = True)
class CircularGenomeMapData:
    """
    Normalized data required to draw one sample's circular
    mitochondrial genome map.

    Attributes:
        sample_id (str): The sample id.
        genome_length (int): The linear FASTA record's own length, in base pairs.
        sequence_id (str): The FASTA record's own id.
        gc_content_percent (float): The genome's GC content, as a percentage.
        features (tuple[DrawableFeature, ...]): Every drawable feature, origin-crossing features already split, sorted by start coordinate.
    """
    sample_id: str
    genome_length: int
    sequence_id: str
    gc_content_percent: float
    features: tuple[DrawableFeature, ...]

def _read_genome_sequence(fasta_path: Path) -> tuple[str, int, float]:
    """
    Reads the genome FASTA's own sequence id, length, and GC content.
    A sample with more than one record (a fragmented assembly) uses
    the first, deterministically -- matching contig_selection's own
    "one representative sequence" choice for this sample.

    Args:
        fasta_path (Path): The genome FASTA to read.

    Returns:
        tuple[str, int, float]: The sequence id, its length, and its GC content percentage.

    Raises:
        FileNotFoundError: If fasta_path doesn't exist.
        ValueError: If fasta_path has no records or no bases.
    """
    fasta_path = Path(fasta_path)
    if not fasta_path.is_file():
        raise FileNotFoundError(f"Genome FASTA not found: {fasta_path}.")

    records = list(SeqIO.parse(fasta_path, "fasta"))
    if not records:
        raise ValueError(f"Genome FASTA contains no records: {fasta_path}.")

    record = records[0]
    sequence = str(record.seq).upper()
    if not sequence:
        raise ValueError(f"Genome FASTA record {record.id!r} has no bases: {fasta_path}.")

    gc_count = sequence.count("G") + sequence.count("C")
    gc_content_percent = gc_count / len(sequence) * 100.0

    return record.id, len(sequence), gc_content_percent

def split_origin_crossing_feature(feature: DrawableFeature, genome_length: int) -> list[DrawableFeature]:
    """
    Splits a feature crossing the circular genome's origin into two
    drawable segments. MITOS2 represents an origin-crossing feature
    with an end coordinate greater than the linear FASTA's own
    length (not start > end) -- confirmed against MITOS2's own real
    gffstr(), which adds the genome size to the wrapped stop
    coordinate before writing it.

    Args:
        feature (DrawableFeature): The feature to check/split.
        genome_length (int): The linear FASTA record's own length.

    Returns:
        list[DrawableFeature]: [feature] unchanged if it doesn't cross the origin, otherwise its two drawable segments.

    Raises:
        ValueError: If the feature's coordinates are invalid for this genome length.
    """
    gene = feature.feature
    if gene.start > genome_length:
        raise ValueError(f"Feature {gene.gene_id} starts beyond the genome length "
                         f"({genome_length}bp): {gene.start}-{gene.end}.")

    if gene.end <= genome_length:
        return [feature]

    wrapped_end = gene.end - genome_length
    if wrapped_end < 1 or wrapped_end > genome_length:
        raise ValueError(f"Feature {gene.gene_id} has invalid circular coordinates "
                         f"{gene.start}-{gene.end} for a {genome_length}bp genome.")

    ending_segment = replace(feature, feature = replace(gene, end = genome_length))
    beginning_segment = replace(feature, feature = replace(gene, start = 1, end = wrapped_end))
    return [beginning_segment, ending_segment]

def build_circular_map_data(sample_id: str, gff_path: str | Path, fasta_path: str | Path,
                             logger: Logger | None = None) -> CircularGenomeMapData:
    """
    Builds normalized data for one sample's circular mitochondrial
    genome map -- every real protein-coding gene, tRNA, and rRNA
    (the shared gene-feature parser, called once per drawable type), with any
    origin-crossing feature split into its drawable segments. A GFF
    with no drawable features at all (e.g. a failed annotation) is
    valid, empty input, not an error.

    Args:
        sample_id (str): The sample id.
        gff_path (str | Path): MITOS2's own result.gff.
        fasta_path (str | Path): The genome FASTA MITOS2 annotated.
        logger (Logger | None): Logger to use. Defaults to None.

    Returns:
        CircularGenomeMapData: The normalized map data.

    Raises:
        FileNotFoundError: If gff_path or fasta_path doesn't exist.
        ValueError: If the GFF/FASTA is malformed, or a feature has invalid coordinates.
    """
    sequence_id, genome_length, gc_content_percent = _read_genome_sequence(Path(fasta_path))

    normalized: list[DrawableFeature] = []
    for feature_type in DRAWABLE_FEATURE_TYPES:
        for gene in parse_mitos_gff(gff_path, logger = logger, feature_type = feature_type):
            segments = split_origin_crossing_feature(DrawableFeature(gene, feature_type), genome_length)
            if len(segments) > 1 and logger is not None:
                logger.info(f"{{{sample_id}}} Split origin-crossing feature {gene.gene_id} "
                           f"({gene.start}-{gene.end}) into 2 drawable segments.")
            normalized.extend(segments)

    normalized.sort(key = lambda d: (d.feature.start, d.feature.end, d.feature.gene_id))

    return CircularGenomeMapData(
        sample_id = sample_id, genome_length = genome_length, sequence_id = sequence_id,
        gc_content_percent = gc_content_percent, features = tuple(normalized),
    )

def _coordinate_to_degrees(coordinate: int, genome_length: int) -> float:
    """Converts a 1-based genome coordinate to a Matplotlib angle in degrees."""
    return (coordinate - 1) / genome_length * 360.0

def _feature_midpoint_angle(feature: GeneFeature, genome_length: int) -> float:
    """Returns a feature's midpoint angle, in radians."""
    midpoint = (feature.start + feature.end) / 2.0
    return midpoint / genome_length * 2.0 * math.pi

def _draw_feature_arc(axes, drawable: DrawableFeature, genome_length: int, radius: float, width: float) -> None:
    """Draws one gene feature as a circular arc."""
    gene = drawable.feature
    theta1 = _coordinate_to_degrees(gene.start, genome_length)
    theta2 = _coordinate_to_degrees(gene.end, genome_length)
    if theta2 <= theta1:
        theta2 = theta1 + 0.5

    color = FEATURE_COLORS.get(drawable.feature_type, "#777777")
    wedge = patches.Wedge(center = (0, 0), r = radius, theta1 = theta1, theta2 = theta2,
                          width = width, facecolor = color, edgecolor = "white",
                          linewidth = 0.6, alpha = 0.95)
    axes.add_patch(wedge)

def _draw_feature_label(axes, drawable: DrawableFeature, genome_length: int, radius: float) -> None:
    """Draws a gene's name at its midpoint angle, flipped upright on the circle's left half."""
    gene = drawable.feature
    angle = _feature_midpoint_angle(gene, genome_length)
    x, y = math.cos(angle) * radius, math.sin(angle) * radius
    angle_degrees = math.degrees(angle)

    rotation = angle_degrees
    horizontal_alignment = "left"
    if 90 < angle_degrees < 270:
        rotation = angle_degrees + 180
        horizontal_alignment = "right"

    axes.text(x, y, gene.gene_id, fontsize = 7, rotation = rotation, rotation_mode = "anchor",
              ha = horizontal_alignment, va = "center")

def _draw_tick_labels(axes, genome_length: int, radius: float, tick_count: int = 8) -> None:
    """Draws coordinate tick marks and labels evenly around the genome's circumference."""
    for index in range(tick_count):
        coordinate = int(genome_length * index / tick_count)
        label = "0 bp" if coordinate == 0 else f"{coordinate / 1000:.1f} kb"
        angle = index / tick_count * 2.0 * math.pi

        inner = (math.cos(angle) * (radius - 0.02), math.sin(angle) * (radius - 0.02))
        outer = (math.cos(angle) * (radius + 0.04), math.sin(angle) * (radius + 0.04))
        text = (math.cos(angle) * (radius + 0.13), math.sin(angle) * (radius + 0.13))

        axes.plot([inner[0], outer[0]], [inner[1], outer[1]], linewidth = 0.6, color = "black")
        axes.text(text[0], text[1], label, fontsize = 6, ha = "center", va = "center")

def _draw_legend(axes) -> None:
    """Draws the feature-type color legend."""
    handles = [patches.Patch(color = color, label = feature_type)
              for feature_type, color in FEATURE_COLORS.items()]
    axes.legend(handles = handles, loc = "lower center", bbox_to_anchor = (0.5, -0.08),
               ncol = len(handles), frameon = False, fontsize = 8)

def render_circular_map_data(map_data: CircularGenomeMapData, output_png: str | Path,
                             dpi: int = 600, logger: Logger | None = None) -> None:
    """
    Renders one sample's circular mitochondrial genome map to a PNG.
    An empty feature set still renders -- a bare, backbone-only
    circle is a valid outcome for a sample whose annotation found
    nothing, not an error.

    Args:
        map_data (CircularGenomeMapData): The normalized map data to draw.
        output_png (str | Path): Where to write the rendered PNG.
        dpi (int): The PNG's resolution. Defaults to 600.
        logger (Logger | None): Logger to use. Defaults to None.

    Raises:
        ValueError: If dpi or map_data.genome_length is not greater than zero.
    """
    if dpi <= 0:
        raise ValueError(f"Visualization DPI must be greater than zero, got {dpi}.")
    if map_data.genome_length <= 0:
        raise ValueError(f"Genome length must be greater than zero, got {map_data.genome_length}.")

    output_png = Path(output_png)
    output_png.parent.mkdir(parents = True, exist_ok = True)

    figure, axes = plt.subplots(figsize = (9, 9))
    try:
        axes.set_aspect("equal")
        axes.axis("off")

        outer_radius, forward_radius, reverse_radius, feature_width = 1.0, 1.0, 0.86, 0.10
        axes.add_patch(patches.Circle((0, 0), radius = outer_radius, fill = False,
                                      linewidth = 1.0, color = "black"))

        for drawable in map_data.features:
            radius = forward_radius if drawable.feature.strand != "-" else reverse_radius
            _draw_feature_arc(axes, drawable, map_data.genome_length, radius, feature_width)

        labeled_gene_ids: set[str] = set()
        for drawable in map_data.features:
            if drawable.feature.gene_id in labeled_gene_ids:
                continue
            label_radius = 1.18 if drawable.feature.strand != "-" else 0.72
            _draw_feature_label(axes, drawable, map_data.genome_length, label_radius)
            labeled_gene_ids.add(drawable.feature.gene_id)

        _draw_tick_labels(axes, map_data.genome_length, outer_radius)
        _draw_legend(axes)

        axes.text(0, 0.08, map_data.sample_id, fontsize = 14, fontweight = "bold",
                  ha = "center", va = "center")
        axes.text(0, -0.05, f"{map_data.genome_length:,} bp\nGC {map_data.gc_content_percent:.2f}%",
                  fontsize = 9, ha = "center", va = "center")

        axes.set_xlim(-1.45, 1.45)
        axes.set_ylim(-1.45, 1.45)
        figure.tight_layout()
        figure.savefig(output_png, dpi = dpi, bbox_inches = "tight")
    finally:
        plt.close(figure)

    if logger is not None:
        logger.info(f"{{{map_data.sample_id}}} Rendered circular mitochondrial genome map to {output_png}.")

def render_circular_genome_map(sample_id: str, gff_path: str | Path, fasta_path: str | Path,
                               output_png: str | Path, dpi: int = 600,
                               logger: Logger | None = None) -> CircularGenomeMapData:
    """
    Parses MITOS2's annotation output and renders one sample's
    circular mitochondrial genome map in a single call.

    Args:
        sample_id (str): The sample id.
        gff_path (str | Path): MITOS2's own result.gff.
        fasta_path (str | Path): The genome FASTA MITOS2 annotated.
        output_png (str | Path): Where to write the rendered PNG.
        dpi (int): The PNG's resolution. Defaults to 600.
        logger (Logger | None): Logger to use. Defaults to None.

    Returns:
        CircularGenomeMapData: The normalized map data that was rendered.
    """
    map_data = build_circular_map_data(sample_id, gff_path, fasta_path, logger = logger)
    render_circular_map_data(map_data, output_png, dpi = dpi, logger = logger)
    return map_data
