"""
tree_figure.py

Renders IQ-TREE's own real Newick output as a labeled PNG
figure, optionally midpoint-rooted -- a tree with no inherent root
(IQ-TREE's own unrooted maximum-likelihood output) is otherwise drawn
arbitrarily rooted at whichever tip happened to come first in the
file, which reads as a misleading claim about ancestry.

Each real tip carries a long, exact identifier (sample_dataset.py's
own "ScientificName__Accession"/"sample_id|assembled" convention) --
readable for a machine, cluttered on a plot with several tips. Tips
are relabeled to a short "name (N)" form here, numbered simply top to
bottom as drawn (not by tree position/relatedness), with a legend
below the plot mapping each number back to its own real accession/
identity -- confirmed directly that the original long labels made a
real 7-tip tree hard to read at a glance.
"""

# Imports
from __future__ import annotations
from logging import Logger
from pathlib import Path
from Bio import Phylo
import matplotlib.pyplot as plt

def _tip_label_and_detail(name: str) -> tuple[str, str]:
    """
    Splits one tip's real, exact identifier into a short display
    label and a fuller legend detail -- recognizes sample_dataset.py's
    own two naming conventions ("sample_id|assembled" for the query
    genome, "ScientificName__Accession" for a BLAST reference); any
    other tip name (e.g. a tree built outside this pipeline) is shown
    verbatim rather than guessed at.

    Args:
        name (str): The tip's real Newick/Clade name.

    Returns:
        tuple[str, str]: (short display label, full legend detail).
    """
    if name.endswith("|assembled"):
        sample_id = name[: -len("|assembled")]
        return sample_id, f"{sample_id} -- this sample's own assembled genome"
    if "__" in name:
        scientific_part, accession_part = name.rsplit("__", 1)
        scientific_name = scientific_part.replace("_", " ")
        return scientific_name, f"{scientific_name} -- {accession_part}"
    return name, name

def render_tree_figure(newick_path: str | Path, output_png: str | Path, dpi: int = 600,
                        midpoint_root: bool = True, logger: Logger | None = None) -> None:
    """
    Renders a Newick tree as a labeled PNG figure, with each tip
    shown as a short "name (N)" label and a legend below the plot
    mapping each number back to its own full identity (scientific
    name/accession, or "this sample's own assembled genome").

    Args:
        newick_path (str | Path): The real Newick tree file (IQ-TREE's own .treefile, or a stable copy of it).
        output_png (str | Path): Where to write the rendered PNG.
        dpi (int): The rendered PNG's resolution. Defaults to 600.
        midpoint_root (bool): Whether to midpoint-root the tree before drawing -- IQ-TREE's own output is unrooted, so the real tip furthest from the tree's own center is otherwise drawn first, misleadingly presented as a root. Defaults to True.
        logger (Logger | None): The logger to use for this step of the pipeline.

    Returns:
        None

    Raises:
        ValueError: If dpi is not greater than zero.
    """
    if dpi <= 0:
        raise ValueError(f"Tree figure DPI must be greater than zero, got {dpi}.")

    newick_path = Path(newick_path)
    output_png = Path(output_png)

    tree = Phylo.read(newick_path, "newick")
    if midpoint_root:
        tree.root_at_midpoint()

    # Numbered top to bottom in the order Bio.Phylo itself draws them
    # (get_terminals()'s own order, which Phylo.draw() preserves) --
    # not reordered by clade/relatedness, exactly a plain reading
    # position so the number carries no implied meaning of its own.
    terminals = tree.get_terminals()
    legend_lines: list[str] = []
    for index, clade in enumerate(terminals, start = 1):
        short_label, detail = _tip_label_and_detail(clade.name or "")
        clade.name = f"{short_label} ({index})"
        legend_lines.append(f"{index}. {detail}")

    tip_count = len(terminals)
    tree_height = max(4.0, tip_count * 0.3)
    legend_height = 0.22 * len(legend_lines) + 0.5
    figure, (tree_axis, legend_axis) = plt.subplots(
        2, 1, figsize = (10.0, tree_height + legend_height),
        gridspec_kw = {"height_ratios": [tree_height, legend_height]},
    )

    Phylo.draw(tree, axes = tree_axis, do_show = False, branch_labels = _branch_support_label)
    tree_axis.set_title(f"Maximum-likelihood phylogeny ({tip_count} tip(s))")

    legend_axis.axis("off")
    if legend_lines:
        legend_axis.text(
            0.0, 1.0, "Key\n" + "\n".join(legend_lines), fontsize = 8, family = "monospace",
            va = "top", ha = "left", transform = legend_axis.transAxes,
        )

    output_png.parent.mkdir(parents = True, exist_ok = True)
    figure.savefig(output_png, dpi = dpi, bbox_inches = "tight")
    plt.close(figure)

    if logger is not None:
        logger.info(f"Rendered tree figure ({tip_count} tip(s), midpoint_root={midpoint_root}) to {output_png}.")

def _branch_support_label(clade) -> str:
    """
    Labels an internal branch with its own confidence value (IQ-TREE's
    SH-aLRT/ultrafast bootstrap support, when present), leaving a tip
    branch and any branch with no recorded support unlabeled.

    Args:
        clade: The Bio.Phylo clade to label.

    Returns:
        str: The branch label, or an empty string when there is nothing real to show.
    """
    if clade.confidence is None:
        return ""
    return f"{clade.confidence:.0f}"
