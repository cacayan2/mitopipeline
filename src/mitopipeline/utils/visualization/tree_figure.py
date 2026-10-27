"""
tree_figure.py

Renders IQ-TREE's own real Newick output (#436) as a labeled PNG
figure, optionally midpoint-rooted -- a tree with no inherent root
(IQ-TREE's own unrooted maximum-likelihood output) is otherwise drawn
arbitrarily rooted at whichever tip happened to come first in the
file, which reads as a misleading claim about ancestry.
"""

# Imports
from __future__ import annotations
from logging import Logger
from pathlib import Path
from Bio import Phylo
import matplotlib.pyplot as plt

def render_tree_figure(newick_path: str | Path, output_png: str | Path, dpi: int = 600,
                        midpoint_root: bool = True, logger: Logger | None = None) -> None:
    """
    Renders a Newick tree as a labeled PNG figure.

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

    tip_count = tree.count_terminals()
    figure_height = max(4.0, tip_count * 0.3)
    figure, axis = plt.subplots(figsize = (10.0, figure_height))

    Phylo.draw(tree, axes = axis, do_show = False, branch_labels = _branch_support_label)
    axis.set_title(f"Maximum-likelihood phylogeny ({tip_count} tip(s))")

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
