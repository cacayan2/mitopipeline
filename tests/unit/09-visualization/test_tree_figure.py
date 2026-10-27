"""
test_tree_figure.py

Unit tests for rendering IQ-TREE's own real Newick output as a PNG
tree figure.
"""

# Imports
from pathlib import Path
import pytest
from mitopipeline.utils.visualization.tree_figure import render_tree_figure

REAL_NEWICK = "((A:0.1,B:0.2)90:0.05,(C:0.15,D:0.25)75:0.08,E:0.3);"

def test_render_tree_figure_writes_a_real_png(tmp_path: Path) -> None:
    newick_path = tmp_path / "tree.nwk"
    newick_path.write_text(REAL_NEWICK)
    output_png = tmp_path / "tree.png"

    render_tree_figure(newick_path, output_png)

    assert output_png.is_file()
    assert output_png.stat().st_size > 0

def test_render_tree_figure_rejects_a_non_positive_dpi(tmp_path: Path) -> None:
    newick_path = tmp_path / "tree.nwk"
    newick_path.write_text(REAL_NEWICK)

    with pytest.raises(ValueError):
        render_tree_figure(newick_path, tmp_path / "tree.png", dpi = 0)

def test_render_tree_figure_works_without_midpoint_rooting(tmp_path: Path) -> None:
    newick_path = tmp_path / "tree.nwk"
    newick_path.write_text(REAL_NEWICK)
    output_png = tmp_path / "tree.png"

    render_tree_figure(newick_path, output_png, midpoint_root = False)

    assert output_png.is_file()
    assert output_png.stat().st_size > 0

def test_render_tree_figure_creates_the_output_directory(tmp_path: Path) -> None:
    newick_path = tmp_path / "tree.nwk"
    newick_path.write_text(REAL_NEWICK)
    output_png = tmp_path / "nested" / "tree.png"

    render_tree_figure(newick_path, output_png)

    assert output_png.is_file()
