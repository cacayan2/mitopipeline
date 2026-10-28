"""
test_tree_figure.py

Unit tests for rendering IQ-TREE's own real Newick output as a PNG
tree figure.
"""

# Imports
from pathlib import Path
import pytest
from mitopipeline.utils.visualization.tree_figure import _tip_label_and_detail, render_tree_figure

REAL_NEWICK = "((A:0.1,B:0.2)90:0.05,(C:0.15,D:0.25)75:0.08,E:0.3);"
SAMPLE_DATASET_NEWICK = (
    "((gargammel_circular|assembled:0.01,Danio_rerio__NC_002333.2:0.02)100:0.05,"
    "Xenocypris_yunnanensis__NC_035954.1:0.08,unknown_species__NC_015838.1:0.1);"
)

def test_tip_label_and_detail_for_the_assembled_genome() -> None:
    label, detail = _tip_label_and_detail("gargammel_circular|assembled")
    assert label == "gargammel_circular"
    assert "this sample's own assembled genome" in detail

def test_tip_label_and_detail_for_a_blast_reference() -> None:
    label, detail = _tip_label_and_detail("Danio_rerio__NC_002333.2")
    assert label == "Danio rerio"
    assert detail == "Danio rerio -- NC_002333.2"

def test_tip_label_and_detail_falls_back_to_the_raw_name() -> None:
    """A tip name that matches neither naming convention (e.g. a tree
    built outside this pipeline) is shown verbatim, never guessed at."""
    assert _tip_label_and_detail("some_other_tip") == ("some_other_tip", "some_other_tip")

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

def test_render_tree_figure_handles_real_sample_dataset_tip_names(tmp_path: Path) -> None:
    """Regression: the real tip names sample_dataset.py produces
    (long, exact, with a `|` and `__` in them) must render without
    crashing and with the new numbered-label/legend layout."""
    newick_path = tmp_path / "tree.nwk"
    newick_path.write_text(SAMPLE_DATASET_NEWICK)
    output_png = tmp_path / "tree.png"

    render_tree_figure(newick_path, output_png)

    assert output_png.is_file()
    assert output_png.stat().st_size > 0

def test_render_tree_figure_creates_the_output_directory(tmp_path: Path) -> None:
    newick_path = tmp_path / "tree.nwk"
    newick_path.write_text(REAL_NEWICK)
    output_png = tmp_path / "nested" / "tree.png"

    render_tree_figure(newick_path, output_png)

    assert output_png.is_file()
