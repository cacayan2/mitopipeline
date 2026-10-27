"""
typst_style.py

Raw Typst blocks embedded directly in generated Markdown -- the same
mechanism v0.1.0's real, hand-written zmisc/docs/adr.md and
developer-manual.md use (pandoc passes a ```{=typst} fence straight
through to its own --pdf-engine typst, untouched) -- reproducing
their real page/font/paragraph preamble, title-page layout, and gray
callout box styling, confirmed directly against both documents'
source rather than guessed.
"""

# Imports
from datetime import datetime, timezone

def render_typst_preamble() -> str:
    """
    The page/font/paragraph/link preamble, matching
    developer-manual.md's own real settings exactly (0.5in margins,
    us-letter, New Computer Modern 11pt, justified paragraphs, blue
    links) -- the fuller of the two real preambles this project's
    documents use.

    Returns:
        str: The raw Typst preamble block, ready to prepend to a report's Markdown source.
    """
    return (
        "```{=typst}\n"
        "#set page(\n"
        '  paper: "us-letter",\n'
        "  margin: (x: 0.5in, y: 0.5in)\n"
        ")\n"
        "#set text(\n"
        '  font: "New Computer Modern",\n'
        "  size: 11pt\n"
        ")\n"
        "#set par(justify: true)\n"
        "#show link: set text(fill: blue)\n"
        "```\n"
    )

def render_title_block(title: str, subtitle: str, prepared_for: str = "") -> str:
    """
    A centered title page, matching developer-manual.md's own real
    stacked #text/#v spacing layout.

    Args:
        title (str): The report's main title (e.g. "Mitochondrial Genome Report: SAMPLE_01").
        subtitle (str): A subtitle line underneath the title.
        prepared_for (str): An optional "prepared for" line (e.g. a job id), omitted when blank.

    Returns:
        str: The raw Typst title block.
    """
    prepared_for_block = ""
    if prepared_for:
        prepared_for_block = (
            "#v(0.5cm)\n"
            "#text(size: 12pt, fill: gray)[\n"
            f"{prepared_for}\n"
            "]\n"
        )
    return (
        "```{=typst}\n"
        "#align(center)[\n\n"
        "#v(2.5cm)\n\n"
        '#text(size: 28pt, weight: "bold")[\n'
        f"{title}\n"
        "]\n\n"
        "#v(0.3cm)\n\n"
        "#text(size: 16pt)[\n"
        f"{subtitle}\n"
        "]\n\n"
        f"{prepared_for_block}"
        "#v(1.5cm)\n\n"
        "#text(size: 11pt, fill: gray)[\n"
        'Generated #datetime.today().display("[month repr:long] [day], [year]")\n'
        "]\n\n"
        "#v(2cm)\n\n"
        "]\n"
        "```\n"
    )

def render_callout_box_open(heading: str) -> str:
    """
    Opens a gray bordered callout box, matching
    developer-manual.md's own real "Key Takeaways"/"Purpose" box
    style -- rgb("#f2f2f2") fill, rgb("#cccccc") stroke, 4pt radius,
    centered. The box's own body is plain Markdown written between
    this and render_callout_box_close(), the same split-fence
    mechanism the real documents use (pandoc compiles the Markdown
    body to Typst markup in between, and the two raw fences supply
    only the surrounding #box(...)[ ... ] call).

    Args:
        heading (str): The box's own bold heading line.

    Returns:
        str: The raw Typst block that opens the box.
    """
    return (
        "```{=typst}\n"
        "#align(center)[\n"
        "#box(\n"
        "  inset: 12pt,\n"
        "  radius: 4pt,\n"
        '  fill: rgb("#f2f2f2"),\n'
        '  stroke: rgb("#cccccc"),\n'
        "[\n"
        f"#strong[{heading}]\n\n"
        "```\n"
    )

def render_callout_box_close() -> str:
    """Closes a callout box opened by render_callout_box_open()."""
    return "```{=typst}\n])\n]\n```\n"
