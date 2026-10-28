"""
typst_style.py

Raw Typst blocks embedded directly in generated Markdown -- the same
mechanism v0.1.0's real, hand-written zmisc/docs/adr.md and
developer-manual.md use (pandoc passes a ```{=typst} fence straight
through to its own --pdf-engine typst, untouched) -- reproducing
their real page/font/paragraph preamble, title-page layout, and gray
callout box styling, confirmed directly against both documents'
source rather than guessed, extended with table/figure/badge/footer
styling for the report pipeline.
"""

# Imports

# A fixed, colorblind-safe status palette (Okabe-Ito-adjacent hues,
# each tier also distinct in LIGHTNESS so grayscale printing keeps
# them apart -- the badge's own text label always carries the actual
# meaning, color is a supplementary signal, never the only one).
# (fill, text) as hex strings, keyed by the raw status value.
STATUS_BADGE_STYLES: dict[str, tuple[str, str]] = {
    "PASS": ("#d4edda", "#155724"),
    "READY": ("#d4edda", "#155724"),
    "REVIEW": ("#fff3cd", "#856404"),
    "MANUAL_REVIEW_REQUIRED": ("#fff3cd", "#856404"),
    "WARN": ("#fff3cd", "#856404"),
    "FAIL": ("#f8d7da", "#721c24"),
    "NOT ASSESSED": ("#e2e3e5", "#383d41"),
}

def render_typst_preamble(sample_id: str | None = None) -> str:
    """
    The page/font/paragraph/link/table preamble, matching
    developer-manual.md's own real settings (0.5in margins,
    us-letter, New Computer Modern 11pt, justified body paragraphs,
    blue links), extended with: a page footer, consistent heading
    sizing, and table cells that never overflow the page (pandoc's
    own Typst tables use bare auto-sized columns with no wrapping,
    confirmed directly to overlap when a cell's content is long --
    real wrapping is applied at the Python layer via wrap_cell_text,
    this just keeps header rows non-orphaned and cell text readable).

    Args:
        sample_id (str | None): The sample id to show in the page footer, or None for a run-level report (job id is passed by the caller as part of the footer text instead).

    Returns:
        str: The raw Typst preamble block, ready to prepend to a report's Markdown source.
    """
    footer_prefix = f"{sample_id} · " if sample_id else ""
    return (
        "```{=typst}\n"
        "#set page(\n"
        '  paper: "us-letter",\n'
        "  margin: (x: 0.5in, y: 0.5in),\n"
        "  footer: context {\n"
        '    set text(size: 8.5pt, fill: rgb("#888888"))\n'
        "    line(length: 100%, stroke: 0.4pt + rgb(\"#cccccc\"))\n"
        "    v(4pt)\n"
        f"    [{footer_prefix}page #counter(page).display() of #counter(page).final().first()]\n"
        "  }\n"
        ")\n"
        "#set text(\n"
        '  font: "New Computer Modern",\n'
        "  size: 11pt\n"
        ")\n"
        "#set par(justify: true, leading: 0.6em)\n"
        "#show link: set text(fill: blue)\n"
        "#set figure.caption(separator: [. ])\n"
        "#show figure: set block(breakable: false, above: 0.8em, below: 0.8em)\n"
        "#show heading.where(level: 1): set text(size: 20pt)\n"
        "#show heading.where(level: 2): set text(size: 14pt)\n"
        "#show heading.where(level: 2): set block(above: 0.9em, below: 0.5em)\n"
        "#show heading.where(level: 2): it => block(breakable: false, sticky: true, it)\n"
        "#show heading.where(level: 3): set text(size: 12pt)\n"
        # A genuine ruled/shaded table, not pandoc's own bare,
        # borderless default -- a light full grid, a shaded header
        # row, and subtle zebra striping on the body rows, applied
        # once here as table-wide defaults so every pandoc-generated
        # #table(...) call picks it up without needing its own
        # per-call styling (a #set rule only changes default
        # constructor arguments -- unlike the #show table: it => ...
        # rebuild this project already tried and abandoned earlier
        # this session for recursing/blanking the page, this never
        # touches table's own content, only its defaults).
        "#set table(\n"
        '  stroke: 0.5pt + rgb("#bbbbbb"),\n'
        "  inset: 6pt,\n"
        "  fill: (x, y) => if y == 0 { "
        'rgb("#e8e8e8") } else if calc.odd(y) { rgb("#f6f6f6") } else { white },\n'
        ")\n"
        "#show table.cell.where(y: 0): set text(weight: \"bold\")\n"
        "```\n"
    )

def render_title_block(title: str, species: str, sample_id: str, generated_at: str) -> str:
    """
    A centered title page, matching developer-manual.md's own real
    stacked #text/#v spacing layout -- each of title/species/sample
    id/timestamp shown exactly once (previously duplicated by a
    second, separately-written set of bold lines immediately
    following this block).

    Args:
        title (str): The report's main title (e.g. "Mitochondrial Genome Report").
        species (str): The sample's scientific name, rendered in italics.
        sample_id (str): The sample id.
        generated_at (str): A pre-formatted "DD Mon YYYY, HH:MM UTC" timestamp -- computed in Python, not Typst's own compile-time date, since the report's real generation time (not just date) matters.

    Returns:
        str: The raw Typst title block.
    """
    return (
        "```{=typst}\n"
        "#align(center)[\n\n"
        "#v(1.0cm)\n\n"
        '#text(size: 26pt, weight: "bold")[\n'
        f"{title}\n"
        "]\n\n"
        "#v(0.15cm)\n\n"
        f"#text(size: 16pt)[#emph[{species}]]\n\n"
        "#v(0.3cm)\n\n"
        "#text(size: 12pt, fill: gray)[\n"
        f"Sample ID: {sample_id}\n"
        "]\n\n"
        "#v(0.15cm)\n\n"
        "#text(size: 11pt, fill: gray)[\n"
        f"Generated {generated_at}\n"
        "]\n\n"
        "#v(0.8cm)\n\n"
        "]\n"
        "```\n"
    )

def render_run_title_block(title: str, job_id: str, sample_count_label: str, generated_at: str) -> str:
    """The run-level report's own title page, mirroring render_title_block()'s single-instance-per-field layout."""
    return (
        "```{=typst}\n"
        "#align(center)[\n\n"
        "#v(1.0cm)\n\n"
        '#text(size: 26pt, weight: "bold")[\n'
        f"{title}\n"
        "]\n\n"
        "#v(0.2cm)\n\n"
        f"#text(size: 16pt)[{sample_count_label}]\n\n"
        "#v(0.3cm)\n\n"
        "#text(size: 12pt, fill: gray)[\n"
        f"Job ID: {job_id}\n"
        "]\n\n"
        "#v(0.15cm)\n\n"
        "#text(size: 11pt, fill: gray)[\n"
        f"Generated {generated_at}\n"
        "]\n\n"
        "#v(0.8cm)\n\n"
        "]\n"
        "```\n"
    )

def render_callout_box_open() -> str:
    """
    Opens a gray bordered callout box, matching
    developer-manual.md's own real "Key Takeaways"/"Purpose" box
    style -- rgb("#f2f2f2") fill, rgb("#cccccc") stroke, 4pt radius,
    centered. The box's own body is plain Markdown written between
    this and render_callout_box_close(), the same split-fence
    mechanism the real documents use. Deliberately does not print its
    own heading -- the enclosing Markdown section's own "##" heading
    is the only heading shown, previously duplicated by this box
    printing a second, bold one inside itself.

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
        "  width: 100%,\n"
        "[\n"
        "```\n"
    )

def render_callout_box_close() -> str:
    """Closes a callout box opened by render_callout_box_open()."""
    return "```{=typst}\n])\n]\n```\n"

def render_figure(image_path: str, caption: str) -> str:
    """
    Wraps an image in a proper numbered Typst figure (never altering
    the image itself): scaled to the available text width, kept on
    the same page as whatever heading precedes it. Pandoc's own plain
    Markdown image syntax produces a bare, uncaptioned, unscaled
    image that Typst is free to push onto its own page, leaving the
    heading on the page before it mostly blank -- confirmed directly.

    Args:
        image_path (str): The image path, relative to the report's own directory (matching markdown_report.py's existing relative_path()).
        caption (str): The figure's full caption text, already composed (numbering, dimensions, any circularity caveat).

    Returns:
        str: The raw Typst figure block.
    """
    escaped_caption = caption.replace('"', '\\"')
    return (
        "```{=typst}\n"
        "#figure(\n"
        f'  image("{image_path}", width: 100%),\n'
        f'  caption: [{escaped_caption}],\n'
        "  placement: none,\n"
        ")\n"
        "```\n"
    )
