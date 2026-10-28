"""
pdf_export.py

Exports a Markdown report (already carrying its own embedded Typst
preamble/title-page/callout blocks, see typst_style.py) to PDF via
pandoc with Typst as the PDF engine (#361/#442) -- both confirmed
already installed on the real execution host.
"""

# Imports
import shutil
import subprocess
from dataclasses import dataclass
from logging import Logger
from pathlib import Path

@dataclass(frozen = True, slots = True)
class PdfExportResult:
    """
    The outcome of one Markdown-to-PDF export.

    Attributes:
        markdown_path (Path): The source Markdown report.
        pdf_path (Path): Where the PDF was written.
        success (bool): Whether the export produced a usable PDF.
        return_code (int): Pandoc's own process exit code.
        stderr (str): Pandoc's captured stderr, for diagnosing a failure.
    """
    markdown_path: Path
    pdf_path: Path
    success: bool
    return_code: int
    stderr: str

def _resolve_executable(name: str) -> str:
    """
    Resolves an executable name to a real path on $PATH.

    Args:
        name (str): The executable name to resolve.

    Returns:
        str: The resolved path.

    Raises:
        FileNotFoundError: If the executable isn't found.
    """
    path = shutil.which(name)
    if path is None:
        raise FileNotFoundError(f"Required executable not found on PATH: {name}.")
    return path

def export_markdown_to_pdf(markdown_path: str | Path, pdf_path: str | Path | None = None,
                            pandoc_bin: str = "pandoc", pdf_engine: str = "typst",
                            logger: Logger | None = None) -> PdfExportResult:
    """
    Exports one Markdown report to PDF. The PDF is written through a
    temporary file and only replaces the real destination once pandoc
    succeeds, so a failed export never leaves a partial/corrupt PDF
    behind.

    Args:
        markdown_path (str | Path): The Markdown report to export.
        pdf_path (str | Path | None): Where to write the PDF. Defaults to markdown_path with a .pdf suffix.
        pandoc_bin (str): The pandoc executable name. Defaults to "pandoc".
        pdf_engine (str): The PDF engine pandoc should use. Defaults to "typst".
        logger (Logger | None): Logger to use. Defaults to None.

    Returns:
        PdfExportResult: The export outcome.

    Raises:
        FileNotFoundError: If markdown_path, pandoc, or pdf_engine can't be found.
    """
    markdown_path = Path(markdown_path)
    if not markdown_path.is_file():
        raise FileNotFoundError(f"Markdown report not found: {markdown_path}.")

    pdf_path = Path(pdf_path) if pdf_path is not None else markdown_path.with_suffix(".pdf")
    pdf_path.parent.mkdir(parents = True, exist_ok = True)

    pandoc_path = _resolve_executable(pandoc_bin)
    _resolve_executable(pdf_engine)

    temporary_pdf_path = pdf_path.with_suffix(pdf_path.suffix + ".tmp")
    temporary_pdf_path.unlink(missing_ok = True)

    command = [
        pandoc_path, str(markdown_path), "--standalone",
        # gfm alone does not enable the raw_attribute extension, so a
        # ```{=typst} fence (the title page/callout box markup in
        # typst_style.py) would otherwise render as a literal code
        # block instead of being passed through to the PDF engine.
        # smart converts "--" into a real en dash -- confirmed
        # directly that gfm alone leaves a literal double hyphen in
        # the compiled PDF's own text layer (visually resembling a
        # dash via font kerning, but not one), even though gfm itself
        # already protects intraword underscores from being read as
        # emphasis without needing this extension.
        "--from", "gfm+raw_attribute+smart", "--to", "pdf", "--pdf-engine", pdf_engine,
        "--resource-path", str(markdown_path.parent),
        # Typst compiles through an intermediate file in its own temp
        # directory (pandoc's "toPdfViaTempFile"), so its default
        # sandbox root (the temp file's own directory) rejects any
        # absolute figure path outside it ("would escape the project
        # root", confirmed via a real compile). Widening the root to
        # the filesystem root lets Typst resolve real absolute image
        # paths anywhere on disk.
        "--pdf-engine-opt=--root=/",
        "--output", str(temporary_pdf_path),
    ]

    if logger is not None: logger.info(f"Running PDF export: {' '.join(command)}")
    result = subprocess.run(command, capture_output = True, text = True, check = False)

    success = result.returncode == 0 and temporary_pdf_path.is_file()
    if success:
        temporary_pdf_path.replace(pdf_path)
        if logger is not None: logger.info(f"Wrote PDF report to {pdf_path}.")
    else:
        temporary_pdf_path.unlink(missing_ok = True)
        if logger is not None: logger.error(f"PDF export failed for {markdown_path}: {result.stderr.strip()}")

    return PdfExportResult(
        markdown_path = markdown_path, pdf_path = pdf_path, success = success,
        return_code = result.returncode, stderr = result.stderr,
    )
