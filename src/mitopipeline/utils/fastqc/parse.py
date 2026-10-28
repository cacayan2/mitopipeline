"""
parse.py

Parses FastQC's own real output format directly out of its own zip
archive -- nothing in this project read FastQC's output before this
(only fastp's own JSON was ever parsed for the report), so raw/
trimmed/recruited-reads FastQC runs were always produced and then
silently never looked at again.
"""

# Imports
from __future__ import annotations
import zipfile
from dataclasses import dataclass
from logging import Logger
from pathlib import Path

@dataclass(frozen = True, slots = True)
class FastQCSummary:
    """
    One FastQC run's own real summary, read back from its zip output.

    Attributes:
        total_sequences (int | None): FastQC's own "Total Sequences" count.
        sequence_length (str | None): FastQC's own "Sequence length" value (a single number, or a range like "35-100").
        percent_gc (float | None): FastQC's own "%GC" value.
        module_verdicts (dict[str, str]): Every module's own PASS/WARN/FAIL verdict, keyed by module name.
        overall_verdict (str): The worst verdict across every module (FAIL > WARN > PASS) -- FastQC's own real per-module calls, never recomputed.
    """
    total_sequences: int | None
    sequence_length: str | None
    percent_gc: float | None
    module_verdicts: dict[str, str]
    overall_verdict: str

def _find_entry(names: list[str], suffix: str) -> str | None:
    """
    Finds the one zip entry ending in suffix. FastQC's own internal
    directory name is derived from the *input* FASTQ filename it was
    given, not the external zip filename a caller may have since
    renamed it to (confirmed directly: run_fastqc_raw.py's own
    postprocess_outputs() renames the zip file itself, but never
    touches what's inside it) -- so the exact internal path can't be
    assumed, only its own real suffix.

    Args:
        names (list[str]): A zip archive's own namelist().
        suffix (str): The entry's own expected suffix (e.g. "fastqc_data.txt").

    Returns:
        str | None: The matching entry name, or None if not found.
    """
    for name in names:
        if name.endswith(suffix):
            return name
    return None

def _parse_basic_statistics(data_text: str) -> tuple[int | None, str | None, float | None]:
    """Reads FastQC's own real "Basic Statistics" module out of a real fastqc_data.txt."""
    in_module = False
    total_sequences = sequence_length = percent_gc = None
    for line in data_text.splitlines():
        if line.startswith(">>Basic Statistics"):
            in_module = True
            continue
        if in_module and line.startswith(">>END_MODULE"):
            break
        if not in_module or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) != 2:
            continue
        measure, value = parts
        if measure == "Total Sequences":
            try: total_sequences = int(value)
            except ValueError: pass
        elif measure == "Sequence length":
            sequence_length = value
        elif measure == "%GC":
            try: percent_gc = float(value)
            except ValueError: pass
    return total_sequences, sequence_length, percent_gc

def _parse_module_verdicts(summary_text: str) -> dict[str, str]:
    """Reads FastQC's own real per-module PASS/WARN/FAIL calls out of a real summary.txt."""
    verdicts: dict[str, str] = {}
    for line in summary_text.splitlines():
        parts = line.split("\t")
        if len(parts) >= 2:
            verdicts[parts[1]] = parts[0]
    return verdicts

def parse_fastqc_zip(zip_path: str | Path, logger: Logger | None = None) -> FastQCSummary | None:
    """
    Reads one FastQC run's own real summary back from its zip output.

    Args:
        zip_path (str | Path): The FastQC zip to read.
        logger (Logger | None): Logger to use. Defaults to None.

    Returns:
        FastQCSummary | None: The parsed summary, or None if the zip is missing, unreadable, or missing its own expected contents.
    """
    zip_path = Path(zip_path)
    if not zip_path.is_file():
        return None

    try:
        with zipfile.ZipFile(zip_path) as archive:
            names = archive.namelist()
            data_entry = _find_entry(names, "fastqc_data.txt")
            summary_entry = _find_entry(names, "summary.txt")
            if data_entry is None or summary_entry is None:
                if logger is not None: logger.warning(f"FastQC zip {zip_path} is missing fastqc_data.txt/summary.txt.")
                return None
            total_sequences, sequence_length, percent_gc = _parse_basic_statistics(
                archive.read(data_entry).decode("utf-8", errors = "replace")
            )
            module_verdicts = _parse_module_verdicts(archive.read(summary_entry).decode("utf-8", errors = "replace"))
    except (zipfile.BadZipFile, OSError) as e:
        if logger is not None: logger.warning(f"Could not read FastQC zip {zip_path}: {e}")
        return None

    verdict_rank = {"PASS": 0, "WARN": 1, "FAIL": 2}
    overall_verdict = "PASS"
    for verdict in module_verdicts.values():
        if verdict_rank.get(verdict, 0) > verdict_rank.get(overall_verdict, 0):
            overall_verdict = verdict

    return FastQCSummary(
        total_sequences = total_sequences, sequence_length = sequence_length, percent_gc = percent_gc,
        module_verdicts = module_verdicts, overall_verdict = overall_verdict,
    )

def extract_fastqc_image(zip_path: str | Path, image_name: str, output_path: str | Path,
                          logger: Logger | None = None) -> Path | None:
    """
    Extracts one named image (e.g. "per_base_quality.png") from a real
    FastQC zip's own Images/ folder to a real file on disk, for
    embedding as a report figure -- FastQC's own real chart, never
    regenerated or approximated.

    Args:
        zip_path (str | Path): The FastQC zip to read.
        image_name (str): The image's own filename inside FastQC's Images/ folder (e.g. "per_base_quality.png").
        output_path (str | Path): Where to write the extracted image.
        logger (Logger | None): Logger to use. Defaults to None.

    Returns:
        Path | None: output_path, or None if the image couldn't be extracted.
    """
    zip_path = Path(zip_path)
    output_path = Path(output_path)
    if not zip_path.is_file():
        return None

    try:
        with zipfile.ZipFile(zip_path) as archive:
            entry = _find_entry(archive.namelist(), f"Images/{image_name}")
            if entry is None:
                if logger is not None: logger.warning(f"FastQC zip {zip_path} has no Images/{image_name}.")
                return None
            output_path.parent.mkdir(parents = True, exist_ok = True)
            output_path.write_bytes(archive.read(entry))
    except (zipfile.BadZipFile, OSError) as e:
        if logger is not None: logger.warning(f"Could not extract {image_name} from {zip_path}: {e}")
        return None
    return output_path
