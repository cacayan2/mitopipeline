"""
package.py

Packages the real SRA/BioSample/GenBank archival tables plus the
validation report into one submission-ready directory (#375),
matching the layout the project's own NCBI archival requirements
research proposed -- reproducible MD5 checksums over every real file
collected, and a human-readable submission_checklist.md summarizing
exactly what still blocks each sample, never a silent "done."
"""

# Imports
import csv
import hashlib
from logging import Logger
from pathlib import Path

def _file_checksum(path: Path) -> str:
    """Computes a real MD5 checksum of a file's own contents, streamed to handle large FASTAs."""
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()

def build_checksum_manifest(file_paths: list[str | Path], base_directory: str | Path) -> list[dict[str, str]]:
    """
    Computes real MD5 checksums for every real file in the archival
    package, relative to the package's own base directory.

    Args:
        file_paths (list[str | Path]): The real files to checksum (only existing files are included).
        base_directory (str | Path): The package's own base directory, used to record relative paths.

    Returns:
        list[dict[str, str]]: One {"path", "md5"} record per existing file, sorted by path.
    """
    base_directory = Path(base_directory)
    records = []
    for file_path in file_paths:
        file_path = Path(file_path)
        if not file_path.is_file():
            continue
        try:
            relative_path = file_path.relative_to(base_directory)
        except ValueError:
            relative_path = file_path
        records.append({"path": str(relative_path), "md5": _file_checksum(file_path)})
    return sorted(records, key = lambda r: r["path"])

def write_checksum_manifest(records: list[dict[str, str]], output_path: str | Path, logger: Logger | None = None) -> None:
    """Writes the checksum manifest as a TSV."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents = True, exist_ok = True)
    with output_path.open("w", encoding = "utf-8", newline = "") as handle:
        writer = csv.DictWriter(handle, fieldnames = ["path", "md5"], delimiter = "\t")
        writer.writeheader()
        writer.writerows(records)
    if logger is not None: logger.info(f"Wrote checksum manifest ({len(records)} file(s)) to {output_path}.")

def render_submission_checklist(validation_records: list[dict[str, str]], job_id: str) -> str:
    """
    Renders a human-readable submission checklist -- the real
    "Submission Blockers" report the project's own NCBI archival
    requirements research called for, so a PI sees exactly what
    still needs resolving per sample before any upload.

    Args:
        validation_records (list[dict[str, str]]): The records built by validate_archival_readiness()/build_archival_validation().
        job_id (str): The pipeline job id.

    Returns:
        str: The checklist, as Markdown text.
    """
    ready_count = sum(1 for r in validation_records if r["ready_for_review"] == "True")
    lines = [
        f"# NCBI Submission Checklist: {job_id}", "",
        f"{ready_count} of {len(validation_records)} sample(s) ready for submission review.", "",
        (
            "This pipeline never submits to NCBI directly. These files are prepared for "
            "human review and upload through the NCBI Submission Portal."
        ), "",
    ]

    for record in validation_records:
        lines.append(f"## {record['sample_id']} ({record['organism'] or 'organism not set'})")
        lines.append("")
        checks = [
            ("Raw FASTQ files", record["raw_fastq_ready"]),
            ("Assembly", record["assembly_ready"]),
            ("Annotation", record["annotation_ready"]),
            ("SRA metadata", record["sra_ready"]),
            ("BioSample metadata", record["biosample_ready"]),
            ("GenBank submission folder", record["genbank_ready"]),
        ]
        for label, ready in checks:
            lines.append(f"- [{'x' if ready == 'True' else ' '}] {label}")
        if record["submission_blockers"]:
            lines.append("")
            lines.append(f"**Blockers:** {record['submission_blockers']}")
        lines.append("")

    return "\n".join(lines)

def write_submission_checklist(validation_records: list[dict[str, str]], job_id: str, output_path: str | Path,
                                logger: Logger | None = None) -> Path:
    """
    Writes the submission checklist to disk.

    Args:
        validation_records (list[dict[str, str]]): The records built by build_archival_validation().
        job_id (str): The pipeline job id.
        output_path (str | Path): Where to write the checklist.
        logger (Logger | None): Logger to use. Defaults to None.

    Returns:
        Path: The written checklist's path.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents = True, exist_ok = True)
    output_path.write_text(render_submission_checklist(validation_records, job_id), encoding = "utf-8")
    if logger is not None: logger.info(f"Wrote submission checklist to {output_path}.")
    return output_path
