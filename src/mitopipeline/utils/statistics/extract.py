"""
extract.py

Extracts per-sample statistics (#353-356) from the real result JSON
each earlier stage already wrote -- assembly from contig_selection's
own selection record, BLAST from the ranked top-hit TSV plus the
biological validation gate's own verdict, MITOS2 annotation from its
own result JSON -- then merges them by sample_id into one record.
"""

# Imports
import json
from logging import Logger
from pathlib import Path
from Bio import SeqIO
from mitopipeline.utils.blast.blast_hits import read_top_blast_hits, resolve_scientific_name_from_fasta
from mitopipeline.utils.statistics.models import (
    AnnotationStatistics,
    AssemblyStatistics,
    BlastStatistics,
    SampleStatistics,
)
from mitopipeline.utils.validation.gates import CORE_ANIMAL_MT_PROTEIN_GENES

def _read_json(path: Path) -> dict | None:
    """Reads a JSON file, returning None when it's absent or unreadable -- a missing earlier stage, not an error here."""
    if not path.is_file():
        return None
    try:
        return json.loads(path.read_text(encoding = "utf-8"))
    except (OSError, json.JSONDecodeError):
        return None

def extract_assembly_statistics(sample_id: str, contig_selection_json: str | Path,
                                 recruited_reads_qc_json: str | Path | None = None,
                                 runtime_seconds: float | None = None,
                                 logger: Logger | None = None) -> AssemblyStatistics | None:
    """
    Extracts assembly statistics from contig_selection's own real
    selection record -- length/status/source_attempt directly, and
    contig_count/GC content computed from the selected FASTA itself.

    Args:
        sample_id (str): The sample id.
        contig_selection_json (str | Path): The selection JSON written by run_contig_selection.py.
        recruited_reads_qc_json (str | Path | None): GetOrganelle's own recruited-reads QC record for the selected attempt, if available.
        runtime_seconds (float | None): The selected attempt's own assembly runtime, if available.
        logger (Logger | None): Logger to use. Defaults to None.

    Returns:
        AssemblyStatistics | None: None if contig_selection found nothing usable, or its record is missing/unreadable.
    """
    selection = _read_json(Path(contig_selection_json))
    if selection is None or not selection.get("selected"):
        if logger is not None: logger.info(f"{{{sample_id}}} No selected assembly; skipping assembly statistics.")
        return None

    fasta_path = Path(selection["fasta_path"])
    records = list(SeqIO.parse(fasta_path, "fasta"))
    total_sequence = "".join(str(r.seq) for r in records).upper()
    gc_count = total_sequence.count("G") + total_sequence.count("C")
    gc_content_percent = (gc_count / len(total_sequence) * 100.0) if total_sequence else 0.0

    recruited_read_count = None
    qc_record = _read_json(Path(recruited_reads_qc_json)) if recruited_reads_qc_json else None
    if qc_record is not None:
        r1, r2 = qc_record.get("r1_read_count"), qc_record.get("r2_read_count")
        if r1 is not None and r2 is not None:
            recruited_read_count = r1 + r2

    return AssemblyStatistics(
        sample_id = sample_id, status = selection["status"], source_attempt = selection["source_attempt"],
        length_bp = selection["length_bp"], contig_count = len(records), gc_content_percent = gc_content_percent,
        recruited_read_count = recruited_read_count, runtime_seconds = runtime_seconds,
    )

def extract_blast_statistics(sample_id: str, top_hits_tsv: str | Path,
                              biological_validation_json: str | Path | None = None,
                              circular_coverage_json: str | Path | None = None,
                              reference_fasta: str | Path | None = None,
                              logger: Logger | None = None) -> BlastStatistics | None:
    """
    Extracts BLAST statistics from the ranked top-hit TSV's own
    highest-ranked hit, plus the biological validation gate's own
    "blast_identity" verdict for whichever attempt was selected.

    Args:
        sample_id (str): The sample id.
        top_hits_tsv (str | Path): The ranked top-hit TSV written by run_blast_hit_selection.py.
        biological_validation_json (str | Path | None): The validation result JSON written by run_biological_validation.py, if available.
        circular_coverage_json (str | Path | None): The result JSON written by run_circular_coverage.py, for a circular selected attempt only -- when present and successful, this replaces the raw top-hit query coverage in the reported statistics too, so the report never shows the same misleading number the validation gate itself no longer uses.
        reference_fasta (str | Path | None): The local BLAST reference FASTA, used to resolve a scientific name from the hit's own real header description when BLAST's own sscinames column came back empty/unresolved.
        logger (Logger | None): Logger to use. Defaults to None.

    Returns:
        BlastStatistics | None: None if the top-hit TSV is missing or empty (BLAST found nothing, or didn't run).
    """
    top_hits_tsv = Path(top_hits_tsv)
    if not top_hits_tsv.is_file():
        if logger is not None: logger.info(f"{{{sample_id}}} No top-hit TSV; skipping BLAST statistics.")
        return None

    hits = read_top_blast_hits(top_hits_tsv, sample_id = sample_id, logger = logger)
    if not hits:
        if logger is not None: logger.info(f"{{{sample_id}}} No BLAST hits found; skipping BLAST statistics.")
        return None

    top_hit = hits[0]

    decision_status = None
    validation = _read_json(Path(biological_validation_json)) if biological_validation_json else None
    if validation is not None:
        selected_attempt = validation.get("selected_attempt")
        for result in validation.get("results", []):
            if result.get("attempt") == selected_attempt:
                for check in result.get("checks", []):
                    if check.get("name") == "blast_identity":
                        decision_status = check.get("verdict")
                break

    top_query_coverage_percent = top_hit.qcovs
    circular_coverage = _read_json(Path(circular_coverage_json)) if circular_coverage_json else None
    if circular_coverage is not None and circular_coverage.get("success") and circular_coverage.get("corrected_query_coverage") is not None:
        top_query_coverage_percent = circular_coverage["corrected_query_coverage"]

    scientific_name = top_hit.sscinames
    if (not scientific_name or scientific_name.strip().upper() in {"N/A", ""}) and reference_fasta is not None:
        resolved = resolve_scientific_name_from_fasta(reference_fasta, top_hit.sseqid, logger = logger)
        if resolved:
            scientific_name = resolved
            if logger is not None: logger.info(f"{{{sample_id}}} Resolved BLAST scientific name from the local reference FASTA: {resolved}.")

    return BlastStatistics(
        sample_id = sample_id, top_accession = top_hit.sseqid, top_scientific_name = scientific_name,
        top_identity_percent = top_hit.pident, top_query_coverage_percent = top_query_coverage_percent,
        top_evalue = top_hit.evalue, decision_status = decision_status,
    )

def extract_annotation_statistics(sample_id: str, annotation_result_json: str | Path,
                                   logger: Logger | None = None) -> AnnotationStatistics | None:
    """
    Extracts annotation statistics from MITOS2's own result JSON --
    total gene count directly, core-gene coverage/gaps computed
    against CORE_ANIMAL_MT_PROTEIN_GENES (#414's own reference set).

    Args:
        sample_id (str): The sample id.
        annotation_result_json (str | Path): The result JSON written by run_mitos2_annotation.py.
        logger (Logger | None): Logger to use. Defaults to None.

    Returns:
        AnnotationStatistics | None: None if MITOS2 didn't run successfully, or its record is missing/unreadable.
    """
    result = _read_json(Path(annotation_result_json))
    if result is None or not result.get("success"):
        if logger is not None: logger.info(f"{{{sample_id}}} No successful annotation; skipping annotation statistics.")
        return None

    gene_names = set(result.get("gene_names") or [])
    found_core = gene_names & CORE_ANIMAL_MT_PROTEIN_GENES
    missing_core = tuple(sorted(CORE_ANIMAL_MT_PROTEIN_GENES - gene_names))

    return AnnotationStatistics(
        sample_id = sample_id, gene_count = result.get("feature_count") or 0,
        core_gene_count = len(found_core), missing_core_genes = missing_core,
        trna_count = result.get("trna_count"), rrna_count = result.get("rrna_count"),
    )

def merge_sample_statistics(sample_id: str, assembly: AssemblyStatistics | None,
                             blast: BlastStatistics | None,
                             annotation: AnnotationStatistics | None) -> SampleStatistics:
    """
    Joins one sample's assembly/BLAST/annotation statistics into a
    single record (#356). Any stage that didn't produce a usable
    result stays None -- never backfilled with a fabricated value.

    Args:
        sample_id (str): The sample id.
        assembly (AssemblyStatistics | None): This sample's assembly statistics, if available.
        blast (BlastStatistics | None): This sample's BLAST statistics, if available.
        annotation (AnnotationStatistics | None): This sample's annotation statistics, if available.

    Returns:
        SampleStatistics: The merged record.
    """
    return SampleStatistics(sample_id = sample_id, assembly = assembly, blast = blast, annotation = annotation)
