"""
run_biological_validation.py

Execution layer for biological validation gates: runs the same
PASS/REVIEW/FAIL checks (gates.validate_assembly) against every
assembly attempt considered for a sample -- the primary attempt and
every rescue attempt, identically -- and records which one contig selection's
contig_selection actually chose. BLAST and MITOS2 evidence only ever
apply to the selected attempt, since those downstream stages only
ever run once, on whichever candidate contig_selection picked.

Never equates a clean tool exit with biological validity -- every
claim here traces back to one of gates.py's own checks.
"""

# Imports
import argparse
import json
from pathlib import Path
from mitopipeline.logging.logger_factory import make_logger
from mitopipeline.utils.contig_selection.models import ContigCandidate
from mitopipeline.utils.validation.gates import validate_assembly

def parse_args() -> argparse.Namespace:
    """
    Parses the command line arguments for this script.

    Returns:
        argparse.Namespace: The parsed command line arguments.
    """
    parser = argparse.ArgumentParser(description = "Runs biological validation gates against every assembly attempt considered for one sample.")

    parser.add_argument("--sample-id", type = str, required = True, help = "The sample id.")
    parser.add_argument("--attempt-dir", type = Path, required = True, action = "append",
                        help = "An assembly attempt's output directory, same order as passed to run_contig_selection.py (first = primary).")
    parser.add_argument("--contig-selection-json", type = Path, required = True,
                        help = "The selection JSON written by run_contig_selection.py -- tells this step which attempt was selected.")
    parser.add_argument("--min-length", type = int, required = True, help = "The configured minimum complete-genome length (taxon-aware or default).")
    parser.add_argument("--max-length", type = int, required = True, help = "The configured maximum complete-genome length (taxon-aware or default).")
    parser.add_argument("--seed-file", type = str, default = None, help = "The seed file used for the selected attempt, if any (provenance only).")
    parser.add_argument("--blast-result-json", type = Path, default = None, help = "The result JSON written by run_blast.py, for the selected attempt only.")
    parser.add_argument("--blast-top-hits-tsv", type = Path, default = None, help = "The ranked top-hit TSV (write_top_blast_hits' own output), for the conflicting-taxa check.")
    parser.add_argument("--annotation-result-json", type = Path, default = None, help = "The result JSON written by run_mitos2_annotation.py, for the selected attempt only.")
    parser.add_argument("--circular-coverage-json", type = Path, default = None,
                        help = "The result JSON written by run_circular_coverage.py, for a circular selected attempt only -- when present and successful, this overrides the raw top-hit query coverage (linear BLAST can't see across a circular genome's own origin).")
    parser.add_argument("--output-dir", type = Path, required = True, help = "The directory to write this step's own output to.")
    parser.add_argument("--log-file", type = Path, required = True, help = "The path to the log file.")
    parser.add_argument("--global-log-file", type = Path, required = True, help = "The path to the global log file.")

    return parser.parse_args()

def _load_candidate(attempt: int, sample_id: str, attempt_dir: Path) -> ContigCandidate:
    """
    Reads one assembly attempt's recorded status and standardized
    output path -- the same contract run_contig_selection.py's own
    _load_candidate() reads, duplicated here deliberately so this
    step never depends on contig_selection's internals, only the
    documented file contract both read.

    Args:
        attempt (int): The attempt number.
        sample_id (str): The sample id.
        attempt_dir (Path): The attempt's own output directory.

    Returns:
        ContigCandidate: The candidate read from this attempt.
    """
    status_path = attempt_dir / f"{sample_id}.assembly_status.json"
    try:
        status = json.loads(status_path.read_text(encoding = "utf-8"))["status"]
    except (OSError, json.JSONDecodeError, KeyError):
        return ContigCandidate(attempt = attempt, status = "missing", fasta_path = None, gfa_path = None)

    fasta_path = attempt_dir / "data.fasta"
    if status not in ("circular", "fragmented") or not fasta_path.is_file():
        return ContigCandidate(attempt = attempt, status = "missing", fasta_path = None, gfa_path = None)

    return ContigCandidate(attempt = attempt, status = status, fasta_path = fasta_path, gfa_path = None)

def main() -> int:
    """
    Runs biological validation against every assembly attempt
    considered for one sample.

    Returns:
        int: The return code for the process (0 unless every attempt was unselectable/missing).
    """
    args = parse_args()

    logger = make_logger(
        name = "biological_validation",
        log_file_path = args.log_file,
        global_log_file_path = args.global_log_file,
    )
    context = f"{{{args.sample_id}}}"

    try:
        selection = json.loads(args.contig_selection_json.read_text(encoding = "utf-8"))
    except (OSError, json.JSONDecodeError) as e:
        logger.error(f"{context} Could not read contig selection {args.contig_selection_json}: {e}")
        return 1

    selected_attempt = selection.get("source_attempt") if selection.get("selected") else None

    blast_top_hit_identity = blast_top_hit_coverage = None
    blast_top_hits_for_conflict = None
    if args.blast_result_json is not None and args.blast_result_json.is_file():
        blast_result = json.loads(args.blast_result_json.read_text(encoding = "utf-8"))
        if blast_result.get("success") and args.blast_top_hits_tsv and args.blast_top_hits_tsv.is_file():
            from mitopipeline.utils.blast.blast_hits import read_top_blast_hits
            hits = read_top_blast_hits(args.blast_top_hits_tsv, sample_id = args.sample_id, logger = logger)
            if hits:
                blast_top_hit_identity = hits[0].pident
                blast_top_hit_coverage = hits[0].qcovs
                blast_top_hits_for_conflict = [{"sscinames": h.sscinames, "pident": h.pident} for h in hits]

    # Linear BLAST can't see across a circular genome's own origin, so
    # a correctly-assembled circular genome whose rotation differs
    # from its reference's own rotation reads as low coverage on the
    # raw top hit above. When a circular-coverage correction is
    # available, it replaces that raw value here -- never the other
    # way around, since the correction is only ever computed for a
    # circular selected attempt in the first place.
    if args.circular_coverage_json is not None and args.circular_coverage_json.is_file():
        circular_coverage = json.loads(args.circular_coverage_json.read_text(encoding = "utf-8"))
        if circular_coverage.get("success") and circular_coverage.get("corrected_query_coverage") is not None:
            logger.info(f"{context} Using circularity-corrected query coverage "
                       f"({circular_coverage['corrected_query_coverage']:.2f}%) in place of the raw top-hit value.")
            blast_top_hit_coverage = circular_coverage["corrected_query_coverage"]

    gene_names = None
    if args.annotation_result_json is not None and args.annotation_result_json.is_file():
        annotation_result = json.loads(args.annotation_result_json.read_text(encoding = "utf-8"))
        if annotation_result.get("success"):
            gene_names = annotation_result.get("gene_names") or []

    results = []
    for i, attempt_dir in enumerate(args.attempt_dir):
        candidate = _load_candidate(i, args.sample_id, attempt_dir)
        is_selected = (i == selected_attempt)
        source = "primary" if i == 0 else "rescue"

        if candidate.fasta_path is None:
            # Nothing to validate for an attempt that produced no
            # usable assembly at all -- recorded as such, not skipped
            # silently.
            logger.info(f"{context} attempt {i} ({source}): no usable assembly ('{candidate.status}'); not validated.")
            continue

        result = validate_assembly(
            sample_id = args.sample_id, attempt = i, source = source, fasta_path = candidate.fasta_path,
            circularization_status = candidate.status, min_length = args.min_length, max_length = args.max_length,
            seed_provenance = args.seed_file if is_selected else None,
            blast_top_hit_identity = blast_top_hit_identity if is_selected else None,
            blast_top_hit_coverage = blast_top_hit_coverage if is_selected else None,
            blast_top_hits_for_conflict = blast_top_hits_for_conflict if is_selected else None,
            gene_names = gene_names if is_selected else None,
            selection_rationale = selection.get("reason", "") if is_selected else "not selected",
            logger = logger,
        )
        results.append(result)

    args.output_dir.mkdir(parents = True, exist_ok = True)
    output_path = args.output_dir / f"{args.sample_id}.biological_validation.json"
    with output_path.open("w", encoding = "utf-8") as handle:
        json.dump({
            "sample_id": args.sample_id,
            "results": [r.to_dict() for r in results],
            "selected_attempt": selected_attempt,
        }, handle, indent = 4)

    if not results:
        logger.warning(f"{context} No attempt had a usable assembly to validate.")
        return 1

    return 0

if __name__ == "__main__":
    raise SystemExit(main())
