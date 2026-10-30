"""
gates.py

Biological validation gates for a selected mitochondrial assembly
: a set of independent, individually-reasoned checks -- FASTA
integrity, taxon-aware length bounds, circularization evidence,
ambiguous-base content, fragment status, BLAST identity/coverage,
conflicting-taxa/contamination signal, and (when annotation is
available) expected gene complement -- rolled up into one explicit
PASS/REVIEW/FAIL verdict with a reason attached to every individual
check. A clean tool exit is never treated as proof of biological
validity on its own; every claim here traces to one of these checks.

The exact same checks apply to a primary assembly attempt and every
rescue attempt -- callers distinguish them only via the "source"
label on the result, never via a different code path.
"""

# Imports
from __future__ import annotations
from logging import Logger
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from mitopipeline.utils.validation.models import FAIL, PASS, REVIEW, BiologicalValidationResult, CheckResult, Verdict

# The standard metazoan mitochondrial protein-coding gene complement
# (13 genes) -- confirmed against MITOS2's own real reference feature
# set (resources/mitos/refseq89m/featureProt/), excluding non-animal
# entries that reference set also carries (atp9, dpo, lagli, msh1,
# mttb -- fungal/plant/other, never expected in an animal mitogenome).
# Deliberately a SET, not an ordered list -- gene arrangement varies
# across taxa; only presence/absence/duplication is checked here.
CORE_ANIMAL_MT_PROTEIN_GENES = frozenset({
    "atp6", "atp8", "cob", "cox1", "cox2", "cox3",
    "nad1", "nad2", "nad3", "nad4", "nad4l", "nad5", "nad6",
})

# Length-bounds tolerance: a length just outside the configured
# complete-genome interval is flagged for review, not failed outright
# -- real mitogenomes vary in control-region length across taxa, and
# a single global/taxon default can't capture every real case exactly.
_LENGTH_REVIEW_TOLERANCE = 0.10

# Ambiguous-base (N) content thresholds.
_AMBIGUOUS_BASE_REVIEW_FRACTION = 0.01
_AMBIGUOUS_BASE_FAIL_FRACTION = 0.05

# BLAST identity/coverage thresholds for the top hit.
_BLAST_PASS_IDENTITY = 95.0
_BLAST_PASS_COVERAGE = 90.0
_BLAST_REVIEW_IDENTITY = 85.0
_BLAST_REVIEW_COVERAGE = 70.0

# Two top hits are flagged as a possible taxon conflict when they're
# from different scientific names AND close enough in identity that
# neither clearly wins.
_CONFLICTING_TAXA_IDENTITY_GAP = 2.0

# Gene-complement thresholds, as a fraction of CORE_ANIMAL_MT_PROTEIN_GENES found.
_GENE_COMPLEMENT_REVIEW_FRACTION = 0.85
_GENE_COMPLEMENT_FAIL_FRACTION = 0.70

def load_fasta_records(fasta_path: Path) -> tuple[list[SeqRecord], CheckResult]:
    """
    Parses a FASTA file and checks its basic integrity -- a real,
    nonempty file with at least one parseable record.

    Args:
        fasta_path (Path): The assembly FASTA to check.

    Returns:
        tuple[list[SeqRecord], CheckResult]: The parsed records (empty list on failure) and this check's own result.
    """
    fasta_path = Path(fasta_path)

    if not fasta_path.is_file():
        return [], CheckResult(name = "fasta_integrity", verdict = FAIL,
                               reason = f"assembly FASTA not found: {fasta_path}")

    try:
        records = list(SeqIO.parse(fasta_path, "fasta"))
    except Exception as error:
        return [], CheckResult(name = "fasta_integrity", verdict = FAIL,
                               reason = f"assembly FASTA could not be parsed: {error}")

    if not records or not any(len(record.seq) > 0 for record in records):
        return [], CheckResult(name = "fasta_integrity", verdict = FAIL,
                               reason = "assembly FASTA contains no real sequence data")

    return records, CheckResult(name = "fasta_integrity", verdict = PASS,
                                reason = f"{len(records)} record(s), {sum(len(r.seq) for r in records)}bp total")

def check_length_bounds(length_bp: int, min_length: int, max_length: int) -> CheckResult:
    """
    Checks the assembly's total length against a (taxon-aware)
    complete-genome length interval.

    Args:
        length_bp (int): The assembly's total length.
        min_length (int): The configured minimum complete-genome length for this taxon (or the default).
        max_length (int): The configured maximum complete-genome length for this taxon (or the default).

    Returns:
        CheckResult: This check's own result.
    """
    if min_length <= length_bp <= max_length:
        return CheckResult(name = "length_bounds", verdict = PASS,
                           reason = f"{length_bp}bp is within the configured interval [{min_length}, {max_length}]bp")

    tolerant_min = min_length * (1 - _LENGTH_REVIEW_TOLERANCE)
    tolerant_max = max_length * (1 + _LENGTH_REVIEW_TOLERANCE)
    if tolerant_min <= length_bp <= tolerant_max:
        return CheckResult(name = "length_bounds", verdict = REVIEW,
                           reason = f"{length_bp}bp is outside [{min_length}, {max_length}]bp but within a {_LENGTH_REVIEW_TOLERANCE:.0%} tolerance")

    return CheckResult(name = "length_bounds", verdict = FAIL,
                       reason = f"{length_bp}bp is far outside the configured interval [{min_length}, {max_length}]bp")

def check_circularization(status: str) -> CheckResult:
    """
    Checks the assembly's own circularization classification (from
    discover_assembly_artifacts()).

    Args:
        status (str): "circular", "fragmented", or "missing".

    Returns:
        CheckResult: This check's own result.
    """
    if status == "circular":
        return CheckResult(name = "circularization", verdict = PASS, reason = "assembly closed into a circular genome")
    if status == "fragmented":
        return CheckResult(name = "circularization", verdict = REVIEW, reason = "assembly did not close into a circle; fragmented/draft output")
    return CheckResult(name = "circularization", verdict = FAIL, reason = f"no usable assembly ('{status}')")

def check_ambiguous_bases(records: list[SeqRecord]) -> CheckResult:
    """
    Checks the fraction of ambiguous (N) bases across the assembly.

    Args:
        records (list[SeqRecord]): The assembly's parsed records.

    Returns:
        CheckResult: This check's own result.
    """
    total_length = sum(len(record.seq) for record in records)
    n_count = sum(str(record.seq).upper().count("N") for record in records)
    fraction = (n_count / total_length) if total_length > 0 else 1.0

    if fraction < _AMBIGUOUS_BASE_REVIEW_FRACTION:
        return CheckResult(name = "ambiguous_bases", verdict = PASS,
                           reason = f"{fraction:.2%} ambiguous (N) bases")
    if fraction < _AMBIGUOUS_BASE_FAIL_FRACTION:
        return CheckResult(name = "ambiguous_bases", verdict = REVIEW,
                           reason = f"{fraction:.2%} ambiguous (N) bases, above the {_AMBIGUOUS_BASE_REVIEW_FRACTION:.0%} review threshold")
    return CheckResult(name = "ambiguous_bases", verdict = FAIL,
                       reason = f"{fraction:.2%} ambiguous (N) bases, above the {_AMBIGUOUS_BASE_FAIL_FRACTION:.0%} fail threshold")

def check_fragment_status(records: list[SeqRecord]) -> CheckResult:
    """
    Checks how many separate sequence records the selected assembly
    actually contains -- a single record is a clean result; several
    suggests an unresolved/fragmented assembly graph, independent of
    whatever discover_assembly_artifacts() itself classified.

    Args:
        records (list[SeqRecord]): The assembly's parsed records.

    Returns:
        CheckResult: This check's own result.
    """
    count = len(records)
    if count <= 1:
        return CheckResult(name = "fragment_status", verdict = PASS, reason = "a single contiguous sequence record")
    if count <= 3:
        return CheckResult(name = "fragment_status", verdict = REVIEW, reason = f"{count} separate sequence records, not fully resolved")
    return CheckResult(name = "fragment_status", verdict = FAIL, reason = f"{count} separate sequence records, assembly graph highly unresolved")

def check_blast_identity(top_hit_identity: float | None, top_hit_coverage: float | None) -> CheckResult:
    """
    Checks the selected assembly's own top BLAST hit, when available.
    Identity alone is never accepted -- coverage must also clear the
    same tier. Missing BLAST evidence is REVIEW, not an automatic
    pass or fail -- the absence of corroborating evidence is itself
    something a reviewer should know about.

    Args:
        top_hit_identity (float | None): The top hit's percent identity, or None if no BLAST evidence is available.
        top_hit_coverage (float | None): The top hit's percent query coverage, or None if no BLAST evidence is available.

    Returns:
        CheckResult: This check's own result.
    """
    if top_hit_identity is None or top_hit_coverage is None:
        return CheckResult(name = "blast_identity", verdict = REVIEW, reason = "no BLAST evidence available to corroborate this assembly")

    if top_hit_identity >= _BLAST_PASS_IDENTITY and top_hit_coverage >= _BLAST_PASS_COVERAGE:
        return CheckResult(name = "blast_identity", verdict = PASS,
                           reason = f"top hit {top_hit_identity:.1f}% identity, {top_hit_coverage:.1f}% query coverage")
    if top_hit_identity >= _BLAST_REVIEW_IDENTITY and top_hit_coverage >= _BLAST_REVIEW_COVERAGE:
        return CheckResult(name = "blast_identity", verdict = REVIEW,
                           reason = f"top hit {top_hit_identity:.1f}% identity, {top_hit_coverage:.1f}% query coverage -- below the high-confidence tier")
    return CheckResult(name = "blast_identity", verdict = FAIL,
                       reason = f"top hit {top_hit_identity:.1f}% identity, {top_hit_coverage:.1f}% query coverage -- too weak to corroborate")

def check_conflicting_taxa(top_hits: list[dict] | None) -> CheckResult:
    """
    Checks whether the top two BLAST hits represent different
    taxa at close enough identity that neither clearly wins --
    a real contamination/misidentification signal, not something
    identity alone would surface.

    Args:
        top_hits (list[dict] | None): Ranked BLAST hits (each with "sscinames" and "pident" keys), strongest first, or None/empty if unavailable.

    Returns:
        CheckResult: This check's own result.
    """
    if not top_hits or len(top_hits) < 2:
        return CheckResult(name = "conflicting_taxa", verdict = PASS, reason = "no competing high-ranked hit to conflict with")

    first, second = top_hits[0], top_hits[1]
    same_taxon = (first.get("sscinames") or "").strip() == (second.get("sscinames") or "").strip()
    identity_gap = abs(float(first.get("pident", 0.0)) - float(second.get("pident", 0.0)))

    if not same_taxon and identity_gap < _CONFLICTING_TAXA_IDENTITY_GAP:
        return CheckResult(name = "conflicting_taxa", verdict = REVIEW,
                           reason = (f"top two hits are different taxa ({first.get('sscinames')!r} vs "
                                    f"{second.get('sscinames')!r}) within {identity_gap:.1f}% identity of each other"))
    return CheckResult(name = "conflicting_taxa", verdict = PASS,
                       reason = "no conflicting-taxa signal among the top hits")

def check_gene_complement(gene_names: list[str], expected_core_genes: frozenset[str] = CORE_ANIMAL_MT_PROTEIN_GENES) -> CheckResult:
    """
    Checks the annotated gene complement against the expected core
    set, when annotation is available. Only presence/absence/
    duplication is checked -- never arrangement or order, since real
    taxa vary there.

    Args:
        gene_names (list[str]): Every gene name MITOS2 actually annotated (duplicates included, as found).
        expected_core_genes (frozenset[str]): The expected core gene set. Defaults to the standard 13 animal mitochondrial protein-coding genes.

    Returns:
        CheckResult: This check's own result.
    """
    found = set(gene_names)
    missing = expected_core_genes - found
    duplicated = {name for name in expected_core_genes if gene_names.count(name) > 1}
    found_fraction = (len(expected_core_genes) - len(missing)) / len(expected_core_genes)

    details = f"{len(expected_core_genes) - len(missing)}/{len(expected_core_genes)} expected core genes found"
    if missing:
        details += f"; missing: {sorted(missing)}"
    if duplicated:
        details += f"; duplicated: {sorted(duplicated)}"

    if not missing and not duplicated:
        return CheckResult(name = "gene_complement", verdict = PASS, reason = details)
    if found_fraction >= _GENE_COMPLEMENT_REVIEW_FRACTION:
        return CheckResult(name = "gene_complement", verdict = REVIEW, reason = details)
    return CheckResult(name = "gene_complement", verdict = FAIL, reason = details)

def aggregate_verdict(checks: list[CheckResult]) -> Verdict:
    """
    Rolls up every individual check into one overall verdict: the
    worst verdict actually reached wins (FAIL beats REVIEW beats
    PASS). An empty check list is itself a REVIEW, never a silent
    PASS with no evidence behind it.

    Args:
        checks (list[CheckResult]): Every check actually run.

    Returns:
        Verdict: The aggregated overall verdict.
    """
    if not checks:
        return REVIEW
    verdicts = {check.verdict for check in checks}
    if FAIL in verdicts:
        return FAIL
    if REVIEW in verdicts:
        return REVIEW
    return PASS

def validate_assembly(sample_id: str, attempt: int, source: str, fasta_path: Path,
                       circularization_status: str, min_length: int, max_length: int,
                       seed_provenance: str | None = None,
                       blast_top_hit_identity: float | None = None, blast_top_hit_coverage: float | None = None,
                       blast_top_hits_for_conflict: list[dict] | None = None,
                       gene_names: list[str] | None = None,
                       expected_core_genes: frozenset[str] = CORE_ANIMAL_MT_PROTEIN_GENES,
                       tool_versions: dict[str, str] | None = None, metrics: dict | None = None,
                       selection_rationale: str = "", logger: Logger | None = None) -> BiologicalValidationResult:
    """
    Runs every applicable biological validation check against one
    selected assembly attempt and rolls them up into one overall
    verdict. The exact same checks apply whether source is "primary"
    or "rescue" -- callers never take a different code path based on
    that, only a different label on the result.

    BLAST and annotation checks are included only when their inputs
    are actually provided -- their absence is never silently treated
    as either a pass or a fail of the WHOLE assessment; it simply
    means those particular checks didn't run.

    Args:
        sample_id (str): The sample this assessment is for.
        attempt (int): Which assembly attempt this is.
        source (str): "primary" or "rescue".
        fasta_path (Path): The selected assembly's FASTA.
        circularization_status (str): "circular", "fragmented", or "missing".
        min_length (int): The configured minimum complete-genome length (taxon-aware or default).
        max_length (int): The configured maximum complete-genome length (taxon-aware or default).
        seed_provenance (str | None): The seed file used, if any.
        blast_top_hit_identity (float | None): The top BLAST hit's percent identity, if available.
        blast_top_hit_coverage (float | None): The top BLAST hit's percent query coverage, if available.
        blast_top_hits_for_conflict (list[dict] | None): The ranked top hits (for the conflicting-taxa check), if available.
        gene_names (list[str] | None): Every annotated gene name, if MITOS2 annotation is available.
        expected_core_genes (frozenset[str]): The expected core gene set for the gene-complement check.
        tool_versions (dict[str, str] | None): Installed tool versions relevant to this attempt, when known.
        metrics (dict | None): Runtime/compute metrics for this attempt, when known.
        selection_rationale (str): Why this candidate was selected (contig_selection's own reason).
        logger (Logger | None): Logger to use. Defaults to None.

    Returns:
        BiologicalValidationResult: The complete validation outcome.
    """
    records, integrity_check = load_fasta_records(fasta_path)
    checks = [integrity_check]

    if integrity_check.verdict == FAIL:
        # A FASTA that can't even be read tells us nothing else is
        # trustworthy -- every other check would just be guessing.
        overall = aggregate_verdict(checks)
        if logger is not None:
            logger.error(f"{{{sample_id}}} attempt {attempt} ({source}): {overall} -- {integrity_check.reason}")
        return BiologicalValidationResult(
            sample_id = sample_id, attempt = attempt, source = source, overall_verdict = overall,
            checks = tuple(checks), seed_provenance = seed_provenance,
            tool_versions = tool_versions or {}, metrics = metrics or {}, selection_rationale = selection_rationale,
        )

    length_bp = sum(len(record.seq) for record in records)
    checks.append(check_length_bounds(length_bp, min_length, max_length))
    checks.append(check_circularization(circularization_status))
    checks.append(check_ambiguous_bases(records))
    checks.append(check_fragment_status(records))

    if blast_top_hit_identity is not None or blast_top_hit_coverage is not None:
        checks.append(check_blast_identity(blast_top_hit_identity, blast_top_hit_coverage))
        checks.append(check_conflicting_taxa(blast_top_hits_for_conflict))

    if gene_names is not None:
        checks.append(check_gene_complement(gene_names, expected_core_genes))

    overall = aggregate_verdict(checks)
    if logger is not None:
        logger.info(f"{{{sample_id}}} attempt {attempt} ({source}): {overall} "
                    f"({sum(1 for c in checks if c.verdict == FAIL)} FAIL, "
                    f"{sum(1 for c in checks if c.verdict == REVIEW)} REVIEW, "
                    f"{sum(1 for c in checks if c.verdict == PASS)} PASS)")

    return BiologicalValidationResult(
        sample_id = sample_id, attempt = attempt, source = source, overall_verdict = overall,
        checks = tuple(checks), seed_provenance = seed_provenance,
        tool_versions = tool_versions or {}, metrics = metrics or {}, selection_rationale = selection_rationale,
    )
