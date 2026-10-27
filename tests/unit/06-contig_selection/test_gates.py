"""
test_gates.py

Unit tests for the biological validation gates (#414): every
individual check's own PASS/REVIEW/FAIL boundary, the aggregation
rule, and the end-to-end validate_assembly() orchestration.
"""

# Imports
from pathlib import Path
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from mitopipeline.utils.validation.gates import (
    CORE_ANIMAL_MT_PROTEIN_GENES,
    aggregate_verdict,
    check_ambiguous_bases,
    check_blast_identity,
    check_circularization,
    check_conflicting_taxa,
    check_fragment_status,
    check_gene_complement,
    check_length_bounds,
    load_fasta_records,
    validate_assembly,
)
from mitopipeline.utils.validation.models import FAIL, PASS, REVIEW, CheckResult

def _record(sequence: str, id: str = "r1") -> SeqRecord:
    return SeqRecord(Seq(sequence), id = id)

def _write_fasta(path: Path, sequences: list[str]) -> None:
    with path.open("w") as handle:
        for i, seq in enumerate(sequences):
            handle.write(f">seq{i}\n{seq}\n")

# ---- load_fasta_records() ----

def test_load_fasta_records_rejects_missing_file(tmp_path: Path) -> None:
    records, check = load_fasta_records(tmp_path / "does_not_exist.fasta")
    assert check.verdict == FAIL
    assert records == []

def test_load_fasta_records_rejects_empty_file(tmp_path: Path) -> None:
    path = tmp_path / "empty.fasta"
    path.write_text("")
    records, check = load_fasta_records(path)
    assert check.verdict == FAIL

def test_load_fasta_records_passes_for_real_content(tmp_path: Path) -> None:
    path = tmp_path / "real.fasta"
    _write_fasta(path, ["ACGT" * 100])
    records, check = load_fasta_records(path)
    assert check.verdict == PASS
    assert len(records) == 1

# ---- check_length_bounds() ----

def test_length_within_bounds_passes() -> None:
    assert check_length_bounds(16000, 14000, 22000).verdict == PASS

def test_length_just_outside_bounds_is_reviewed() -> None:
    # 13000 is below 14000 but within 10% tolerance (12600)
    assert check_length_bounds(13000, 14000, 22000).verdict == REVIEW

def test_length_far_outside_bounds_fails() -> None:
    assert check_length_bounds(500, 14000, 22000).verdict == FAIL

def test_length_far_above_bounds_fails() -> None:
    assert check_length_bounds(50000, 14000, 22000).verdict == FAIL

# ---- check_circularization() ----

def test_circular_status_passes() -> None:
    assert check_circularization("circular").verdict == PASS

def test_fragmented_status_is_reviewed() -> None:
    assert check_circularization("fragmented").verdict == REVIEW

def test_missing_status_fails() -> None:
    assert check_circularization("missing").verdict == FAIL

# ---- check_ambiguous_bases() ----

def test_no_ambiguous_bases_passes() -> None:
    assert check_ambiguous_bases([_record("ACGT" * 1000)]).verdict == PASS

def test_moderate_ambiguous_bases_is_reviewed() -> None:
    # 2% N content
    sequence = "N" * 20 + "ACGT" * 245
    assert check_ambiguous_bases([_record(sequence)]).verdict == REVIEW

def test_excessive_ambiguous_bases_fails() -> None:
    sequence = "N" * 100 + "ACGT" * 225
    assert check_ambiguous_bases([_record(sequence)]).verdict == FAIL

# ---- check_fragment_status() ----

def test_single_record_passes() -> None:
    assert check_fragment_status([_record("ACGT")]).verdict == PASS

def test_a_few_records_is_reviewed() -> None:
    assert check_fragment_status([_record("ACGT", "a"), _record("ACGT", "b")]).verdict == REVIEW

def test_many_records_fails() -> None:
    records = [_record("ACGT", f"r{i}") for i in range(10)]
    assert check_fragment_status(records).verdict == FAIL

# ---- check_blast_identity() ----

def test_blast_high_identity_and_coverage_passes() -> None:
    assert check_blast_identity(98.0, 95.0).verdict == PASS

def test_blast_moderate_identity_is_reviewed() -> None:
    assert check_blast_identity(90.0, 80.0).verdict == REVIEW

def test_blast_weak_identity_fails() -> None:
    assert check_blast_identity(60.0, 40.0).verdict == FAIL

def test_blast_missing_evidence_is_reviewed_not_failed() -> None:
    assert check_blast_identity(None, None).verdict == REVIEW

def test_blast_identity_alone_is_not_accepted_without_coverage() -> None:
    # High identity but very low coverage (e.g. a short conserved
    # fragment matching by chance) must not pass on identity alone.
    assert check_blast_identity(99.0, 10.0).verdict != PASS

# ---- check_conflicting_taxa() ----

def test_no_second_hit_passes() -> None:
    assert check_conflicting_taxa([{"sscinames": "Danio rerio", "pident": 99.0}]).verdict == PASS

def test_none_hits_passes() -> None:
    assert check_conflicting_taxa(None).verdict == PASS

def test_same_taxon_top_two_hits_passes() -> None:
    hits = [{"sscinames": "Danio rerio", "pident": 99.0}, {"sscinames": "Danio rerio", "pident": 98.5}]
    assert check_conflicting_taxa(hits).verdict == PASS

def test_different_taxa_close_identity_is_reviewed() -> None:
    hits = [{"sscinames": "Danio rerio", "pident": 90.0}, {"sscinames": "Danio albolineatus", "pident": 89.0}]
    assert check_conflicting_taxa(hits).verdict == REVIEW

def test_different_taxa_far_apart_identity_passes() -> None:
    hits = [{"sscinames": "Danio rerio", "pident": 99.0}, {"sscinames": "Cyprinus carpio", "pident": 70.0}]
    assert check_conflicting_taxa(hits).verdict == PASS

# ---- check_gene_complement() ----

def test_full_gene_complement_passes() -> None:
    genes = list(CORE_ANIMAL_MT_PROTEIN_GENES)
    assert check_gene_complement(genes).verdict == PASS

def test_one_missing_gene_is_reviewed() -> None:
    genes = [g for g in CORE_ANIMAL_MT_PROTEIN_GENES if g != "nad6"]
    assert check_gene_complement(genes).verdict == REVIEW

def test_many_missing_genes_fails() -> None:
    genes = list(CORE_ANIMAL_MT_PROTEIN_GENES)[:5]
    assert check_gene_complement(genes).verdict == FAIL

def test_duplicated_gene_is_flagged_even_with_full_complement() -> None:
    genes = list(CORE_ANIMAL_MT_PROTEIN_GENES) + ["cox1"]   # cox1 duplicated
    result = check_gene_complement(genes)
    assert result.verdict != PASS
    assert "duplicated" in result.reason

def test_gene_complement_never_checks_arrangement_order() -> None:
    forward = list(CORE_ANIMAL_MT_PROTEIN_GENES)
    reversed_order = list(reversed(forward))
    assert check_gene_complement(forward).verdict == check_gene_complement(reversed_order).verdict == PASS

# ---- aggregate_verdict() ----

def test_aggregate_all_pass_is_pass() -> None:
    assert aggregate_verdict([CheckResult("a", PASS, ""), CheckResult("b", PASS, "")]) == PASS

def test_aggregate_any_review_without_fail_is_review() -> None:
    assert aggregate_verdict([CheckResult("a", PASS, ""), CheckResult("b", REVIEW, "")]) == REVIEW

def test_aggregate_any_fail_wins_over_review() -> None:
    assert aggregate_verdict([CheckResult("a", REVIEW, ""), CheckResult("b", FAIL, "")]) == FAIL

def test_aggregate_empty_checks_is_review_not_pass() -> None:
    assert aggregate_verdict([]) == REVIEW

# ---- validate_assembly() -- end-to-end orchestration ----

def test_validate_assembly_bad_fasta_short_circuits_to_fail_without_other_checks(tmp_path: Path) -> None:
    result = validate_assembly(
        sample_id = "S1", attempt = 0, source = "primary", fasta_path = tmp_path / "missing.fasta",
        circularization_status = "circular", min_length = 14000, max_length = 22000,
    )
    assert result.overall_verdict == FAIL
    assert len(result.checks) == 1   # only fasta_integrity ran

def test_validate_assembly_real_pass_case(tmp_path: Path) -> None:
    fasta = tmp_path / "assembly.fasta"
    _write_fasta(fasta, ["ACGT" * 4000])   # 16000bp, within bounds
    result = validate_assembly(
        sample_id = "S1", attempt = 0, source = "primary", fasta_path = fasta,
        circularization_status = "circular", min_length = 14000, max_length = 22000,
        blast_top_hit_identity = 99.0, blast_top_hit_coverage = 98.0,
        gene_names = list(CORE_ANIMAL_MT_PROTEIN_GENES),
    )
    assert result.overall_verdict == PASS
    assert {c.name for c in result.checks} == {
        "fasta_integrity", "length_bounds", "circularization", "ambiguous_bases",
        "fragment_status", "blast_identity", "conflicting_taxa", "gene_complement",
    }

def test_validate_assembly_without_blast_or_annotation_skips_those_checks(tmp_path: Path) -> None:
    fasta = tmp_path / "assembly.fasta"
    _write_fasta(fasta, ["ACGT" * 4000])
    result = validate_assembly(
        sample_id = "S1", attempt = 0, source = "primary", fasta_path = fasta,
        circularization_status = "circular", min_length = 14000, max_length = 22000,
    )
    check_names = {c.name for c in result.checks}
    assert "blast_identity" not in check_names
    assert "gene_complement" not in check_names

def test_validate_assembly_applies_identically_to_rescue_source(tmp_path: Path) -> None:
    fasta = tmp_path / "assembly.fasta"
    _write_fasta(fasta, ["ACGT" * 4000])
    primary = validate_assembly(sample_id = "S1", attempt = 0, source = "primary", fasta_path = fasta,
                                circularization_status = "circular", min_length = 14000, max_length = 22000)
    rescue = validate_assembly(sample_id = "S1", attempt = 1, source = "rescue", fasta_path = fasta,
                               circularization_status = "circular", min_length = 14000, max_length = 22000)
    # Same checks, same inputs -- same verdicts, only the label differs.
    assert [c.verdict for c in primary.checks] == [c.verdict for c in rescue.checks]
    assert primary.source == "primary"
    assert rescue.source == "rescue"

def test_validate_assembly_records_provenance(tmp_path: Path) -> None:
    fasta = tmp_path / "assembly.fasta"
    _write_fasta(fasta, ["ACGT" * 4000])
    result = validate_assembly(
        sample_id = "S1", attempt = 0, source = "primary", fasta_path = fasta,
        circularization_status = "circular", min_length = 14000, max_length = 22000,
        seed_provenance = "seed.fasta", tool_versions = {"getorganelle": "1.7.7"},
        metrics = {"runtime_seconds": 120.5}, selection_rationale = "best by status then length",
    )
    assert result.seed_provenance == "seed.fasta"
    assert result.tool_versions == {"getorganelle": "1.7.7"}
    assert result.metrics == {"runtime_seconds": 120.5}
    assert result.selection_rationale == "best by status then length"

def test_result_to_dict_is_json_serializable(tmp_path: Path) -> None:
    import json
    fasta = tmp_path / "assembly.fasta"
    _write_fasta(fasta, ["ACGT" * 4000])
    result = validate_assembly(sample_id = "S1", attempt = 0, source = "primary", fasta_path = fasta,
                               circularization_status = "circular", min_length = 14000, max_length = 22000)
    json.dumps(result.to_dict())  # must not raise
