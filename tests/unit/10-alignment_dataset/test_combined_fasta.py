"""
test_combined_fasta.py

Unit tests for building the combined, batch-level phylogeny FASTA
(#344/#345): duplicate IDs, rejected/missing assemblies, and
deterministic ordering.
"""

# Imports
from pathlib import Path
import pytest
from mitopipeline.utils.phylogeny.combined_fasta import CombinedFastaEntry, build_combined_fasta

def _write_fasta(path: Path, record_id: str, sequence: str = "ACGT" * 10) -> Path:
    path.write_text(f">{record_id}\n{sequence}\n")
    return path

def test_included_samples_are_written_in_deterministic_sorted_order(tmp_path: Path) -> None:
    entries = [
        CombinedFastaEntry("zebra", _write_fasta(tmp_path / "zebra.fasta", "zebra_native"), True),
        CombinedFastaEntry("alpha", _write_fasta(tmp_path / "alpha.fasta", "alpha_native"), True),
        CombinedFastaEntry("mid", _write_fasta(tmp_path / "mid.fasta", "mid_native"), True),
    ]
    result = build_combined_fasta(entries, tmp_path / "combined.fasta")
    assert result.included_sample_ids == ("alpha", "mid", "zebra")

def test_order_is_deterministic_regardless_of_input_order(tmp_path: Path) -> None:
    a = [
        CombinedFastaEntry("s1", _write_fasta(tmp_path / "s1.fasta", "s1_native"), True),
        CombinedFastaEntry("s2", _write_fasta(tmp_path / "s2.fasta", "s2_native"), True),
    ]
    b = list(reversed(a))
    result_a = build_combined_fasta(a, tmp_path / "combined_a.fasta")
    result_b = build_combined_fasta(b, tmp_path / "combined_b.fasta")
    assert result_a.included_sample_ids == result_b.included_sample_ids

def test_rejected_assembly_is_excluded_and_recorded(tmp_path: Path) -> None:
    entries = [
        CombinedFastaEntry("good", _write_fasta(tmp_path / "good.fasta", "good_native"), True),
        CombinedFastaEntry("bad", _write_fasta(tmp_path / "bad.fasta", "bad_native"), False),
    ]
    result = build_combined_fasta(entries, tmp_path / "combined.fasta")
    assert result.included_sample_ids == ("good",)
    assert result.excluded_sample_ids == ("bad",)
    assert "not biologically accepted" in result.exclusion_reasons["bad"]

def test_missing_fasta_path_is_excluded_and_recorded(tmp_path: Path) -> None:
    entries = [
        CombinedFastaEntry("good", _write_fasta(tmp_path / "good.fasta", "good_native"), True),
        CombinedFastaEntry("none_selected", None, True),
    ]
    result = build_combined_fasta(entries, tmp_path / "combined.fasta")
    assert "none_selected" in result.excluded_sample_ids
    assert "no selected assembly FASTA" in result.exclusion_reasons["none_selected"]

def test_nonexistent_fasta_file_is_excluded_and_recorded(tmp_path: Path) -> None:
    entries = [
        CombinedFastaEntry("good", _write_fasta(tmp_path / "good.fasta", "good_native"), True),
        CombinedFastaEntry("missing_file", tmp_path / "does_not_exist.fasta", True),
    ]
    result = build_combined_fasta(entries, tmp_path / "combined.fasta")
    assert "missing_file" in result.excluded_sample_ids

def test_duplicate_sample_ids_raise(tmp_path: Path) -> None:
    entries = [
        CombinedFastaEntry("dup", _write_fasta(tmp_path / "dup1.fasta", "dup1_native"), True),
        CombinedFastaEntry("dup", _write_fasta(tmp_path / "dup2.fasta", "dup2_native"), True),
    ]
    with pytest.raises(ValueError, match = "Duplicate sample_id"):
        build_combined_fasta(entries, tmp_path / "combined.fasta")

def test_zero_included_samples_raises(tmp_path: Path) -> None:
    entries = [CombinedFastaEntry("bad", _write_fasta(tmp_path / "bad.fasta", "bad_native"), False)]
    with pytest.raises(ValueError, match = "No samples were included"):
        build_combined_fasta(entries, tmp_path / "combined.fasta")

def test_headers_are_the_sample_id_not_the_original_record_id(tmp_path: Path) -> None:
    entries = [CombinedFastaEntry("sample_001", _write_fasta(tmp_path / "f.fasta", "some_native_header"), True)]
    output = tmp_path / "combined.fasta"
    build_combined_fasta(entries, output)
    content = output.read_text()
    assert content.startswith(">sample_001\n")
    assert "some_native_header" not in content

def test_multi_record_fasta_contributes_only_the_first_record(tmp_path: Path) -> None:
    multi = tmp_path / "multi.fasta"
    multi.write_text(">frag1\nACGTACGT\n>frag2\nTTTTTTTT\n")
    entries = [CombinedFastaEntry("fragmented_sample", multi, True)]
    output = tmp_path / "combined.fasta"
    build_combined_fasta(entries, output)
    records_written = output.read_text().count(">")
    assert records_written == 1

def test_empty_fasta_file_is_excluded_and_recorded(tmp_path: Path) -> None:
    empty = tmp_path / "empty.fasta"
    empty.write_text("")
    entries = [
        CombinedFastaEntry("good", _write_fasta(tmp_path / "good.fasta", "good_native"), True),
        CombinedFastaEntry("empty_one", empty, True),
    ]
    result = build_combined_fasta(entries, tmp_path / "combined.fasta")
    assert "empty_one" in result.excluded_sample_ids
