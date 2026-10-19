"""
test_blast_hits.py

Unit tests for parsing, ranking, and selecting BLAST hits.
"""

# Imports
from pathlib import Path
import pytest
from mitopipeline.utils.blast.blast_hits import (
    parse_blast_tsv,
    rank_blast_hits,
    select_top_blast_hits,
    write_top_blast_hits,
)
from mitopipeline.utils.blast.models import BlastHit

def _row(sseqid: str = "subj1", pident: str = "98.5", evalue: str = "1e-50", bitscore: str = "200",
         qcovs: str = "95", length: str = "500") -> str:
    """
    Builds one well-formed tabular BLAST row as a tab-separated string,
    overridable per test.
    """
    return "\t".join([
        "query1", sseqid, pident, length, "2", "0", "1", str(int(length) + 1), "1", str(int(length) + 1),
        evalue, bitscore, qcovs, "9606", "Homo sapiens", "Homo sapiens mitochondrion, complete genome",
    ])

def test_parse_blast_tsv_well_formed_row(tmp_path: Path) -> None:
    """
    A well-formed row parses into a BlastHit with the right types.
    """
    path = tmp_path / "blast.tsv"
    path.write_text(_row() + "\n")
    hits = parse_blast_tsv(path)
    assert len(hits) == 1
    hit = hits[0]
    assert hit.sseqid == "subj1"
    assert hit.pident == 98.5
    assert hit.length == 500
    assert hit.evalue == 1e-50
    assert hit.rank is None

def test_parse_blast_tsv_skips_blank_lines(tmp_path: Path) -> None:
    """
    Blank lines between real rows are skipped, not treated as
    malformed rows.
    """
    path = tmp_path / "blast.tsv"
    path.write_text(_row() + "\n\n" + _row(sseqid = "subj2") + "\n")
    hits = parse_blast_tsv(path)
    assert len(hits) == 2

def test_parse_blast_tsv_empty_file(tmp_path: Path) -> None:
    """
    An empty file (zero BLAST hits) is valid, not an error -- it
    parses to an empty list.
    """
    path = tmp_path / "blast.tsv"
    path.write_text("")
    assert parse_blast_tsv(path) == []

def test_parse_blast_tsv_missing_file(tmp_path: Path) -> None:
    """
    A missing file raises FileNotFoundError.
    """
    with pytest.raises(FileNotFoundError):
        parse_blast_tsv(tmp_path / "does_not_exist.tsv")

def test_parse_blast_tsv_rejects_directory(tmp_path: Path) -> None:
    """
    A path that is a directory, not a file, raises ValueError.
    """
    with pytest.raises(ValueError):
        parse_blast_tsv(tmp_path)

def test_parse_blast_tsv_rejects_wrong_column_count(tmp_path: Path) -> None:
    """
    A row with too few or too many columns raises ValueError naming
    the offending line.
    """
    path = tmp_path / "blast.tsv"
    path.write_text("query1\tsubj1\t98.5\n")
    with pytest.raises(ValueError, match = "row 1"):
        parse_blast_tsv(path)

def test_parse_blast_tsv_rejects_malformed_numeric_value(tmp_path: Path) -> None:
    """
    A row whose numeric column can't be parsed (e.g. "N/A" where a
    float is expected) raises ValueError naming the offending line.
    """
    path = tmp_path / "blast.tsv"
    path.write_text(_row(pident = "N/A") + "\n")
    with pytest.raises(ValueError, match = "row 1"):
        parse_blast_tsv(path)

def test_rank_blast_hits_keeps_best_per_subject(tmp_path: Path) -> None:
    """
    Only the strongest hit per subject is kept; a weaker duplicate
    hit against the same subject is dropped.
    """
    path = tmp_path / "blast.tsv"
    path.write_text(_row(sseqid = "subj1", bitscore = "100") + "\n" + _row(sseqid = "subj1", bitscore = "200") + "\n")
    hits = parse_blast_tsv(path)
    ranked = rank_blast_hits(hits)
    assert len(ranked) == 1
    assert ranked[0].bitscore == 200

def test_rank_blast_hits_orders_strongest_first(tmp_path: Path) -> None:
    """
    Hits against different subjects are sorted strongest match first
    (lower e-value wins).
    """
    path = tmp_path / "blast.tsv"
    path.write_text(_row(sseqid = "weak", evalue = "1e-5") + "\n" + _row(sseqid = "strong", evalue = "1e-80") + "\n")
    hits = parse_blast_tsv(path)
    ranked = rank_blast_hits(hits)
    assert [hit.sseqid for hit in ranked] == ["strong", "weak"]

def test_select_top_blast_hits_assigns_rank(tmp_path: Path) -> None:
    """
    Selected hits are numbered starting at 1, strongest first.
    """
    path = tmp_path / "blast.tsv"
    path.write_text(_row(sseqid = "a", evalue = "1e-80") + "\n" + _row(sseqid = "b", evalue = "1e-5") + "\n")
    hits = parse_blast_tsv(path)
    selected = select_top_blast_hits(hits, maximum_hits = 6)
    assert [(hit.sseqid, hit.rank) for hit in selected] == [("a", 1), ("b", 2)]

def test_select_top_blast_hits_caps_at_maximum(tmp_path: Path) -> None:
    """
    No more than maximum_hits are returned, even with more candidates
    available.
    """
    rows = "\n".join(_row(sseqid = f"subj{i}") for i in range(10))
    path = tmp_path / "blast.tsv"
    path.write_text(rows + "\n")
    hits = parse_blast_tsv(path)
    selected = select_top_blast_hits(hits, maximum_hits = 3)
    assert len(selected) == 3

def test_select_top_blast_hits_rejects_non_positive_maximum() -> None:
    """
    maximum_hits must be greater than 0.
    """
    with pytest.raises(ValueError):
        select_top_blast_hits([], maximum_hits = 0)

def test_write_top_blast_hits_round_trips(tmp_path: Path) -> None:
    """
    Writing selected hits produces a TSV with the sample ID and rank
    columns present and correct.
    """
    hit = BlastHit(qseqid = "q1", sseqid = "s1", pident = 99.0, length = 500, mismatch = 1, gapopen = 0,
                    qstart = 1, qend = 500, sstart = 1, send = 500, evalue = 1e-50, bitscore = 200.0,
                    qcovs = 95.0, staxids = "9606", sscinames = "Homo sapiens", stitle = "title", rank = 1)
    output_path = tmp_path / "top_hits.tsv"
    write_top_blast_hits([hit], output_path, sample_id = "SAMPLE_01")
    lines = output_path.read_text().strip().splitlines()
    assert lines[0].split("\t")[:2] == ["sample_id", "rank"]
    assert lines[1].split("\t")[:2] == ["SAMPLE_01", "1"]
