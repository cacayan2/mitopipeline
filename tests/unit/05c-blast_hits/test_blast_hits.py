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

def test_rank_blast_hits_breaks_evalue_tie_by_bitscore(tmp_path: Path) -> None:
    """
    Equal e-values are broken by bitscore (higher wins).
    """
    path = tmp_path / "blast.tsv"
    path.write_text(_row(sseqid = "low_score", evalue = "1e-50", bitscore = "100") + "\n"
                    + _row(sseqid = "high_score", evalue = "1e-50", bitscore = "200") + "\n")
    ranked = rank_blast_hits(parse_blast_tsv(path))
    assert ranked[0].sseqid == "high_score"

def test_rank_blast_hits_breaks_bitscore_tie_by_qcovs(tmp_path: Path) -> None:
    """
    Equal e-value and bitscore are broken by query coverage (higher wins).
    """
    path = tmp_path / "blast.tsv"
    path.write_text(_row(sseqid = "low_cov", evalue = "1e-50", bitscore = "200", qcovs = "50") + "\n"
                    + _row(sseqid = "high_cov", evalue = "1e-50", bitscore = "200", qcovs = "95") + "\n")
    ranked = rank_blast_hits(parse_blast_tsv(path))
    assert ranked[0].sseqid == "high_cov"

def test_rank_blast_hits_breaks_qcovs_tie_by_pident(tmp_path: Path) -> None:
    """
    Equal e-value, bitscore, and coverage are broken by percent
    identity (higher wins).
    """
    path = tmp_path / "blast.tsv"
    path.write_text(_row(sseqid = "low_ident", evalue = "1e-50", bitscore = "200", qcovs = "95", pident = "90.0") + "\n"
                    + _row(sseqid = "high_ident", evalue = "1e-50", bitscore = "200", qcovs = "95", pident = "99.5") + "\n")
    ranked = rank_blast_hits(parse_blast_tsv(path))
    assert ranked[0].sseqid == "high_ident"

def test_rank_blast_hits_breaks_pident_tie_by_length(tmp_path: Path) -> None:
    """
    Equal e-value, bitscore, coverage, and identity are broken by
    alignment length (longer wins).
    """
    path = tmp_path / "blast.tsv"
    path.write_text(_row(sseqid = "short", evalue = "1e-50", bitscore = "200", qcovs = "95", pident = "99.0", length = "400") + "\n"
                    + _row(sseqid = "long", evalue = "1e-50", bitscore = "200", qcovs = "95", pident = "99.0", length = "800") + "\n")
    ranked = rank_blast_hits(parse_blast_tsv(path))
    assert ranked[0].sseqid == "long"

def test_rank_blast_hits_breaks_full_tie_alphabetically_by_accession(tmp_path: Path) -> None:
    """
    A complete tie across every other criterion falls back to the
    subject accession itself, alphabetically -- fully deterministic,
    never arbitrary dict/set ordering.
    """
    path = tmp_path / "blast.tsv"
    path.write_text(_row(sseqid = "zzz_subject") + "\n" + _row(sseqid = "aaa_subject") + "\n")
    ranked = rank_blast_hits(parse_blast_tsv(path))
    assert [hit.sseqid for hit in ranked] == ["aaa_subject", "zzz_subject"]

def test_rank_blast_hits_deduplicates_repeated_accession_keeping_the_strongest(tmp_path: Path) -> None:
    """
    The same accession appearing in 3 separate rows (e.g. multiple
    local alignments against the same reference) is deduplicated down
    to exactly one hit -- the strongest of the three, not the first
    or last seen.
    """
    path = tmp_path / "blast.tsv"
    path.write_text(
        _row(sseqid = "dup", evalue = "1e-10", bitscore = "80") + "\n"
        + _row(sseqid = "dup", evalue = "1e-90", bitscore = "300") + "\n"
        + _row(sseqid = "dup", evalue = "1e-30", bitscore = "150") + "\n"
    )
    ranked = rank_blast_hits(parse_blast_tsv(path))
    assert len(ranked) == 1
    assert ranked[0].bitscore == 300

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
