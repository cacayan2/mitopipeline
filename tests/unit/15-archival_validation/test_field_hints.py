"""
test_field_hints.py

Unit tests for NCBI archival field hints (Emil's own direct ask):
every real SRA/BioSample/GenBank field this pipeline can report as
missing has a concrete, actionable hint, and an uncatalogued field
name still gets a safe generic fallback rather than a crash.
"""

# Imports
from mitopipeline.utils.archival.field_hints import FIELD_HINTS, field_hint
from mitopipeline.utils.archival.models import BIOSAMPLE_REQUIRED_FIELDS, GENBANK_REQUIRED_FIELDS, SRA_REQUIRED_FIELDS

def test_every_required_field_across_all_three_record_types_has_a_real_hint() -> None:
    for field in (*SRA_REQUIRED_FIELDS, *BIOSAMPLE_REQUIRED_FIELDS, *GENBANK_REQUIRED_FIELDS):
        assert field in FIELD_HINTS, f"{field} has no catalogued hint"
        assert FIELD_HINTS[field].strip()

def test_field_hint_returns_the_real_catalogued_text() -> None:
    assert "Illumina NovaSeq 6000" in field_hint("instrument_model")
    assert "YYYY-MM-DD" in field_hint("collection_date")

def test_field_hint_falls_back_gracefully_for_an_uncatalogued_field() -> None:
    hint = field_hint("some_future_ncbi_field")
    assert hint
    assert "stage_options.archival" in hint
