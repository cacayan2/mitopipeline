"""
test_oom_monitor.py

Unit tests for launcher.oom_monitor -- specifically the sample-id
extraction regex, the one piece of real, novel parsing logic this
module depends on (everything else is thin /proc reads already
covered elsewhere).
"""

# Imports
from mitopipeline.launcher.oom_monitor import _SAMPLE_ID_PATTERN


def test_sample_id_pattern_matches_a_real_primary_assembly_cmdline() -> None:
    cmdline = ("/home/fishl/.../get_organelle_from_reads.py -1 r1.fq -2 r2.fq -o "
              "/home/fishl/mitogenome_project/outputs/final_output_v020/04d-assembly/"
              "culaea_inconstans_run_may2026_093/primary -F animal_mt")
    match = _SAMPLE_ID_PATTERN.search(cmdline)
    assert match is not None
    assert match.group(1) == "culaea_inconstans_run_may2026_093"


def test_sample_id_pattern_matches_a_real_rescue_reassembly_cmdline() -> None:
    cmdline = ("/home/fishl/.../get_organelle_from_reads.py -1 r1.fq -2 r2.fq -o "
              "/home/fishl/mitogenome_project/outputs/final_output_v020/07-assembly_rescue/"
              "aplodinotus_grunniens_run_may2026_066/attempt_3/reassembly -F animal_mt")
    match = _SAMPLE_ID_PATTERN.search(cmdline)
    assert match is not None
    assert match.group(1) == "aplodinotus_grunniens_run_may2026_066"


def test_sample_id_pattern_returns_none_for_an_unrelated_cmdline() -> None:
    cmdline = "/usr/bin/python3 -m some.other.tool --flag value"
    assert _SAMPLE_ID_PATTERN.search(cmdline) is None
