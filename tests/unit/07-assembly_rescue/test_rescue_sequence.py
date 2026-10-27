"""
test_rescue_sequence.py

Integration-level tests for rescue termination and retry limits (#427):
decide_rescue() is already unit-tested per call, but Snakemake drives
the real multi-attempt sequence by calling it again each time with the
PREVIOUS attempt's own outcome -- these tests simulate that sequence
directly to confirm the three real termination paths actually chain
together correctly: no improvement across every configured profile,
successful circularization partway through, and a tool failure during
an attempt.
"""

# Imports
from mitopipeline.utils.rescue.trigger import decide_rescue

CONFIG = {"trigger_statuses": ["fragmented", "partial"], "getorganelle_profiles": [None, None, None]}

def _simulate_sequence(statuses: list[str]) -> list:
    """
    Drives decide_rescue() across a simulated sequence of reassembly
    outcomes, one call per attempt -- statuses[i] is what the PREVIOUS
    attempt (or the primary assembly, for i=0) actually produced.
    Stops calling decide_rescue() the moment a decision is terminal,
    exactly as Snakemake would stop requesting another attempt.

    Args:
        statuses (list[str]): The assembly status feeding each successive decide_rescue() call.

    Returns:
        list: The decisions actually made, in order, up to and including the first terminal one.
    """
    decisions = []
    for attempt, status in enumerate(statuses):
        decision = decide_rescue(status, attempt, CONFIG)
        decisions.append(decision)
        if decision.terminal:
            break
    return decisions

def test_no_improvement_exhausts_every_configured_profile() -> None:
    # Every attempt stays fragmented -- rescue should keep trying
    # through all 3 configured profiles, then stop.
    decisions = _simulate_sequence(["fragmented", "fragmented", "fragmented", "fragmented"])
    assert [d.should_rescue for d in decisions] == [True, True, True, False]
    assert decisions[-1].terminal is True
    assert "exhausted" in decisions[-1].reason

def test_successful_circularization_stops_further_attempts() -> None:
    # Attempt 0: fragmented, rescue triggers. Attempt 1's reassembly
    # actually succeeds and produces a circular genome -- the next
    # decide_rescue() call must see that and stop immediately, never
    # consuming a 3rd profile it no longer needs.
    decisions = _simulate_sequence(["fragmented", "circular"])
    assert len(decisions) == 2
    assert decisions[0].should_rescue is True
    assert decisions[1].should_rescue is False
    assert decisions[1].terminal is True
    assert "already circular" in decisions[1].reason

def test_tool_failure_mid_sequence_does_not_fabricate_improvement() -> None:
    # A rescue attempt's own tool chain fails outright (e.g. GetOrganelle
    # itself crashes) -- run_rescue_reassembly.py classifies that as
    # "missing" (see _write_status's success=False path), never as an
    # improvement. The next decide_rescue() call must still see
    # "missing" (a trigger status is NOT required here -- "missing" is
    # its own real classification, distinct from "fragmented"/"partial")
    # and, since "missing" isn't configured as a trigger status in this
    # sequence, rescue correctly stops rather than retrying forever on a
    # tool that keeps crashing.
    decisions = _simulate_sequence(["fragmented", "missing"])
    assert decisions[0].should_rescue is True
    assert decisions[1].should_rescue is False
    assert decisions[1].terminal is True
    assert "not in trigger_statuses" in decisions[1].reason

def test_missing_configured_as_a_trigger_keeps_retrying_through_a_failure() -> None:
    # The opposite configuration: if "missing" IS a configured trigger
    # status, a tool failure classified as "missing" is treated the
    # same as "fragmented" -- rescue keeps retrying rather than giving
    # up on the first tool failure, up to the same profile-count bound.
    config = {"trigger_statuses": ["fragmented", "partial", "missing"], "getorganelle_profiles": [None, None]}
    decisions = [decide_rescue(s, i, config) for i, s in enumerate(["fragmented", "missing", "fragmented"])]
    assert [d.should_rescue for d in decisions] == [True, True, False]
    assert decisions[-1].terminal is True
    assert "exhausted" in decisions[-1].reason
