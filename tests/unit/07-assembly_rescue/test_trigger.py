"""
test_trigger.py

Unit tests for the assembly-rescue trigger decision function (#405).
"""

# Imports
import pytest
from mitopipeline.utils.rescue.trigger import decide_rescue

CONFIG = {
    "trigger_statuses": ["fragmented", "missing"],
    "getorganelle_profiles": [
        {"name": "reduced_coverage"},
        {"name": "alternate_kmers"},
    ],
}

def test_circular_never_rescued_even_if_configured_as_a_trigger() -> None:
    config = dict(CONFIG, trigger_statuses = ["circular", "fragmented"])
    decision = decide_rescue("circular", attempt = 0, config = config)
    assert decision.should_rescue is False
    assert decision.terminal is True
    assert "complete" in decision.reason

def test_fragmented_triggers_rescue_on_first_attempt() -> None:
    decision = decide_rescue("fragmented", attempt = 0, config = CONFIG)
    assert decision.should_rescue is True
    assert decision.terminal is False
    assert decision.max_attempts == 2

def test_missing_triggers_rescue_same_as_fragmented() -> None:
    decision = decide_rescue("missing", attempt = 0, config = CONFIG)
    assert decision.should_rescue is True
    assert decision.terminal is False

def test_status_not_in_trigger_statuses_is_terminal_and_not_rescued() -> None:
    config = dict(CONFIG, trigger_statuses = ["fragmented"])
    decision = decide_rescue("missing", attempt = 0, config = config)
    assert decision.should_rescue is False
    assert decision.terminal is True
    assert "not in trigger_statuses" in decision.reason

def test_second_attempt_still_within_bounds_rescues_again() -> None:
    decision = decide_rescue("fragmented", attempt = 1, config = CONFIG)
    assert decision.should_rescue is True
    assert decision.terminal is False
    assert "2 of 2" in decision.reason

def test_attempts_exhausted_at_profile_count_is_terminal() -> None:
    decision = decide_rescue("fragmented", attempt = 2, config = CONFIG)
    assert decision.should_rescue is False
    assert decision.terminal is True
    assert "exhausted" in decision.reason

def test_no_profiles_configured_means_rescue_is_effectively_off() -> None:
    config = dict(CONFIG, getorganelle_profiles = [])
    decision = decide_rescue("fragmented", attempt = 0, config = config)
    assert decision.should_rescue is False
    assert decision.terminal is True
    assert decision.max_attempts == 0

def test_empty_trigger_statuses_rescues_nothing() -> None:
    config = dict(CONFIG, trigger_statuses = [])
    for status in ("fragmented", "missing"):
        decision = decide_rescue(status, attempt = 0, config = config)
        assert decision.should_rescue is False
        assert decision.terminal is True

@pytest.mark.parametrize("status", ["fragmented", "missing", "circular"])
def test_decision_always_echoes_back_the_attempt_number(status: str) -> None:
    decision = decide_rescue(status, attempt = 1, config = CONFIG)
    assert decision.attempt == 1
