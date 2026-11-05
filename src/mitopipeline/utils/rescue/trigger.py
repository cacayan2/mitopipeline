"""
trigger.py

Decides whether a sample's assembly should go through rescue, given its
classification from discover_assembly_artifacts() (circular/fragmented/
missing) and how many rescue attempts have already been made.
"""

# Imports
from mitopipeline.utils.rescue.models import RescueDecision

# An assembly GetOrganelle itself marked complete/circular is never
# rescued, regardless of what trigger_statuses says -- rescue exists to
# recover a WORSE outcome, never to second-guess a good one.
_NEVER_RESCUED_STATUS = "circular"

# "missing" means GetOrganelle produced no candidate assembly at all --
# so "candidate_mapped" recruitment (which maps reads against the
# existing candidate) has nothing to map against, regardless of what
# trigger_statuses configures. Same kind of hard, structural override
# as _NEVER_RESCUED_STATUS, just narrower: only this one recruitment
# mode is blocked for "missing", not the status entirely -- "full_pool"
# recruitment needs no candidate (it just hands GetOrganelle the
# complete original read pool again, often with a different seed
# locus), so it can still legitimately apply when trigger_statuses
# configures "missing" as eligible.
_MISSING_STATUS = "missing"
_CANDIDATE_REQUIRED_RECRUITMENT = "candidate_mapped"

def decide_rescue(status: str, attempt: int, config: dict,
                   recruitment: str = _CANDIDATE_REQUIRED_RECRUITMENT) -> RescueDecision:
    """
    Decides whether to attempt rescue for one sample, given its current
    assembly status and how many rescue attempts have already run.

    Retry bound: the number of configured getorganelle_profiles -- each
    rescue attempt tries the next untried profile (config.yaml's own
    documented mechanism), so once every profile has been tried once,
    rescue is exhausted. A config with no profiles configured means
    rescue is effectively off (max_attempts = 0), not an error.

    Terminal statuses: "circular" (already complete -- never rescued),
    "missing" under candidate_mapped recruitment specifically (no
    candidate to map against -- see _CANDIDATE_REQUIRED_RECRUITMENT),
    any status not in trigger_statuses (rescue was never going to
    help), and "attempts exhausted" (every configured profile already
    tried). Everything else is non-terminal: rescue should run again
    next time with the same status, until one of the above is reached.

    Args:
        status (str): The assembly's current classification, as returned by discover_assembly_artifacts() -- "circular", "fragmented", or "missing".
        attempt (int): How many rescue attempts have already been made for this sample (0-indexed; 0 means none yet).
        config (dict): This sample's assembly_rescue config section (trigger_statuses, getorganelle_profiles, etc. -- see config.yaml).
        recruitment (str): This attempt's own recruitment strategy ("candidate_mapped" or "full_pool") -- see run_rescue_mapping.py's module docstring. Defaults to "candidate_mapped" (attempt 1's own default) for callers that don't pass it.

    Returns:
        RescueDecision: Whether to rescue now, and whether the decision is terminal.
    """
    trigger_statuses = set(config.get("trigger_statuses") or [])
    max_attempts = len(config.get("getorganelle_profiles") or [])

    if status == _NEVER_RESCUED_STATUS:
        return RescueDecision(should_rescue = False, terminal = True,
                               reason = "assembly is already circular/complete; rescue not needed",
                               attempt = attempt, max_attempts = max_attempts)

    if status == _MISSING_STATUS and recruitment == _CANDIDATE_REQUIRED_RECRUITMENT:
        return RescueDecision(should_rescue = False, terminal = True,
                               reason = f"status '{status}' has no existing candidate for "
                                        f"{_CANDIDATE_REQUIRED_RECRUITMENT} recruitment to map against",
                               attempt = attempt, max_attempts = max_attempts)

    if status not in trigger_statuses:
        return RescueDecision(should_rescue = False, terminal = True,
                               reason = f"status '{status}' is not in trigger_statuses {sorted(trigger_statuses)}",
                               attempt = attempt, max_attempts = max_attempts)

    if attempt >= max_attempts:
        return RescueDecision(should_rescue = False, terminal = True,
                               reason = f"rescue attempts exhausted ({attempt} of {max_attempts} configured profiles tried)",
                               attempt = attempt, max_attempts = max_attempts)

    return RescueDecision(should_rescue = True, terminal = False,
                           reason = f"status '{status}' triggers rescue via {recruitment} recruitment; "
                                    f"attempt {attempt + 1} of {max_attempts}",
                           attempt = attempt, max_attempts = max_attempts)
