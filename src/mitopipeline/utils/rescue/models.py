"""
models.py

Data models for assembly-rescue decision-making.
"""

# Imports
from dataclasses import dataclass

@dataclass(frozen = True, slots = True)
class RescueDecision:
    """
    The outcome of deciding whether a sample's assembly should go
    through rescue, given its current classification and how many
    rescue attempts have already been made.

    Attributes:
        should_rescue (bool): Whether a rescue attempt should be made now.
        terminal (bool): Whether this decision is final -- no further rescue attempts should ever be made for this sample, regardless of future calls (a successful circular assembly, an untriggered status, or exhausted retry bounds all terminate).
        reason (str): A short, human-readable explanation of the decision, for logging.
        attempt (int): The rescue attempt number this decision corresponds to (0-indexed, same value passed in).
        max_attempts (int): The configured retry bound this decision was evaluated against.
    """
    should_rescue: bool
    terminal: bool
    reason: str
    attempt: int
    max_attempts: int
