"""
models.py

Data models for biological validation of a selected mitochondrial
assembly: typed PASS/REVIEW/FAIL verdicts with reason codes,
never a bare "the tool exited 0" standing in for biological validity.
"""

# Imports
from __future__ import annotations
from dataclasses import dataclass, field

# The three verdicts a single check, or the overall assessment, can
# reach. PASS: the evidence actively supports validity. REVIEW: the
# evidence is ambiguous, missing, or borderline -- a person should
# look, this is not an automatic accept. FAIL: the evidence actively
# contradicts validity.
Verdict = str
PASS: Verdict = "PASS"
REVIEW: Verdict = "REVIEW"
FAIL: Verdict = "FAIL"

@dataclass(frozen = True, slots = True)
class CheckResult:
    """
    One individual biological check's own outcome.

    Attributes:
        name (str): A short, stable identifier for this check (e.g. "length_bounds"), used as the reason code's namespace.
        verdict (Verdict): This check's own verdict, independent of every other check.
        reason (str): A short, human-readable explanation of why this check reached that verdict.
    """
    name: str
    verdict: Verdict
    reason: str

@dataclass(frozen = True, slots = True)
class BiologicalValidationResult:
    """
    The complete biological validation outcome for one assembly
    attempt (primary or rescue) -- every individual check's result,
    the aggregated overall verdict, and the provenance a reviewer
    needs to judge a REVIEW/FAIL case without re-deriving it.

    Attributes:
        sample_id (str): The sample this result is for.
        attempt (int): Which assembly attempt this is (0 = primary, 1+ = rescue attempts).
        source (str): "primary" or "rescue" -- the same checks apply to both, this just labels which.
        overall_verdict (Verdict): The worst verdict among every check actually run (FAIL beats REVIEW beats PASS).
        checks (tuple[CheckResult, ...]): Every individual check's own result, in the order they were run.
        seed_provenance (str | None): The seed file used for this attempt, if any, or None for an unseeded run.
        tool_versions (dict[str, str]): Installed tool versions relevant to this attempt, when known.
        metrics (dict): Whatever runtime/compute metrics were available for this attempt (e.g. assembly runtime).
        selection_rationale (str): Why this candidate was selected over any others considered (contig_selection's own reason).
    """
    sample_id: str
    attempt: int
    source: str
    overall_verdict: Verdict
    checks: tuple[CheckResult, ...]
    seed_provenance: str | None
    tool_versions: dict[str, str] = field(default_factory = dict)
    metrics: dict = field(default_factory = dict)
    selection_rationale: str = ""

    def to_dict(self) -> dict:
        """
        Serializes this result to a JSON-safe dictionary.

        Args:
            None

        Returns:
            dict: The JSON-serializable representation.
        """
        return {
            "sample_id": self.sample_id,
            "attempt": self.attempt,
            "source": self.source,
            "overall_verdict": self.overall_verdict,
            "checks": [{"name": c.name, "verdict": c.verdict, "reason": c.reason} for c in self.checks],
            "seed_provenance": self.seed_provenance,
            "tool_versions": self.tool_versions,
            "metrics": self.metrics,
            "selection_rationale": self.selection_rationale,
        }
