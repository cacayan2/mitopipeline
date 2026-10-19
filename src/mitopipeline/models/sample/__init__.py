"""
__init__.py

Re-exports Sample so "from mitopipeline.models.sample import Sample"
works -- manifest_parser.py (and anything else importing it this way)
already expected this; the package __init__ was just empty. Pre-
existing, unrelated to any of todays other changes, found only
because it blocked verifying that this session's edits import
cleanly (2026-10-19).
"""
from mitopipeline.models.sample.sample import Sample

__all__ = ["Sample"]

