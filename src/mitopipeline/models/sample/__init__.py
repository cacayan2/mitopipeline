"""
__init__.py

Re-exports Sample so "from mitopipeline.models.sample import Sample"
works, matching how manifest_parser.py (and anything else importing
it this way) expects to import it.
"""
from mitopipeline.models.sample.sample import Sample

__all__ = ["Sample"]

