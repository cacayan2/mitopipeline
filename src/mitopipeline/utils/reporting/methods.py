"""
methods.py

Generates a paper-style "Methods" paragraph (Emil's own explicit
request) describing the real tools, settings, and rescue history
actually used for one sample -- built live from the same resolved
config JSON and statistics every other report section reads, never
boilerplate text standing in for what actually ran.
"""

# Imports
from mitopipeline.utils.reporting.models import SampleReportData
from mitopipeline.utils.reporting.typst_style import render_callout_box_close, render_callout_box_open

def _assembly_sentence(report_data: SampleReportData) -> str:
    """Describes the real GetOrganelle invocation and whether rescue was needed."""
    config = report_data.assembly_config
    assembly = report_data.assembly

    if config is None:
        return "Mitochondrial genome assembly was attempted with GetOrganelle, though its resolved configuration was not recorded for this run."

    seed_clause = (
        f"seeded with a reference locus (`{config.seed_file}`)" if config.seeded
        else "run unseeded, against GetOrganelle's own default reference database"
    )
    sentence = (
        f"Mitochondrial genome assembly was performed with GetOrganelle "
        f"(organelle type `{config.organelle_type}`), {seed_clause}."
    )

    if assembly is not None and assembly.source_attempt > 0:
        sentence += (
            f" The primary assembly attempt did not yield a usable genome; a rescue reassembly "
            f"(attempt {assembly.source_attempt}, mapping recruited reads with minimap2/samtools before "
            f"re-running GetOrganelle) succeeded and was selected."
        )
    elif assembly is not None:
        sentence += " The primary assembly attempt succeeded and was selected without requiring rescue."

    return sentence

def _identification_sentence(report_data: SampleReportData) -> str:
    """Describes the real BLAST invocation and its outcome."""
    config = report_data.blast_config
    blast = report_data.blast

    if config is None:
        return "Taxonomic identification via BLAST was not configured for this run."

    sentence = f"Taxonomic identification was performed with BLAST (`{config.mode}` against `{config.database}`)."
    if blast is not None:
        sentence += (
            f" The top match was `{blast.top_scientific_name or blast.top_accession}` "
            f"at {blast.top_identity_percent:.2f}% identity and {blast.top_query_coverage_percent:.2f}% query coverage."
        )
    return sentence

def _annotation_sentence(report_data: SampleReportData) -> str:
    """Describes the real MITOS2 invocation and its outcome."""
    config = report_data.annotation_config
    annotation = report_data.annotation

    if config is None:
        return "Gene annotation via MITOS2 was not configured for this run."

    sentence = (
        f"Gene annotation was performed with MITOS2 "
        f"(genetic code {config.genetic_code}, reference set `{config.refseqver}`)."
    )
    if annotation is not None:
        sentence += (
            f" {annotation.gene_count} gene feature(s) were identified, including "
            f"{annotation.core_gene_count} of the 13 core animal mitochondrial protein-coding genes."
        )
    return sentence

def _validation_sentence(report_data: SampleReportData) -> str:
    """Describes the real biological validation outcome."""
    if report_data.overall_verdict is None:
        return ""
    return (
        f" Biological validation against assembly, annotation, and identification evidence reached an "
        f"overall verdict of {report_data.overall_verdict}."
    )

def render_methods_text(report_data: SampleReportData) -> str:
    """
    Builds the full live methods paragraph for one sample.

    Args:
        report_data (SampleReportData): The assembled report data.

    Returns:
        str: The methods paragraph, as plain prose.
    """
    sentences = [
        _assembly_sentence(report_data),
        _identification_sentence(report_data),
        _annotation_sentence(report_data),
    ]
    text = " ".join(sentences) + _validation_sentence(report_data)
    return text

def render_methods_section(report_data: SampleReportData) -> list[str]:
    """
    Renders the "Methods" section as Markdown lines, with the
    paragraph itself set inside a gray callout box matching the
    project's own real document styling.

    Args:
        report_data (SampleReportData): The assembled report data.

    Returns:
        list[str]: The Markdown lines for this section.
    """
    return [
        "## Methods", "",
        render_callout_box_open("Methods"),
        render_methods_text(report_data),
        "",
        render_callout_box_close(),
        "",
    ]
