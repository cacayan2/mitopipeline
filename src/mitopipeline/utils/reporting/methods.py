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

def _assembly_sentence(report_data: SampleReportData) -> str:
    """Describes the real GetOrganelle invocation and whether rescue was needed."""
    # Deferred import avoids a circular import: markdown_report.py
    # itself imports render_methods_text from this module.
    from mitopipeline.utils.reporting.markdown_report import basename_code

    config = report_data.assembly_config
    assembly = report_data.assembly

    if config is None:
        return "Mitochondrial genome assembly was attempted with GetOrganelle, though its resolved configuration was not recorded for this run."

    seed_clause = (
        f"seeded with a reference locus ({basename_code(config.seed_file)})" if config.seeded
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
    from mitopipeline.utils.reporting.markdown_report import basename_code

    config = report_data.blast_config
    blast = report_data.blast

    # Distinguishing "never configured" (the stage itself is off)
    # from "skipped because this sample had no assembly to identify"
    # -- both read as config is None, since blast_config.smk never
    # even calls the exec script without a selected FASTA, but only
    # the first is genuinely "not configured." assembly_config being
    # present confirms GetOrganelle's own stage did run for this
    # sample (it resolves unconditionally, success or not).
    if config is None:
        if report_data.assembly_config is not None and report_data.assembly is None:
            return "Taxonomic identification via BLAST was skipped: this sample had no selected assembly to identify."
        return "Taxonomic identification via BLAST was not configured for this run."

    sentence = f"Taxonomic identification was performed with BLAST (`{config.mode}` against {basename_code(config.database)})."
    if blast is not None:
        name = blast.top_scientific_name
        if name and name.strip().upper() != "N/A":
            sentence += (
                f" The top match was `{name}` "
                f"at {blast.top_identity_percent:.2f}% identity and {blast.top_query_coverage_percent:.2f}% query coverage."
            )
        else:
            sentence += (
                f" The top match was `{blast.top_accession}` (taxon unresolved) "
                f"at {blast.top_identity_percent:.2f}% identity and {blast.top_query_coverage_percent:.2f}% query coverage."
            )
    return sentence

def _annotation_sentence(report_data: SampleReportData) -> str:
    """Describes the real MITOS2 invocation and its outcome."""
    config = report_data.annotation_config
    annotation = report_data.annotation

    if config is None:
        if report_data.assembly_config is not None and report_data.assembly is None:
            return "Gene annotation via MITOS2 was skipped: this sample had no selected assembly to annotate."
        return "Gene annotation via MITOS2 was not configured for this run."

    sentence = (
        f"Gene annotation was performed with MITOS2 "
        f"(genetic code {config.genetic_code}, reference set `{config.refseqver}`)."
    )
    if annotation is not None:
        sentence += f" {annotation.core_gene_count} of the 13 core animal mitochondrial protein-coding genes were identified"
        if annotation.trna_count is not None and annotation.rrna_count is not None:
            sentence += f", along with {annotation.trna_count} tRNA(s) and {annotation.rrna_count} rRNA(s)."
        else:
            sentence += "."
    return sentence

def _validation_sentence(report_data: SampleReportData) -> str:
    """Describes the real biological validation outcome."""
    from mitopipeline.utils.reporting.markdown_report import humanize_enum

    if report_data.overall_verdict is None:
        return ""
    return (
        f" Biological validation against assembly, annotation, and identification evidence reached a "
        f"plausibility verdict of {humanize_enum(report_data.overall_verdict)}."
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
