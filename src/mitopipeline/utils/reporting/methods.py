"""
methods.py

Generates a paper-style "Methods" paragraph (Emil's own explicit
request) describing the real tools, settings, databases, and rescue
history actually used for one sample -- built live from the same
resolved config JSON and statistics every other report section
reads, never boilerplate text standing in for what actually ran.
Concise by design (a handful of sentences), but covers every real
stage: trimming/QC, assembly (including every real rescue attempt),
identification, annotation, visualization, phylogeny, and validation.
"""

# Imports
import re
from mitopipeline.utils.reporting.compute_metrics import tool_version_for_stage
from mitopipeline.utils.reporting.models import SampleReportData

def _short_version(raw: str) -> str:
    """
    Extracts a clean version-number-like substring (e.g. "1.4.0",
    "2.31-r1302") from a tool's own real raw --version output, for
    compact inline citation next to a tool's own name. Falls back to
    the full raw string when no such pattern is found -- the
    extraction only ever narrows what's shown, it never invents a
    version that wasn't really there.
    """
    match = re.search(r"\d+(?:\.\d+)+(?:-\S+)?\+?", raw)
    return match.group(0) if match else raw

def _version_clause(report_data: SampleReportData, stage: str, tool_name: str) -> str:
    """Returns " vX.Y.Z" (leading space, ready to splice right after a tool name) for a real detected tool version, or "" when none was ever recorded for that real stage/tool."""
    raw = tool_version_for_stage(report_data.compute_summary, stage, tool_name)
    return f" {_short_version(raw)}" if raw else ""

def _trimming_sentence(report_data: SampleReportData) -> str:
    """Describes the real fastp settings used and the real FastQC verdicts recorded at each stage."""
    from mitopipeline.utils.reporting.markdown_report import format_percent, humanize_enum

    stats = report_data.fastp_stats
    if stats is None:
        return "Read trimming and quality filtering via fastp was not configured for this run."

    settings = []
    if "qualified_quality_phred" in stats:
        settings.append(f"minimum base quality Q{stats['qualified_quality_phred']}")
    if "length_required" in stats:
        settings.append(f"minimum length {stats['length_required']}bp")
    if stats.get("adapter_detection"):
        settings.append("automatic paired-end adapter detection")
    settings_clause = f" ({'; '.join(settings)})" if settings else ""
    fastp_version = _version_clause(report_data, "trimming", "fastp")

    sentence = (
        f"Raw reads were quality-trimmed with fastp{fastp_version}{settings_clause}, retaining "
        f"{format_percent(stats['read_retention_rate'], fraction = True)} of input reads."
    )

    qc_clauses = [
        f"{label} {humanize_enum(qc.overall_verdict).lower()}"
        for label, qc in (
            ("raw", report_data.raw_read_qc),
            ("trimmed", report_data.trimmed_read_qc),
            ("recruited", report_data.recruited_read_qc),
        )
        if qc is not None
    ]
    if qc_clauses:
        fastqc_version = (
            _version_clause(report_data, "qc/raw", "fastqc")
            or _version_clause(report_data, "qc/trimmed", "fastqc")
        )
        sentence += f" FastQC{fastqc_version} quality control was performed at each stage ({', '.join(qc_clauses)})."

    return sentence

def _assembly_sentence(report_data: SampleReportData) -> str:
    """Describes the real GetOrganelle invocation, the real seed used (or why not), and every real rescue attempt made."""
    # Deferred import avoids a circular import: markdown_report.py
    # itself imports render_methods_text from this module.
    from mitopipeline.utils.reporting.markdown_report import basename_code, humanize_enum

    config = report_data.assembly_config
    assembly = report_data.assembly

    if config is None:
        return "Mitochondrial genome assembly was attempted with GetOrganelle, though its resolved configuration was not recorded for this run."

    if config.seeded and config.seed_gene and config.seed_accession and config.seed_source == "mitofish":
        seed_clause = f"seeded with a real MitoFish reference locus ({config.seed_gene}, accession `{config.seed_accession}`)"
    elif config.seeded and config.seed_file:
        seed_clause = f"seeded with a reference locus ({basename_code(config.seed_file)})"
    elif config.seed_fallback_reason:
        seed_clause = (
            f"run unseeded, against GetOrganelle's own default reference database "
            f"(no usable MitoFish seed: {humanize_enum(config.seed_fallback_reason).lower()})"
        )
    else:
        seed_clause = "run unseeded, against GetOrganelle's own default reference database"

    primary_version = _version_clause(report_data, "04d-assembly", "getorganelle")
    sentence = (
        f"Mitochondrial genome assembly was performed with GetOrganelle{primary_version} "
        f"(organelle type `{config.organelle_type}`), {seed_clause}."
    )

    # Rescue mapping/reassembly genuinely runs as its own invocation
    # under 07-assembly_rescue -- cited separately from the primary
    # attempt's own 04d-assembly invocation above, since the two can
    # in principle resolve to different installed versions even
    # though they're almost always identical in practice.
    minimap2_version = _version_clause(report_data, "07-assembly_rescue", "minimap2")
    samtools_version = _version_clause(report_data, "07-assembly_rescue", "samtools")
    rescue_getorganelle_version = _version_clause(report_data, "07-assembly_rescue", "getorganelle")

    attempts = report_data.rescue_attempts
    if attempts and any(a.attempted for a in attempts):
        descriptions = []
        for a in attempts:
            if a.attempted:
                outcome = f"produced a {humanize_enum(a.status).lower()} result" if a.status else "did not complete"
                profile_clause = f" (profile `{a.profile_name}`)" if a.profile_name else ""
                descriptions.append(f"attempt {a.attempt}{profile_clause} {outcome}")
            else:
                descriptions.append(f"attempt {a.attempt} was not run ({a.reason})")
        sentence += (
            " The primary attempt did not yield a usable circular or fragmented genome, triggering rescue "
            f"reassembly (mapping recruited reads with minimap2{minimap2_version}/samtools{samtools_version} "
            f"before re-running GetOrganelle{rescue_getorganelle_version}): "
            f"{'; '.join(descriptions)}."
        )
        if assembly is not None and assembly.source_attempt > 0:
            sentence += f" Rescue attempt {assembly.source_attempt} was selected."
        elif assembly is not None:
            sentence += " None improved on the primary result, so the primary assembly was retained."
    elif assembly is not None and assembly.source_attempt > 0:
        sentence += (
            f" The primary assembly attempt did not yield a usable genome; a rescue reassembly "
            f"(attempt {assembly.source_attempt}, mapping recruited reads with minimap2{minimap2_version}/"
            f"samtools{samtools_version} before re-running GetOrganelle{rescue_getorganelle_version}) "
            f"succeeded and was selected."
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

    blast_version = _version_clause(report_data, "05b-blast", "blast")
    sentence = f"Taxonomic identification was performed with BLAST{blast_version} (`{config.mode}` against {basename_code(config.database)})."
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

    assembly_config = report_data.assembly_config
    biomarker_clause = ""
    if assembly_config is not None and assembly_config.seed_gene and assembly_config.seed_source == "mitofish":
        biomarker_clause = f", MitoFish biomarker {assembly_config.seed_gene}"

    mitos2_version = _version_clause(report_data, "08b-annotation", "mitos2")
    sentence = (
        f"Gene annotation was performed with MITOS2{mitos2_version} "
        f"(genetic code {config.genetic_code}{biomarker_clause}, reference set `{config.refseqver}`)."
    )
    if annotation is not None:
        sentence += f" {annotation.core_gene_count} of the 13 core animal mitochondrial protein-coding genes were identified"
        if annotation.trna_count is not None and annotation.rrna_count is not None:
            sentence += f", along with {annotation.trna_count} tRNA(s) and {annotation.rrna_count} rRNA(s)."
        else:
            sentence += "."
    return sentence

def _visualization_sentence(report_data: SampleReportData) -> str:
    """Notes whether an annotated circular genome map was rendered."""
    if report_data.visualization_png is None:
        return ""
    return " An annotated circular genome map was rendered from the real MITOS2 gene coordinates."

def _phylogeny_sentence(report_data: SampleReportData) -> str:
    """Describes the real phylogenetic placement, or why one wasn't built."""
    phylogeny = report_data.phylogeny
    if phylogeny is None:
        return ""
    if not phylogeny.success:
        return f" Phylogenetic placement was not completed: {phylogeny.reason}"
    model_clause = f" under the `{phylogeny.model}` substitution model" if phylogeny.model else ""
    tip_clause = f" ({phylogeny.tip_count} total tips)" if phylogeny.tip_count is not None else ""
    mafft_version = _version_clause(report_data, "11-phylogeny", "mafft")
    iqtree_version = _version_clause(report_data, "11-phylogeny", "iqtree")
    return (
        f" Sequences were aligned with MAFFT{mafft_version}, and a maximum-likelihood phylogeny was inferred "
        f"with IQ-TREE{iqtree_version}{model_clause}, placing the assembled genome among its own top BLAST "
        f"reference matches{tip_clause}."
    )

def _validation_sentence(report_data: SampleReportData) -> str:
    """Describes the real biological validation outcome."""
    from mitopipeline.utils.reporting.markdown_report import humanize_enum

    if report_data.overall_verdict is None:
        return ""
    return (
        f" Biological validation against assembly, annotation, and identification evidence reached a "
        f"plausibility verdict of {humanize_enum(report_data.overall_verdict)}."
    )

def _environment_sentence(report_data: SampleReportData) -> str:
    """Describes the real compute environment this sample was processed on (non-identifying specs only)."""
    info = report_data.machine_info
    if info is None:
        return ""

    cpu_clause = f"{info.cpu_model}" if info.cpu_model else f"a {info.architecture} CPU"
    if info.cpu_physical_cores is not None and info.cpu_logical_threads is not None:
        cpu_clause += f" ({info.cpu_physical_cores} physical core(s), {info.cpu_logical_threads} logical thread(s)"
        cpu_clause += f", {info.cpu_clock_mhz:.0f} MHz)" if info.cpu_clock_mhz is not None else ")"
    memory_clause = f", {info.total_memory_gb} GiB RAM" if info.total_memory_gb is not None else ""
    if info.memory_type or info.memory_speed_mhz is not None:
        memory_clause += " (" + " ".join(filter(None, [
            info.memory_type, f"{info.memory_speed_mhz} MHz" if info.memory_speed_mhz is not None else None,
        ])) + ")"
    gpu_clause = (
        f" A {info.gpu_model} GPU was present but is not utilized by any part of this pipeline "
        f"(reported for reproducibility and any future GPU-accelerated tooling)."
        if info.gpu_model else ""
    )
    return (
        f" Computation was performed on a {info.os_name} {info.architecture} host running {cpu_clause}{memory_clause}, "
        f"Python {info.python_version}.{gpu_clause}"
    )

def render_methods_text(report_data: SampleReportData) -> str:
    """
    Builds the full live methods paragraph for one sample -- concise,
    but covering every real stage the pipeline ran: trimming/QC,
    assembly (including every real rescue attempt), identification,
    annotation, visualization, phylogeny, and validation.

    Args:
        report_data (SampleReportData): The assembled report data.

    Returns:
        str: The methods paragraph, as plain prose.
    """
    sentences = [
        _trimming_sentence(report_data),
        _assembly_sentence(report_data),
        _identification_sentence(report_data),
        _annotation_sentence(report_data),
    ]
    text = " ".join(s for s in sentences if s)
    text += _visualization_sentence(report_data)
    text += _phylogeny_sentence(report_data)
    text += _validation_sentence(report_data)
    text += _environment_sentence(report_data)
    return text
