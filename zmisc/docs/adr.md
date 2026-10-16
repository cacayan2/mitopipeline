```{=typst}
#set page(margin: (x:0.5in, y:0.5in))
#show link: set text(fill: blue)
#align(center)[

    = Mitogenome Batch Processing Pipeline Architectural Design Record
    == Version 0.2.0
    == Emil Cacayan

]
```

---

# Table of Contents
1. [Executive Summary](#executive-summary)
2. [Architectural Evolution from v0.1.0](#architectural-evolution-from-v010)
3. [Repository Structure and Project Design](#repository-structure-and-project-design)
4. [Python Package Layout](#python-package-layout)
5. [Tool API and the BaseTool Execution Layer](#tool-api-and-the-basetool-execution-layer)
6. [Logging Architecture](#logging-architecture)
7. [Manifest Parsing and the Sample Model](#manifest-parsing-and-the-sample-model)
8. [Pipeline Job and Output Management](#pipeline-job-and-output-management)
9. [Configuration System](#configuration-system)
10. [Mitochondrial Reference Selection (MitoFish)](#mitochondrial-reference-selection-mitofish)
11. [Assembly Rescue](#assembly-rescue)
12. [NCBI/SRA Archival Preparation](#ncbisra-archival-preparation)
13. [Workflow Orchestration Status](#workflow-orchestration-status)
14. [Testing Strategy](#testing-strategy)
15. [Known Implementation Issues](#known-implementation-issues)
16. [Future Work](#future-work)
17. [References](#references)

# Executive Summary
This document is the version 0.2.0 architectural design record (ADR) for the mitogenome batch processing pipeline. It assumes the reader has read the [version 0.1.0 ADR](adr.md) (`zmisc/docs/adr.md` on the `feature/version-0.1.0` branch), which covers the project's scientific motivation, the Fox River/Stuart Lab research context, and the initial proposed architecture. That background is not repeated here. This document instead records what the architecture actually became once implementation began in earnest: an installable, object-oriented Python package (`mitopipeline`) orchestrated by Snakemake, a much finer-grained numbered-stage execution model, a wholly new local-database-backed seed/reference selection subsystem, and scaffolding (config and output directories, but not yet logic) for assembly rescue and NCBI/SRA archival submission.

```{=typst}
#align(center)[
```
[Repository](https://github.com/cacayan2/mitopipeline.git) \
[Project](https://github.com/users/cacayan2/projects/4)
```{=typst}
]
```

As of the last commit reflected in this document (September 17, 2026, `feature/version-0.2.0`), the pipeline is functional through preprocessing (quality control and trimming) and through MitoFish-backed seed/reference lookup in isolation. It is not yet functional end-to-end: several Snakemake rule files referenced by the top-level `Snakefile` do not yet exist, and the assembly, annotation, phylogeny, reporting, and archival stages are represented by output-directory scaffolding and configuration rather than working code. These gaps are documented plainly in [Workflow Orchestration Status](#workflow-orchestration-status), [Assembly Rescue](#assembly-rescue), and [NCBI/SRA Archival Preparation](#ncbisra-archival-preparation) rather than glossed over, consistent with the project's own test-driven, honesty-first development philosophy.

# Architectural Evolution from v0.1.0
## Context and Technology Choice
The v0.1.0 ADR described an *intended* architecture, written largely in future tense ("will be implemented," "is expected to evolve") before significant code existed. v0.2.0 is the first point at which that architecture can be evaluated against real, running code rather than a plan. This section summarizes what changed and offers an honest assessment of what improved, what regressed, and what remains aspirational.

```{=typst}
#rect(fill: luma(245), width: 100%, inset: 15pt)[
```
**Improvements over v0.1.0:**
* **Resource-instrumented execution.** v0.1.0 described `BaseTool` only as handling "subprocess exec, logging, error handling." The implemented `BaseTool.run()` wraps every tool invocation in `/usr/bin/time`, capturing peak memory and user/system CPU seconds into every `CommandResult` — a concrete reproducibility and performance-auditing capability v0.1.0 never specified.
* **Far more granular stage decomposition.** v0.1.0 described five conceptual stages (QC, trimming, assembly, annotation, visualization/reporting). v0.2.0 decomposes these into 22 numbered sub-stages (`01-qc_raw` through `16-table2asn_validation`, several with `a`/`b`/`c` sub-steps), each with its own output directory, log directory, Snakemake rule, and `exec/` wrapper. This is a substantially finer unit of failure isolation and restart than originally proposed.
* **A real reference-selection subsystem.** v0.1.0 did not describe seed/reference selection as an architectural concern at all. v0.2.0 introduces a local, Parquet-backed MitoFish database and a three-stage name-to-accession resolution path (see [Mitochondrial Reference Selection](#mitochondrial-reference-selection-mitofish)) — new, working code with its own unit test.
* **Immutable, metadata-preserving domain objects.** The v0.2.0 `Sample` class retains every nonblank manifest column (not just the fields v0.1.0 anticipated), frozen after construction via `MappingProxyType`. This is a more flexible and defensive design than the dataclass v0.1.0 only sketched in prose.
* **Job metadata persistence.** `PipelineJob` now writes both a live `current_metadata.json` and an appended `historical_metadata.json` per job directory — a concrete, inspectable execution record that v0.1.0 only gestured at ("would store statistics, manage outputs, etc.").
```{=typst}
]
```

```{=typst}
#rect(fill: luma(245), width: 100%, inset: 15pt)[
```
**Regressions and gaps relative to the v0.1.0 plan:**
* **The three-layer documentation commitment has lapsed.** v0.1.0's own "Future Extensibility" section promised an end-user manual, a developer manual, and this ADR as three permanent layers of documentation. As of this writing, `zmisc/docs/` does not exist at all on `feature/version-0.2.0` — it was dropped at some point after the v0.1.0 branch diverged, and none of the three documents carried forward. This document re-establishes the ADR layer; the other two remain absent.
* **The orchestration layer has outrun the implementation.** The `Snakefile` declares `include:` statements for 22 rule files, of which only 6 exist on disk (see [Workflow Orchestration Status](#workflow-orchestration-status)). As committed, the pipeline cannot run past `assembly` without a Snakemake parse error.
* **The most central tool has no API wrapper.** `GetOrganelle` — the tool the entire project exists to run — has an empty `api/getorganelle.py` (0 bytes). The Tool API layer is real and working for `fastp`, `fastqc`, `blast`, `mitos2`, `iqtree`, `megahit`, `minimap2`, and `samtools`, but not yet for the assembly step itself.
* **Several concrete bugs ship in code that looks finished.** See [Known Implementation Issues](#known-implementation-issues) — these were found by reading the modules directly rather than assuming the Technical Runbook's description of them was current, and they are the kind of defect that will not surface until the affected path is actually exercised (manifest parsing and job metadata writing, both on the critical path of every run).
```{=typst}
]
```

## Justifications
Documenting both columns honestly, rather than only the improvements, is the point of an ADR at this stage of the project: the v0.1.0 document was written to justify choices before they had consequences, while this document is written after some of those choices have had consequences, some good and some not yet resolved. Treating the gaps as findings to fix rather than details to omit keeps this document useful as an engineering reference rather than a status report optimized to look complete. Dr. Stuart's planned departure from day-to-day oversight in mid-November 2026 makes an accurate as-built record more valuable than an aspirational one: whoever picks this codebase up next — including a future version of Emil — needs to know what is load-bearing and what is scaffolding.

# Repository Structure and Project Design
## Context and Technology Choice
The repository retains the source-layout structure proposed in v0.1.0 — `src/mitopipeline/` for the installable package, `ctrl/` for Snakemake orchestration, `tests/` for unit and integration tests, `envs/` for per-stage Conda environments — and these four directories are now populated with real, if incomplete, content rather than a plan for them. `zmisc/` (documentation and legacy scripts) is the one piece of the original structure not currently present on this branch; this ADR's addition under `zmisc/docs/adr.md` restores the first of its intended two subdirectories.

## Justifications
Keeping the same top-level layout v0.1.0 proposed, rather than redesigning it mid-project, preserved the boundary between orchestration (`ctrl/`) and reusable implementation (`src/mitopipeline/`) that the original ADR argued for, and nothing observed in the current codebase suggests that boundary has been compromised — Snakemake rules remain thin, and the Python package has absorbed the actual logic. The absence of `zmisc/` is worth treating as a process gap rather than a deliberate simplification: the v0.1.0 ADR specifically argued that this project's complexity (orchestration, OOP design, data modeling, packaging) justified a three-layer documentation commitment, and that argument has not become less true as the codebase has grown larger and more complex than it was at v0.1.0.

# Python Package Layout
## Context and Technology Choice
The installed package now contains six top-level subpackages, organized by responsibility rather than by pipeline stage, matching the principle v0.1.0 argued for:

```{=typst}
#rect(fill: luma(245), width: 100%, inset: 15pt)[
```
* **`launcher/`** — `pipeline_launcher.py`, a single CLI entry point that loads a YAML config, parses or discovers samples, creates or resumes a `PipelineJob`, writes a runtime manifest and runtime config, and invokes Snakemake.
* **`api/`** — one Tool API class per external program, each inheriting `BaseTool`: `fastqc`, `fastp`, `blast`, `mitos2`, `iqtree`, `megahit`, `minimap2`, `samtools`, and an as-yet-empty `getorganelle`.
* **`exec/`** — 22 numbered directories (`01-qc_raw` through `16-table2asn_validation`), each a thin wrapper package intended to translate Snakemake rule parameters into Tool API calls. Most currently contain only an empty `__init__.py`; two (`01-qc_raw`, `02-trimming`) contain real wrapper scripts.
* **`models/`** — `base_tool/`, `command_result/`, `pipeline_job/`, `sample/`, and `stats/` (currently empty placeholder modules for `annotation_stats`, `assembly_stats`, `fastp_stats`, `run_stats`, `sample_summary`).
* **`logging/`** — `logger_factory.py`, a single factory function used by every stage.
* **`utils/`** — `manifest/` (the sample-manifest parser) and `seed/` (the MitoFish reference-lookup subsystem).
```{=typst}
]
```

Two subpackages anticipated by project notes written against this branch — a dedicated `rescue/` package for assembly-rescue logic and a dedicated `archival/` package for NCBI/SRA metadata generation — do not exist. That functionality is represented only by empty `exec/07-assembly_rescue`, `exec/14a-sra_archival`, `exec/14b-biosample_archival`, and `exec/14c-mitogenome_archival` stage stubs and by configuration keys under `stage_options.assembly_post` and `stage_options.archival` in `ctrl/configs/updated_config.yaml`.

## Justifications
Organizing by responsibility (`api/`, `exec/`, `models/`) rather than by pipeline stage, as v0.1.0 argued, has held up: the two real wrapper scripts (`run_fastqc_raw.py`, `run_fastp.py`) are each under 50 lines because they delegate everything substantive to the corresponding Tool API class, and extending the pipeline with `megahit`, `minimap2`, and `samtools` — none of which existed in the v0.1.0 plan — required adding three `api/` modules and no changes to `launcher/`, `models/`, or `logging/`. That is the extensibility the original design argued for, demonstrated rather than asserted.

The absence of `rescue/` and `archival/` as real packages is not itself a design flaw — the `exec/` stage-stub convention is a reasonable place to eventually put that logic — but it does mean two of the project's eventual scientific deliverables (recovering fragmented assemblies, submitting to NCBI) currently exist only as directory names and config schema. Readers of this document, and of the Technical Runbook that preceded it, should not infer working logic from the presence of these directories.

# Tool API and the BaseTool Execution Layer
## Context and Technology Choice
Every Tool API class inherits `BaseTool` (`models/base_tool/base_tool.py`), an abstract base class following the template-method pattern. Concrete subclasses implement three abstract methods — `validate_inputs()`, `build_command()`, `validate_outputs()` — and an optional `postprocess_outputs()` hook; `BaseTool.run()` is concrete and orchestrates all four around the actual subprocess execution:

1. Logs a sample-scoped context label, then calls `validate_inputs()`.
2. Calls `build_command()` to construct the argv list.
3. Wraps the command in `/usr/bin/time -f '{"user_cpu_seconds": %U, "system_cpu_seconds": %S, "peak_memory_mb": %M}'`, writing the metrics to a temporary file.
4. Executes via `subprocess.run(..., capture_output=True, text=True)`, timing wall-clock runtime with `time.perf_counter()`.
5. Reads the GNU time metrics back, converts peak memory from KB to MB, and deletes the temp file in a `finally` block.
6. Constructs a `CommandResult` carrying the command, return code, stdout/stderr, timestamps, runtime, and resource metrics.
7. On failure (nonzero return code), logs the error and returns immediately. On success, calls `postprocess_outputs()` then `validate_outputs()` before returning.

Eight tools currently have concrete Tool API classes: `fastqc`, `fastp`, `blast`, `mitos2`, `iqtree`, `megahit`, `minimap2`, and `samtools`. `api/getorganelle.py` exists as a file but is empty.

## Justifications
Wrapping every external call in GNU `time` was not part of the v0.1.0 design and is a genuine improvement: because it is implemented once in `BaseTool.run()` rather than per-tool, every current and future Tool API gets runtime, memory, and CPU accounting for free, which directly supports the project's stated goal (per the v0.1.0 ADR) of auditable, reproducible batch execution — assembly runtime and memory use vary enormously across ~140 species, and this gives the eventual batch report a basis for flagging outliers without any additional per-tool work.

The template-method structure (abstract `validate_inputs`/`build_command`/`validate_outputs`, concrete `run`) enforces a consistent contract across tools with genuinely different command-line interfaces, which is the same justification v0.1.0 gave for the `API` object in principle — the difference is that this is now a real `ABC` with `@abstractmethod` decorators enforcing it at instantiation time, not a documented convention developers could silently deviate from.

That `getorganelle.py` is empty is the most consequential open item in the Tool API layer: GetOrganelle is the only assembly tool the pipeline currently targets, and every downstream stage (rescue, BLAST validation, annotation, phylogeny) depends on its output. Until this wrapper exists, the pipeline has quality control and reference selection but no assembly, which is consistent with the Development Status note that the pipeline is "functional through preprocessing and seed/reference lookup" and no further.

# Logging Architecture
## Context and Technology Choice
`logging/logger_factory.py` implements the two-stream design v0.1.0 proposed (concise console, detailed file) as a single `make_logger()` factory used by every stage. Each logger writes to three destinations: the console (filtered via a custom `ExcludeToolOutputFilter` that drops records tagged `event_type="tool_output"`, keyed to a `TOOL_OUTPUT_EVENT` constant), a component-specific log file, and a job-wide global log file. Loggers are keyed by both a logical name and the resolved log file path specifically to prevent two samples using the same tool name from silently sharing a file handler within the same process. `PipelineJob.create_output_directory()` calls `make_logger()` once per job to produce the job-level logger; downstream stages are expected to call it again per component.

## Justifications
Filtering tool stdout/stderr out of the console stream while still capturing it in file handlers delivers exactly the "concise console, detailed file" behavior v0.1.0 specified, implemented as a `logging.Filter` rather than ad hoc string matching, which keeps the filtering logic testable and reusable across every logger the factory produces. Keying loggers by name-plus-filepath rather than name alone is a detail v0.1.0 did not anticipate but addresses a real correctness risk in a pipeline that runs the same tool (e.g., `fastqc`) across ~140 samples: without it, two samples' loggers could resolve to the same Python `logging.Logger` singleton and interleave or misattribute output. Writing to both a component log and a global job log gives a developer the choice of debugging a single stage in isolation or reading one continuous narrative of a run, without requiring two separate logging configurations.

# Manifest Parsing and the Sample Model
## Context and Technology Choice
`utils/manifest/manifest_parser.py` implements the manifest-to-`Sample` pipeline v0.1.0 proposed, using `pandas` to read a tab-separated manifest with `dtype=str` and `keep_default_na=False` (so that blank cells remain empty strings rather than becoming `NaN`). It validates that `sample_id`, `r1`, and `r2` columns are present and non-duplicated, resolves `r1`/`r2` paths relative to either the manifest's own directory or an optional separate input directory, and constructs one `Sample` object per row.

`models/sample/sample.py` defines `Sample` as a frozen `dataclass`. Unlike a schema that enumerates only the fields the pipeline currently uses, every nonblank manifest column is retained in an immutable `metadata` mapping (`MappingProxyType`) in `__post_init__`, so a manifest can carry additional project- or sample-specific columns (e.g., a collection date, or an `organism` override for the scientific name) without requiring a code change to `Sample` itself. `scientific_name` is a computed property: it prefers an explicit `organism` metadata value, falling back to `"{genus} {species} {subspecies}"` constructed from the available parts.

## Justifications
Retaining every manifest column rather than a fixed field list is a direct, positive response to the flexibility v0.1.0 argued for ("supporting optional metadata fields while enforcing a minimal required schema") — archival and reporting stages that need a field never anticipated by this ADR can read it from `sample.get_metadata(...)` without a `Sample` schema migration. Freezing both the dataclass and the metadata mapping after construction (`frozen=True` plus `MappingProxyType`) prevents a later pipeline stage from silently mutating sample data that an earlier stage already validated — a defensive design choice v0.1.0 did not specify but which fits its own stated goal of auditability: a `Sample` object's contents at stage 16 are guaranteed identical to its contents at stage 1.

This path currently cannot execute successfully: see [Known Implementation Issues](#known-implementation-issues) for a keyword-argument typo in `manifest_parser.py` that raises a `TypeError` on every call to construct a `Sample`.

# Pipeline Job and Output Management
## Context and Technology Choice
`models/pipeline_job/pipeline_job.py` implements the job-oriented execution model v0.1.0 proposed. `PipelineJob.__init__` either generates a new job ID (`YYYYMMDD_HMM_<uuid4>`) or resumes an existing one if a `job_id` is passed in, then creates 22 numbered output subdirectories (`01-qc_raw` through `16-table2asn_v`) and a parallel set of per-stage log directories beneath `job_dir/logs/`. `mark_created()`, `mark_failed()`, and `mark_completed()` each write a metadata snapshot via `write_metadata()`: a live `current_metadata.json` is overwritten on every call, and the same payload is appended to a cumulative `historical_metadata.json`.

## Justifications
Generating output and log directories for every defined stage at job creation time — rather than lazily, per stage, as it runs — gives a developer or lab member a complete, predictable directory skeleton to inspect immediately after launching a job, which is useful during the current phase of the project where most of those directories will sit empty until the corresponding stage is implemented: an empty `07-assembly_rescue/` directory is a visible, honest signal of unfinished work rather than a missing path that would otherwise need explaining. Maintaining both a live and a historical metadata file lets a developer answer "what is this job doing right now" and "what has this job directory been used for across every run" from two files rather than reconstructing history from log greps, directly satisfying the auditability goal v0.1.0 set for the `PipelineJob` object.

Two defects currently prevent this metadata writing from functioning as described — see [Known Implementation Issues](#known-implementation-issues).

# Configuration System
## Context and Technology Choice
`ctrl/configs/updated_config.yaml` is the current production configuration, structured as three top-level blocks: `job_configuration` (output root, manifest path, input directory, job ID, input validation toggle), `stages` (a boolean per numbered stage, including nested toggles for `assembly_post.contig_selection`/`assembly_post.assembly_rescue`), and `stage_options` (per-tool parameter blocks for `fastqc`, `fastp`, `getorganelle`, the assembly-assessment and assembly-rescue heuristics, `mitos2`, BLAST, MAFFT, IQ-TREE, Pandoc/Typst report export, and SRA/BioSample/table2asn archival). `tests/fixtures/config/` retains the multiple-profile approach v0.1.0 described (e.g., `config_preprocessing_only.yaml`, `config_qc_raw_only.yaml`, `config.manifest_free.yaml`) for scoped test runs. `pipeline_launcher.py` loads whichever config is passed via `--config` and fails fast (`sys.exit(1)`, logged) if the file is missing or empty.

The `assembly_post.assembly_rescue` block is worth noting in detail because it is more fully specified than any code that currently consumes it: it defines trigger statuses (`fragmented`, `partial`), length and gap thresholds, terminal-overlap parameters, and boolean toggles for four candidate recovery strategies (`run_read_mapping`, `run_graph_resolution`, `run_reference_scaffolding`, `run_megahit`), plus two named alternate GetOrganelle re-run profiles (`reduced_coverage`, `alternate_kmers`) with their own k-mer and coverage-target overrides.

## Justifications
A single structured YAML file with a `stages`/`stage_options` split lets a developer disable everything past preprocessing (as the project currently must, given the gaps described in this document) without touching Python code, which is exactly the decoupling v0.1.0 argued for. That the `assembly_rescue` configuration schema is this detailed — four named strategies, two re-run profiles — while no `rescue/` implementation exists to read it is worth reading as a signal of design-first intent rather than as a documentation error: it suggests the rescue approach has been thought through (likely informed by the fragmentation issues the Technical Runbook records for the `Fundulus` test sample) ahead of implementation, which is consistent with the project's TDD-roadmap-style issue tracker. Whether that config schema survives contact with the actual implementation unchanged is something this ADR cannot verify and a future revision should confirm once `exec/07-assembly_rescue` has real code behind it.

# Mitochondrial Reference Selection (MitoFish)
## Context and Technology Choice
`utils/seed/` is new in v0.2.0; it has no counterpart in the v0.1.0 ADR. `MitoFishReferenceLookup` (`mitofish.py`) resolves a scientific name to candidate mitochondrial reference sequences against a local, Parquet-backed MitoFish database (`mitofish_database.py` handles setup; four tables are queried: `taxonid_name`, `speciesid_lineageid`, `taxonid_speciesid`, `seq_taxonid`, `seq_annotation`). `lookup(scientific_name, gene="12S_rRNA")` performs three sequential resolution steps, each logged and each short-circuiting to an empty list with a warning if it finds nothing:

1. `_lookup_taxon_ids()` — exact scientific-name match against `taxonid_name`, joined through `speciesid_lineageid` and `taxonid_speciesid` to NCBI taxonomy IDs. Also fails loudly (`sys.exit(1)`) if any of the three required Parquet tables is missing, or if a given accession resolves to more than one taxonomy ID (treated as a data-integrity error, not a warning).
2. `_lookup_accessions()` — taxonomy IDs to mitochondrial sequence accessions via `seq_taxonid`.
3. `_build_candidates()` — accessions filtered to the requested `gene`, producing one `SeedReferenceCandidate` (accession, scientific name, NCBI taxon ID, gene, sequence length) per matching annotation row, de-duplicated while preserving order.

`SeedReferenceCandidate` (`models.py`) is a frozen, slotted dataclass. A real unit test exists at `tests/unit/04b-mitofish_config/test_mitofish_reference_lookup.py`.

## Justifications
Resolving through MitoFish's own internal ID hierarchy (name → lineage → species → taxon → accession → annotation) rather than querying NCBI directly at runtime keeps reference selection fast and offline-capable for a batch of ~140 species, which matters for a pipeline intended to run unattended against a large sample set. Treating an accession that maps to more than one taxonomy ID as a hard failure rather than a silent pick is a conservative, defensible choice for a pipeline whose downstream BLAST and phylogenetic validation steps depend on knowing exactly which organism a reference came from — an ambiguous mapping here would be a data-quality problem worth surfacing immediately rather than resolving arbitrarily.

This subsystem is also the clearest example in the current codebase of a feature that is real but not yet integrated: `lookup()` has no callers anywhere else in `src/`. It is not wired into any `exec/` stage or Snakemake rule, it only resolves one gene per call rather than cascading through a priority-ordered list of biomarkers (e.g., 12S → 16S → COI → cytb) when the first yields nothing, and there is no fallback to GetOrganelle's own default `animal_mt` reference when MitoFish has no entry for an organism at all. Each of these is tracked as a distinct, already-filed open issue rather than an unscoped gap (see [Future Work](#future-work)), which is a reasonable way to sequence this kind of work — get the core resolution logic correct and tested in isolation before building the fallback and integration layers on top of it — but it does mean the pipeline today cannot yet select a seed sequence end-to-end for any sample.

# Assembly Rescue
## Context and Technology Choice
Assembly rescue is represented in the current codebase only as an empty `exec/07-assembly_rescue/__init__.py` stub and the detailed `stage_options.assembly_post.assembly_rescue` configuration block described in [Configuration System](#configuration-system). No rescue logic, retry orchestration, or `rescue/` package exists anywhere in `src/mitopipeline` as of this writing.

## Justifications
This section is included — short as it is — because the Technical Runbook maintained alongside this project described a working `rescue/assembly_rescue.py` module, and that description does not match the source tree as directly inspected while writing this ADR. Recording the actual state here, rather than repeating the prior description, matters for the same reason the v0.1.0 ADR argued for ADRs in general: a document that states what the code does should be corrected the moment it is found to be wrong, not left to accumulate as a second, inconsistent source of truth alongside the code itself.

# NCBI/SRA Archival Preparation
## Context and Technology Choice
As with assembly rescue, archival preparation exists in v0.2.0 only as configuration and directory scaffolding: `stage_options.archival.{sra, biosample, table2asn}` in `updated_config.yaml` defines the metadata fields v0.1.0 anticipated (library strategy/source/selection/layout, platform, BioSample isolation source/collection date/geo-location, `table2asn` submission template path), and `exec/14a-sra_archival`, `exec/14b-biosample_archival`, `exec/14c-mitogenome_archival`, `exec/15-archival_validation`, and `exec/16-table2asn_validation` exist as empty stage stubs with corresponding (also currently empty) `tests/unit/` directories. No `archival/` package or metadata-generation logic exists in `src/mitopipeline`.

## Justifications
That the configuration schema here is as detailed as it is — distinct blocks for SRA, BioSample, and `table2asn` validation, matching the three-part submission v0.1.0 described as a stretch goal — suggests the archival requirements have been researched and specified even though implementation has not started, which is a reasonable allocation of effort this early: NCBI's submission requirements are fixed external constraints worth getting right in the schema before writing code against them. As with assembly rescue, this document records that reality plainly rather than inferring progress from the presence of well-named directories and config keys.

# Workflow Orchestration Status
## Context and Technology Choice
`ctrl/Snakefile` declares `include:` statements for 22 rule files under `ctrl/rules/`, one per numbered stage. As committed, only 6 of those 22 files exist: `01-qc_raw.smk`, `02-trimming.smk`, `03-qc_trimmed.smk`, `04a-getorganelle_config.smk`, `04b-assembly.smk`, and `05a-blast_database_setup.smk`. The remaining 16 — including `04c-seed_sequences.smk`/`04d-assembly.smk` (note: the `exec/` package and `PipelineJob` output directories were renumbered to `04c-seed_sequences`/`04d-assembly`, but the one existing assembly rule file is still named `04b-assembly.smk`, and the one existing BLAST rule file is `05a-blast_database_setup.smk` rather than the `05a-blast_config`/`05b-blast`/`05c-blast_hits` split the rest of the codebase now uses) — do not exist. Snakemake resolves `include:` paths at parse time, so running this `Snakefile` as committed fails before any rule executes, once a target past `04b-getorganelle_config` is requested.

## Justifications
This is recorded as a status finding rather than folded into a "future work" list because it is a correctness issue in code that is already committed and already wired into the `Snakefile`'s `include:` list, not a feature that has simply not been started yet. It is also informative on its own: the stage renumbering visible in `PipelineJob` and `exec/` (which use `04c-seed_sequences`/`04d-assembly`) having outrun the `ctrl/rules/` directory (which still uses the older `04b-assembly`/`05a-blast_database_setup` names) indicates the numbering scheme itself was revised mid-development, and that revision was only partially propagated. Anyone resuming work on orchestration should treat the `exec/`/`PipelineJob` numbering as the current source of truth and reconcile `ctrl/rules/` and the `Snakefile`'s `include:` list to match it, rather than assuming the two are already consistent.

# Testing Strategy
## Context and Technology Choice
The unit/integration/fixture structure v0.1.0 proposed is in place: `tests/unit/`, `tests/integration/`, and the SRA-accession-based fixture strategy (the same six species listed in the v0.1.0 ADR — lemon shark, common carp, blue mahseer, Korean loach, and two sea chub species — regenerated deterministically rather than committed as binary FASTQ). `tests/unit/` now mirrors the 22-stage numbering (`00-job` through `16-table2asn_validation`), one directory per stage. As of this writing, only `04b-mitofish_config/test_mitofish_reference_lookup.py` contains a real test; every other stage directory contains only a `.gitkeep` placeholder.

## Justifications
Mirroring the pipeline's own stage numbering in `tests/unit/` makes it immediately legible which stages have test coverage and which do not — a developer or reviewer does not need to cross-reference a separate test plan to see that 21 of 22 stage directories are still placeholders, because the repository structure states it directly. This is consistent with the issue tracker's own pre-planned TDD-roadmap style noted in the project's Development Status record: tests and their directories are scaffolded ahead of implementation, which is a defensible ordering for a short-timeline project where the next engineer (or the same engineer returning to a stage) should write the test before writing the stage logic rather than after.

# Known Implementation Issues
The following were found by reading the modules directly while compiling this document, not by inference from configuration or prior notes. They are recorded here, with file and line references, because each sits on the critical path of a run and none is reflected in any existing project documentation.

```{=typst}
#rect(fill: luma(245), width: 100%, inset: 15pt)[
```
1. **`utils/manifest/manifest_parser.py`, line ~120:** `Sample(..., speices = _optional_value(metadata, "species"), ...)` — the keyword argument is misspelled `speices`. `Sample` is a plain `@dataclass` with a field named `species`, not `speices`, so this raises `TypeError: __init__() got an unexpected keyword argument 'speices'` on every call. As written, manifest parsing cannot currently complete for any manifest.
2. **`models/pipeline_job/pipeline_job.py`, lines ~201 and ~224:** `json.dump(data, handle, ident = 4)` (both the live-metadata write and the historical-metadata write) — the keyword is `indent`, not `ident`. `json.dump()` does not accept an `ident` keyword and raises `TypeError` on every call, meaning `write_metadata()` — and therefore `mark_created()`, `mark_failed()`, and `mark_completed()` — currently cannot succeed.
3. **`models/command_result/command_result.py`, lines ~76–84:** the `cpu_time_seconds` property is defined with `@property` *inside* `__init__()`, after the instance attributes are assigned. A `def` (even a decorated one) written inside `__init__` creates a local function object that is discarded when `__init__` returns; it is never attached to the class or any instance. `cpu_time_seconds` is consequently not accessible on any `CommandResult` instance despite appearing in the class docstring's attribute list.
```{=typst}
]
```

None of these three is architecturally significant on its own — each is a small, mechanical fix — but together they mean that two of the paths this document describes as "implemented" (manifest parsing, job metadata persistence) cannot currently run to completion, and a result field documented as available (`cpu_time_seconds`) is not. They are listed here rather than silently fixed in the course of writing this document because an ADR's job is to describe the architecture as it exists, and because Emil should decide how and when to land the fixes rather than have them arrive bundled, unannounced, inside a documentation commit.

# Future Work
Carried forward from the open-issue backlog and from the gaps identified in this document, roughly in dependency order:

```{=typst}
#rect(fill: luma(245), width: 100%, inset: 15pt)[
```
1. Fix the three defects in [Known Implementation Issues](#known-implementation-issues) — manifest parsing and job metadata persistence are both currently non-functional.
2. Implement `api/getorganelle.py` and its `exec/04d-assembly` wrapper; wire the MitoFish seed result into GetOrganelle's configuration.
3. Extend `MitoFishReferenceLookup` with a priority-ordered biomarker cascade (12S → 16S → COI → cytb) and a fallback to GetOrganelle's default `animal_mt` reference when MitoFish has no entry.
4. Create the 16 missing `ctrl/rules/*.smk` files and reconcile the `04b`/`05a` naming drift described in [Workflow Orchestration Status](#workflow-orchestration-status) so the `Snakefile` parses cleanly end-to-end.
5. Implement assembly rescue logic against the already-specified `assembly_post.assembly_rescue` config schema.
6. Implement BLAST validation, MITOS2 annotation, and MAFFT/IQ-TREE phylogeny wrappers and wire them into the DAG.
7. Produce one verified single-sample PDF report before attempting batch orchestration.
8. Implement the SRA/BioSample/`table2asn` archival stages against the already-specified config schema.
9. Restore the second documentation layer (`zmisc/docs/developer-manual.md`) that the v0.1.0 ADR committed to.
```{=typst}
]
```

# References
This document is an engineering status record rather than a literature review; the scientific background, motivation, and citations remain in the [version 0.1.0 ADR](https://github.com/cacayan2/mitopipeline/blob/feature/version-0.1.0/zmisc/docs/adr.md). Tools referenced in this document:

[Snakemake](https://snakemake.readthedocs.io/) · [Conda](https://docs.conda.io/) · [pytest](https://docs.pytest.org/en/stable/) · [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) · [fastp](https://github.com/opengene/fastp) · [GetOrganelle](https://github.com/kinggerm/getorganelle) · [MITOS2](https://gitlab.com/Bernt/MITOS) · [IQ-TREE](http://www.iqtree.org/) · [MitoFish](https://mitofish.aori.u-tokyo.ac.jp/) · [NCBI BLAST](https://blast.ncbi.nlm.nih.gov/) · [NCBI Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra)

---
