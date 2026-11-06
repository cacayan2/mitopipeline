"""
run_rescue_reassembly.py

Execution layer for assembly rescue's reassembly step: takes
the rescued reads the rescue-mapping chain produced and re-runs
GetOrganelle through the SAME GetOrganelleRunner (BaseTool) the
primary attempt uses, but with this attempt's own
getorganelle_profiles entry (config.yaml's documented retry
mechanism -- a different kmer/target_coverage/max_rounds each
attempt) and its own output directory, so repeated attempts never
collide and contig selection can rank every attempt afterward.

Mirrors run_getorganelle_assembly.py's own run/classify/record
lifecycle deliberately -- this IS that same lifecycle, just fed
rescued reads and a named profile's tuning instead of the originally
trimmed reads and no tuning at all.

Seeding: a profile either reuses the primary attempt's own already-
resolved seed (--seed-file, passed straight through -- rescue used to
silently drop it and always ran unseeded, which threw away a real,
already-paid-for signal for no reason) or, when its own config marks
it reseed=true, derives a genuinely DIFFERENT seed by walking further
down the same documented locus fallback chain select_seed() already
uses (--reseed/--exclude-gene) -- never both for the same attempt.
"""

# Imports
import argparse
import json
from logging import Logger
from pathlib import Path
from mitopipeline.api.getorganelle import GetOrganelleRunner, discover_assembly_artifacts, downgrade_unverified_circular
from mitopipeline.logging.logger_factory import make_logger
from mitopipeline.models.sample.sample import Sample
from mitopipeline.utils.seed.models import DEFAULT_SEED_MODE, LOCUS_FALLBACK_ORDER, SEED_MODES
from mitopipeline.utils.seed.seed_selector import select_seed

def parse_args() -> argparse.Namespace:
    """
    Parses the command line arguments for this script.

    Returns:
        argparse.Namespace: The parsed command line arguments.
    """
    parser = argparse.ArgumentParser(description = "Re-runs GetOrganelle on rescued reads with one getorganelle_profiles entry's tuning, and classifies the result.")

    parser.add_argument("--sample-id", type = str, required = True, help = "The sample id.")
    parser.add_argument("--r1", type = Path, required = True, help = "The rescued R1 reads (the rescue-mapping step's own output).")
    parser.add_argument("--r2", type = Path, required = True, help = "The rescued R2 reads (the rescue-mapping step's own output).")
    parser.add_argument("--organelle-type", type = str, default = "animal_mt", help = "The GetOrganelle organelle type.")
    parser.add_argument("--profile-name", type = str, required = True, help = "This attempt's getorganelle_profiles entry name (config.yaml), for logging/provenance only.")
    parser.add_argument("--kmer", type = str, default = None, help = "This profile's -k value (e.g. '21,33,55,85,115').")
    parser.add_argument("--target-coverage", type = float, default = None, help = "This profile's --reduce-reads-for-coverage value.")
    parser.add_argument("--max-rounds", type = int, default = None, help = "This profile's -R value.")
    parser.add_argument("--seed-file", type = Path, default = None,
                        help = "A resolved seed FASTA to reuse as-is for this attempt (typically the primary attempt's own). Ignored when --reseed is set.")
    parser.add_argument("--reseed", action = "store_true",
                        help = "Derive a genuinely different seed for this attempt instead of reusing --seed-file: walks LOCUS_FALLBACK_ORDER starting AFTER --exclude-gene. Requires --scientific-name and --database-dir.")
    parser.add_argument("--exclude-gene", type = str, default = None,
                        help = "The locus the primary attempt's own seed already used (its real seed_selection.json 'gene' field) -- --reseed starts past this one, never repeating it.")
    parser.add_argument("--scientific-name", type = str, default = None, help = "The sample's scientific name. Required when --reseed is set.")
    parser.add_argument("--database-dir", type = Path, default = None, help = "The local MitoFish reference database directory. Required when --reseed is set.")
    parser.add_argument("--seed-mode", type = str, choices = SEED_MODES, default = DEFAULT_SEED_MODE,
                        help = "Seed mode for the reseed (see 04c-seed_sequences' --seed-mode). Only used with --reseed.")
    parser.add_argument("--output-dir", type = Path, required = True, help = "This attempt's OWN output directory -- must differ from the primary attempt's and every other rescue attempt's.")
    parser.add_argument("--working-dir", type = Path, required = True, help = "The working directory.")
    parser.add_argument("--threads", type = int, required = True, help = "The number of threads to use.")
    parser.add_argument("--log-file", type = Path, required = True, help = "The path to the log file.")
    parser.add_argument("--global-log-file", type = Path, required = True, help = "The path to the global log file.")
    parser.add_argument("--metrics-dir", type = Path, default = None, help = "Where to persist this run's compute metrics as JSON. Omit to skip persistence.")

    return parser.parse_args()

def _resolve_reseed(sample_id: str, scientific_name: str | None, database_dir: Path | None,
                     exclude_gene: str | None, output_dir: Path, logger: Logger,
                     seed_mode: str = DEFAULT_SEED_MODE) -> Path | None:
    """
    Derives a genuinely different seed for this rescue attempt by
    walking LOCUS_FALLBACK_ORDER starting just past exclude_gene --
    the same documented-density-ordered chain select_seed() always
    uses, just not restarted from the top (which would almost always
    just re-resolve the exact same locus the primary attempt already
    tried and found insufficient). Writes its own real provenance
    record under output_dir, same shape as 04c-seed_sequences' own
    {sample}.seed_selection.json, so this attempt's real seed choice
    is on record too, not just folded silently into GetOrganelle's
    own tool_options.

    Args:
        sample_id (str): The sample id.
        scientific_name (str | None): The sample's scientific name. None skips reseeding entirely (logged, not an error).
        database_dir (Path | None): The local MitoFish reference database directory. None skips reseeding entirely.
        exclude_gene (str | None): The locus the primary attempt's seed already used, if known. None tries the full chain from the top.
        output_dir (Path): This attempt's own output directory, to write the provenance record under.
        logger (Logger): The logger to use.
        seed_mode (str): The seed mode for this reseed (models.SEED_MODES). "default" yields no FASTA on purpose.

    Returns:
        Path | None: The newly selected seed FASTA, or None if reseeding was skipped or found nothing usable (GetOrganelle's own default database applies, same as any other unseeded run).
    """
    context = f"{{{sample_id}}}"
    # The rule layer passes the manifest's organism verbatim, so a sample
    # with no known scientific name arrives here as "" rather than None.
    if scientific_name is None or not scientific_name.strip() or database_dir is None:
        logger.warning(f"{context} --reseed requested but --scientific-name/--database-dir missing or blank; "
                       "falling back to unseeded for this attempt.")
        return None

    loci = LOCUS_FALLBACK_ORDER
    if exclude_gene in LOCUS_FALLBACK_ORDER:
        loci = LOCUS_FALLBACK_ORDER[LOCUS_FALLBACK_ORDER.index(exclude_gene) + 1:]
    if not loci:
        logger.info(f"{context} --reseed: no locus remains past '{exclude_gene}' in the fallback chain; "
                    "falling back to unseeded for this attempt.")
        return None

    result = select_seed(sample_id = sample_id, scientific_name = scientific_name,
                         database_dir = database_dir, output_dir = output_dir, logger = logger, loci = loci,
                         seed_mode = seed_mode)

    output_dir.mkdir(parents = True, exist_ok = True)
    (output_dir / f"{sample_id}.rescue_seed_selection.json").write_text(json.dumps({
        "sample_id": result.sample_id, "scientific_name": result.scientific_name, "gene": result.gene,
        "accession": result.accession, "source": result.source,
        "fasta_path": str(result.fasta_path) if result.fasta_path else None,
        "fallback_used": result.fallback_used, "fallback_reason": result.fallback_reason,
        "seed_mode": result.seed_mode, "sequence_length": result.sequence_length,
        "excluded_gene": exclude_gene,
    }, indent = 4), encoding = "utf-8")

    if result.fallback_used or result.fasta_path is None:
        logger.info(f"{context} --reseed: no usable candidate past '{exclude_gene}' ({result.fallback_reason}); "
                    "falling back to unseeded for this attempt.")
        return None
    logger.info(f"{context} --reseed: selected locus={result.gene} (primary used '{exclude_gene}').")
    return result.fasta_path

def main() -> int:
    """
    Re-runs GetOrganelle on rescued reads with this attempt's profile
    tuning, and records its classified outcome.

    Returns:
        int: The return code for the process.
    """
    args = parse_args()

    logger = make_logger(
        name = "rescue_reassembly",
        log_file_path = args.log_file,
        global_log_file_path = args.global_log_file,
    )
    context = f"{{{args.sample_id}}}"

    sample = Sample(sample_id = args.sample_id, r1 = args.r1, r2 = args.r2)
    tool_options = {}
    if args.reseed:
        reseeded_fasta = _resolve_reseed(args.sample_id, args.scientific_name, args.database_dir,
                                         args.exclude_gene, args.output_dir, logger, seed_mode = args.seed_mode)
        if reseeded_fasta is not None:
            tool_options["seed_file"] = str(reseeded_fasta)
    elif args.seed_file is not None:
        tool_options["seed_file"] = str(args.seed_file)
    if args.kmer is not None:
        tool_options["kmer"] = args.kmer
    if args.target_coverage is not None:
        tool_options["target_coverage"] = args.target_coverage
    if args.max_rounds is not None:
        tool_options["max_rounds"] = args.max_rounds

    runner = GetOrganelleRunner(
        working_dir = args.working_dir,
        output_dir = args.output_dir,
        sample = sample,
        organelle_type = args.organelle_type,
        tool_options = tool_options,
        threads = args.threads,
        logger = logger,
        metrics_dir = args.metrics_dir,
    )

    logger.info(f"{context} Initiating rescue reassembly (profile '{args.profile_name}': "
                f"kmer={args.kmer}, target_coverage={args.target_coverage}, max_rounds={args.max_rounds}).")

    seed_file_used = tool_options.get("seed_file")

    try:
        result = runner.run()
    except SystemExit as e:
        logger.warning(f"{context} Rescue reassembly (profile '{args.profile_name}') did not produce a usable output (exit {e.code}); classifying as missing.")
        status = discover_assembly_artifacts(args.output_dir, logger)
        _write_status(args.output_dir, args.sample_id, downgrade_unverified_circular(status), success = False,
                     profile_name = args.profile_name, seed_file = seed_file_used)
        runner.cleanup_failed_output()
        _cleanup_rescued_reads(args.r1, args.r2, logger, context)
        return 1
    except Exception as e:
        logger.error(f"{context} Rescue reassembly (profile '{args.profile_name}') failed unexpectedly: {e}.")
        # Unlike the SystemExit branch above, this one previously never
        # called _write_status() -- any unanticipated exception here
        # (not just the specific one this pattern was first noticed for)
        # left no status record at all, just a `.done` marker touched
        # by the .smk rule's own `|| true`. Writing one here closes
        # that gap the same way: discover_assembly_artifacts() globs
        # output_dir regardless of whether it exists (Path.glob on a
        # missing directory yields nothing, not an error), so this is
        # safe even when the exception struck before output_dir was
        # ever created.
        status = discover_assembly_artifacts(args.output_dir, logger)
        _write_status(args.output_dir, args.sample_id, downgrade_unverified_circular(status), success = False,
                     profile_name = args.profile_name, seed_file = seed_file_used)
        runner.cleanup_failed_output()
        _cleanup_rescued_reads(args.r1, args.r2, logger, context)
        return 1

    if not result.success:
        logger.error(f"{context} Rescue reassembly exited with return code {result.return_code} after {result.runtime_seconds:.1f}s.")
        status = discover_assembly_artifacts(args.output_dir, logger)
        _write_status(args.output_dir, args.sample_id, downgrade_unverified_circular(status), success = False,
                     profile_name = args.profile_name, seed_file = seed_file_used)
        runner.cleanup_failed_output()
        _cleanup_rescued_reads(args.r1, args.r2, logger, context)
        return result.return_code or 1

    # Same reasoning as run_getorganelle_assembly.py: postprocess_outputs()
    # (inside run(), above) already standardized/cleaned up GetOrganelle's
    # native filenames, so runner.assembly_status (captured BEFORE that
    # cleanup) is the real classification -- re-deriving it now would
    # always find nothing.
    status = runner.assembly_status
    logger.info(f"{context} Rescue reassembly (profile '{args.profile_name}') complete after "
                f"{result.runtime_seconds:.1f}s; classified as '{status.status}' ({status.candidate_count} native candidate(s)).")
    _write_status(args.output_dir, args.sample_id, status.status, success = True,
                 profile_name = args.profile_name, seed_file = seed_file_used)
    _cleanup_rescued_reads(args.r1, args.r2, logger, context)
    return 0

def _cleanup_rescued_reads(r1: Path, r2: Path, logger: Logger, context: str) -> None:
    """
    Deletes this attempt's own rescued R1/R2 FASTQ (rescue_mapping's
    real output, this script's own --r1/--r2) once GetOrganelle has
    consumed them -- whether it succeeded or not. Nothing downstream
    ever reads these again (they're deliberately NOT declared as
    Snakemake outputs -- see 07-assembly_rescue.smk's own comment),
    and a full_pool attempt's copy is the ENTIRE sample's trimmed read
    pool duplicated verbatim: confirmed as the dominant real cause of
    the disk filling to 100% (2.5TB total across every attempt/sample
    on the live run, dwarfing everything else combined).

    Args:
        r1 (Path): This attempt's rescued R1 FASTQ.
        r2 (Path): This attempt's rescued R2 FASTQ.
        logger (Logger): The logger to use.
        context (str): The log line's sample-id prefix.

    Returns:
        None
    """
    for path in (r1, r2):
        try:
            path.unlink(missing_ok = True)
        except OSError as e:
            logger.warning(f"{context} Could not remove rescued reads {path}: {e}.")

def _write_status(output_dir: Path, sample_id: str, status: str, success: bool, profile_name: str,
                  seed_file: str | None = None) -> None:
    """
    Records this rescue attempt's classification, same shape as the
    primary attempt's (plus which profile produced it and whether/
    which seed it actually used), for contig selection's cross-attempt
    ranking and the report's own rescue narrative to both read without
    re-deriving either.

    Args:
        output_dir (Path): This attempt's own output directory.
        sample_id (str): The sample id.
        status (str): The classification from discover_assembly_artifacts().
        success (bool): Whether the GetOrganelle command itself completed successfully.
        profile_name (str): Which getorganelle_profiles entry produced this attempt.
        seed_file (str | None): The real seed FASTA this attempt actually ran with (reused or freshly re-derived), or None if it ran unseeded.

    Returns:
        None
    """
    output_dir.mkdir(parents = True, exist_ok = True)
    status_path = output_dir / f"{sample_id}.assembly_status.json"
    with status_path.open("w", encoding = "utf-8") as handle:
        json.dump({"sample_id": sample_id, "status": status, "success": success,
                  "profile_name": profile_name, "seed_file": seed_file}, handle, indent = 4)

if __name__ == "__main__":
    raise SystemExit(main())
