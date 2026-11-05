"""
run_rescue_mapping.py

Execution layer connecting the rescue trigger to one of two real
recruitment strategies, per the current attempt's own configured
profile:

- "candidate_mapped": the original minimap2 -> samtools sort -> index
  -> fastq chain -- maps reads against the sample's own best candidate
  so far and keeps only the properly-paired reads that mapped. Fast,
  and sufficient when the candidate's real problem is coverage/noise,
  not a genuinely missing region.
- "full_pool": deliberately skips that chain. A read whose mate fails
  to map against an already-INCOMPLETE candidate is exactly a read
  near the gap that candidate is missing -- "candidate_mapped" throws
  that read away precisely when it would matter most. GetOrganelle's
  own internal seed-and-extend recruitment is more sensitive than an
  external mapper working from an incomplete reference, so this mode
  just hands it the complete original read pool again, unfiltered,
  and leaves recruitment to GetOrganelle itself.

Doesn't decide HOW to reassemble (kmer/coverage/seed) -- that's the
reassembly stage's job, reading this step's own output. This step's
only job is: should rescue even happen, and if so, produce the reads
the next attempt should see.
"""

# Imports
import argparse
import gzip
import json
import shutil
from pathlib import Path
from mitopipeline.api.minimap2 import Minimap2Runner
from mitopipeline.api.samtools import SamtoolsRunner
from mitopipeline.logging.logger_factory import make_logger
from mitopipeline.models.sample.sample import Sample
from mitopipeline.utils.rescue.trigger import decide_rescue

VALID_RECRUITMENT_MODES = {"candidate_mapped", "full_pool"}

def parse_args() -> argparse.Namespace:
    """
    Parses the command line arguments for this script.

    Returns:
        argparse.Namespace: The parsed command line arguments.
    """
    parser = argparse.ArgumentParser(description = "Decides whether to rescue one sample's assembly, and if so, runs the minimap2/samtools mapping chain.")

    parser.add_argument("--sample-id", type = str, required = True, help = "The sample id.")
    parser.add_argument("--r1", type = Path, required = True, help = "The trimmed R1 reads to map for rescue.")
    parser.add_argument("--r2", type = Path, required = True, help = "The trimmed R2 reads to map for rescue.")
    parser.add_argument("--reference-fasta", type = Path, required = True,
                        help = "The reference to map against -- the sample's own best assembly candidate so far (contig selection's own choice). Unused when --recruitment is full_pool.")
    parser.add_argument("--recruitment", type = str, default = "candidate_mapped", choices = sorted(VALID_RECRUITMENT_MODES),
                        help = "This attempt's profile.recruitment setting: 'candidate_mapped' (default, today's minimap2/samtools chain) or 'full_pool' (skip it, hand GetOrganelle the complete original read pool again).")
    parser.add_argument("--status", type = str, required = True,
                        help = "The current assembly classification (e.g. from contig_selection.json) -- circular, fragmented, partial, or missing.")
    parser.add_argument("--attempt", type = int, required = True, help = "How many rescue attempts have already been made for this sample (0-indexed).")
    parser.add_argument("--trigger-statuses", type = str, default = "fragmented,partial",
                        help = "Comma-separated statuses that warrant rescue (config.yaml's assembly_rescue.trigger_statuses).")
    parser.add_argument("--max-attempts", type = int, required = True,
                        help = "The configured retry bound -- the number of getorganelle_profiles entries (config.yaml's own documented mechanism; see utils.rescue.trigger).")
    parser.add_argument("--output-dir", type = Path, required = True, help = "The directory to write this attempt's mapping artifacts and decision record to.")
    parser.add_argument("--working-dir", type = Path, required = True, help = "The working directory for the mapping tools.")
    parser.add_argument("--threads", type = int, default = 4, help = "The number of threads to use.")
    parser.add_argument("--log-file", type = Path, required = True, help = "The path to the log file.")
    parser.add_argument("--global-log-file", type = Path, required = True, help = "The path to the global log file.")
    parser.add_argument("--metrics-dir", type = Path, default = None, help = "Where to persist this run's compute metrics as JSON. Omit to skip persistence.")

    return parser.parse_args()

def _write_decision(output_dir: Path, sample_id: str, attempted: bool, reason: str, recruitment: str,
                     rescued_r1: Path | None = None, rescued_r2: Path | None = None) -> None:
    """
    Records this sample's rescue-mapping outcome for the reassembly
    stage to read, without re-deriving the decision itself --
    including WHICH recruitment strategy this attempt actually used,
    as a real structured field (not just folded into the free-text
    reason string), so the report can cite it directly.

    Args:
        output_dir (Path): Where to write the decision record.
        sample_id (str): The sample id.
        attempted (bool): Whether a recruitment strategy was actually run.
        reason (str): A short, human-readable explanation.
        recruitment (str): Which recruitment strategy this attempt used ("candidate_mapped" or "full_pool"), regardless of whether it ultimately ran.
        rescued_r1 (Path | None): The rescued R1 FASTQ, if recruitment ran successfully.
        rescued_r2 (Path | None): The rescued R2 FASTQ, if recruitment ran successfully.

    Returns:
        None
    """
    output_dir.mkdir(parents = True, exist_ok = True)
    path = output_dir / f"{sample_id}.rescue_mapping.json"
    with path.open("w", encoding = "utf-8") as handle:
        json.dump({
            "sample_id": sample_id,
            "attempted": attempted,
            "reason": reason,
            "recruitment": recruitment,
            "rescued_r1": str(rescued_r1) if rescued_r1 else None,
            "rescued_r2": str(rescued_r2) if rescued_r2 else None,
        }, handle, indent = 4)

def _full_pool_passthrough(output_dir: Path, sample_id: str, r1_src: Path, r2_src: Path) -> tuple[Path, Path]:
    """
    Copies the sample's complete original read pool into rescue's own
    expected filenames -- same naming convention the candidate_mapped
    chain's samtools fastq step writes, so the reassembly stage (and
    the Snakemake rule's own file-existence check) never has to know
    which recruitment mode actually produced them. Decompresses
    on the way in when the source is gzipped: the candidate_mapped
    chain's own output is always plain-text FASTQ (samtools fastq's
    real default), and GetOrganelleRunner is only ever confirmed
    against that, not a .gz under a misleading .fastq extension.

    Args:
        output_dir (Path): This attempt's own output directory.
        sample_id (str): The sample id.
        r1_src (Path): The complete original R1 reads (trimmed, not recruitment-narrowed).
        r2_src (Path): The complete original R2 reads (trimmed, not recruitment-narrowed).

    Returns:
        tuple[Path, Path]: (rescued_r1, rescued_r2), both real, already-written files.
    """
    output_dir.mkdir(parents = True, exist_ok = True)
    rescued_r1 = output_dir / f"{sample_id}.rescue.R1.fastq"
    rescued_r2 = output_dir / f"{sample_id}.rescue.R2.fastq"
    for src, dst in ((r1_src, rescued_r1), (r2_src, rescued_r2)):
        if src.suffix == ".gz":
            with gzip.open(src, "rb") as compressed, dst.open("wb") as handle:
                shutil.copyfileobj(compressed, handle)
        else:
            shutil.copy2(src, dst)
    return rescued_r1, rescued_r2

def main() -> int:
    """
    Decides whether to rescue this sample and, if so, runs the real
    mapping chain.

    Returns:
        int: The return code for the process.
    """
    args = parse_args()

    logger = make_logger(
        name = "rescue_mapping",
        log_file_path = args.log_file,
        global_log_file_path = args.global_log_file,
    )
    context = f"{{{args.sample_id}}}"

    config = {
        "trigger_statuses": [s.strip() for s in args.trigger_statuses.split(",") if s.strip()],
        # decide_rescue only needs the COUNT of configured profiles, not
        # their content -- a list of that length stands in for it.
        "getorganelle_profiles": [None] * args.max_attempts,
    }
    decision = decide_rescue(args.status, args.attempt, config, recruitment = args.recruitment)
    logger.info(f"{context} Rescue decision: should_rescue={decision.should_rescue}, "
                f"terminal={decision.terminal}, reason: {decision.reason}")

    if not decision.should_rescue:
        _write_decision(args.output_dir, args.sample_id, attempted = False, reason = decision.reason,
                        recruitment = args.recruitment)
        return 0

    if args.recruitment == "full_pool":
        rescued_r1, rescued_r2 = _full_pool_passthrough(args.output_dir, args.sample_id, args.r1, args.r2)
        logger.info(f"{context} Full-pool recruitment: handing GetOrganelle the complete original "
                    f"read pool again, unfiltered (minimap2/samtools narrowing skipped by design).")
        _write_decision(args.output_dir, args.sample_id, attempted = True, reason = decision.reason,
                        recruitment = args.recruitment, rescued_r1 = rescued_r1, rescued_r2 = rescued_r2)
        return 0

    sample = Sample(sample_id = args.sample_id, r1 = args.r1, r2 = args.r2)

    def _run_stage(runner, step_name: str):
        """Runs one BaseTool stage and treats BOTH failure shapes as a
        chain failure: a validate_outputs() rejection (raises
        SystemExit) AND the underlying command itself exiting nonzero
        (BaseTool.run() returns a CommandResult with success=False
        WITHOUT raising -- confirmed missing here against a real
        failure: samtools index/fastq exiting 1 on an unsupported -@
        flag was silently treated as success, writing a decision
        record that pointed at rescued FASTQ files which were never
        actually created)."""
        result = runner.run()
        if not result.success:
            raise RuntimeError(f"{step_name} exited {result.return_code}")
        return result

    try:
        mapper = Minimap2Runner(working_dir = args.working_dir, output_dir = args.output_dir,
                                sample = sample, reference_fasta = args.reference_fasta,
                                threads = args.threads, logger = logger, metrics_dir = args.metrics_dir)
        _run_stage(mapper, "minimap2 mapping")

        sorter = SamtoolsRunner(working_dir = args.working_dir, output_dir = args.output_dir, sample = sample,
                                operation = "sort", input_path = mapper.output_sam,
                                threads = args.threads, logger = logger, metrics_dir = args.metrics_dir)
        _run_stage(sorter, "samtools sort")

        indexer = SamtoolsRunner(working_dir = args.working_dir, output_dir = args.output_dir, sample = sample,
                                 operation = "index", input_path = sorter.sorted_bam,
                                 threads = args.threads, logger = logger, metrics_dir = args.metrics_dir)
        _run_stage(indexer, "samtools index")

        # Extracting from the ORIGINAL mapping output, not sorter.sorted_bam:
        # confirmed against a real run that samtools fastq needs mates
        # adjacent (name/original order) to pair them -- coordinate
        # sorting (needed for index, above) scatters a pair's two reads
        # apart, and this samtools version silently extracts nothing at
        # all from a coordinate-sorted BAM (0-byte R1/R2, no error) rather
        # than failing loudly. Confirmed directly: the same input through
        # samtools fastq straight from minimap2's SAM correctly yields all
        # 1659 read pairs; through the sorted BAM, it yields zero.
        # Sort+index remain for provenance (a standard mapping artifact
        # to have on record); nothing downstream of rescue actually reads
        # them, only the extracted FASTQ.
        extractor = SamtoolsRunner(working_dir = args.working_dir, output_dir = args.output_dir, sample = sample,
                                   operation = "fastq", input_path = mapper.output_sam,
                                   threads = args.threads, logger = logger, metrics_dir = args.metrics_dir)
        _run_stage(extractor, "samtools fastq")
    except (SystemExit, RuntimeError) as e:
        logger.error(f"{context} Rescue mapping chain failed ({e}).")
        _write_decision(args.output_dir, args.sample_id, attempted = True,
                        reason = f"mapping chain failed: {decision.reason}", recruitment = args.recruitment)
        return 1

    logger.info(f"{context} Rescue mapping chain complete: {extractor.rescued_r1.name}, {extractor.rescued_r2.name}")
    _write_decision(args.output_dir, args.sample_id, attempted = True, reason = decision.reason,
                    recruitment = args.recruitment, rescued_r1 = extractor.rescued_r1, rescued_r2 = extractor.rescued_r2)
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
