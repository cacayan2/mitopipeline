"""
run_rescue_mapping.py

Execution layer connecting the rescue trigger (#405) to the proven
minimap2 -> samtools sort -> index -> fastq chain (#423/#424/#425):
decides whether this sample's current assembly classification actually
warrants a rescue attempt, and if so, runs the real chain and records
the rescued paired reads #406's reassembly stage will consume next.

Doesn't decide HOW to reassemble (which getorganelle_profiles entry,
what k-mer/coverage) -- that's #406's job, reading this step's own
output. This step's only job is: should rescue even happen, and if so,
produce the rescued reads.
"""

# Imports
import argparse
import json
from pathlib import Path
from mitopipeline.api.minimap2 import Minimap2Runner
from mitopipeline.api.samtools import SamtoolsRunner
from mitopipeline.logging.logger_factory import make_logger
from mitopipeline.models.sample.sample import Sample
from mitopipeline.utils.rescue.trigger import decide_rescue

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
                        help = "The reference to map against -- the sample's own best assembly candidate so far (#402's selection).")
    parser.add_argument("--status", type = str, required = True,
                        help = "The current assembly classification (e.g. from #402's contig_selection.json) -- circular, fragmented, partial, or missing.")
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

def _write_decision(output_dir: Path, sample_id: str, attempted: bool, reason: str,
                     rescued_r1: Path | None = None, rescued_r2: Path | None = None) -> None:
    """
    Records this sample's rescue-mapping outcome for #406 to read,
    without re-deriving the decision itself.

    Args:
        output_dir (Path): Where to write the decision record.
        sample_id (str): The sample id.
        attempted (bool): Whether the mapping chain was actually run.
        reason (str): A short, human-readable explanation.
        rescued_r1 (Path | None): The rescued R1 FASTQ, if the chain ran successfully.
        rescued_r2 (Path | None): The rescued R2 FASTQ, if the chain ran successfully.

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
            "rescued_r1": str(rescued_r1) if rescued_r1 else None,
            "rescued_r2": str(rescued_r2) if rescued_r2 else None,
        }, handle, indent = 4)

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
    decision = decide_rescue(args.status, args.attempt, config)
    logger.info(f"{context} Rescue decision: should_rescue={decision.should_rescue}, "
                f"terminal={decision.terminal}, reason: {decision.reason}")

    if not decision.should_rescue:
        _write_decision(args.output_dir, args.sample_id, attempted = False, reason = decision.reason)
        return 0

    sample = Sample(sample_id = args.sample_id, r1 = args.r1, r2 = args.r2)

    try:
        mapper = Minimap2Runner(working_dir = args.working_dir, output_dir = args.output_dir,
                                sample = sample, reference_fasta = args.reference_fasta,
                                threads = args.threads, logger = logger, metrics_dir = args.metrics_dir)
        mapper.run()

        sorter = SamtoolsRunner(working_dir = args.working_dir, output_dir = args.output_dir, sample = sample,
                                operation = "sort", input_path = mapper.output_sam,
                                threads = args.threads, logger = logger, metrics_dir = args.metrics_dir)
        sorter.run()

        indexer = SamtoolsRunner(working_dir = args.working_dir, output_dir = args.output_dir, sample = sample,
                                 operation = "index", input_path = sorter.sorted_bam,
                                 threads = args.threads, logger = logger, metrics_dir = args.metrics_dir)
        indexer.run()

        extractor = SamtoolsRunner(working_dir = args.working_dir, output_dir = args.output_dir, sample = sample,
                                   operation = "fastq", input_path = sorter.sorted_bam,
                                   threads = args.threads, logger = logger, metrics_dir = args.metrics_dir)
        extractor.run()
    except SystemExit as e:
        logger.error(f"{context} Rescue mapping chain failed (exit {e.code}).")
        _write_decision(args.output_dir, args.sample_id, attempted = True,
                        reason = f"mapping chain failed: {decision.reason}")
        return 1

    logger.info(f"{context} Rescue mapping chain complete: {extractor.rescued_r1.name}, {extractor.rescued_r2.name}")
    _write_decision(args.output_dir, args.sample_id, attempted = True, reason = decision.reason,
                    rescued_r1 = extractor.rescued_r1, rescued_r2 = extractor.rescued_r2)
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
