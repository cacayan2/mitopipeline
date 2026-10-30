"""
run_getorganelle_config.py

Execution layer that resolves one sample's GetOrganelle assembly
configuration deterministically -- organelle type, seed/no-seed input,
and the output directory the assembly stage will write to -- and
records it as a JSON manifest the assembly exec reads next.

Doesn't run GetOrganelle itself; that's the assembly stage's own job.
This step exists so the seed/no-seed decision and the output paths are
pinned down and logged BEFORE a (slow) assembly attempt starts, not
discovered as a side effect of it.
"""

# Imports
import argparse
import json
from pathlib import Path
from mitopipeline.api.getorganelle import VALID_ORGANELLE_TYPES
from mitopipeline.logging.logger_factory import make_logger

def parse_args() -> argparse.Namespace:
    """
    Parses the command line arguments for this script.

    Returns:
        argparse.Namespace: The parsed command line arguments.
    """
    parser = argparse.ArgumentParser(description = "Resolves one sample's GetOrganelle assembly configuration.")

    parser.add_argument("--sample-id", type = str, required = True, help = "The sample id.")
    parser.add_argument("--organelle-type", type = str, required = True, help = "The GetOrganelle organelle type to assemble (e.g. animal_mt).")
    parser.add_argument("--seed-file", type = Path, default = None, help = "A resolved seed FASTA for this sample (omit for an unseeded, fallback-database run).")
    parser.add_argument("--assembly-output-dir", type = Path, required = True, help = "The directory the assembly stage will write this sample's output to.")
    parser.add_argument("--output-dir", type = Path, required = True, help = "The directory to write this config step's own output (the resolved-config JSON and .done marker) to.")
    parser.add_argument("--log-file", type = Path, required = True, help = "The path to the log file.")
    parser.add_argument("--global-log-file", type = Path, required = True, help = "The path to the global log file.")

    return parser.parse_args()

def resolve_config(sample_id: str, organelle_type: str, seed_file: Path | None,
                    assembly_output_dir: Path) -> dict:
    """
    Resolves one sample's deterministic GetOrganelle configuration.
    Pure function (no I/O beyond what's already in its arguments) --
    the same inputs always resolve to the same configuration, with no
    hidden state or clock/randomness involved.

    Args:
        sample_id (str): The sample id.
        organelle_type (str): The requested GetOrganelle organelle type.
        seed_file (Path | None): A resolved seed FASTA, or None for an unseeded run.
        assembly_output_dir (Path): Where the assembly stage will write this sample's output.

    Returns:
        dict: The resolved configuration, JSON-serializable.

    Raises:
        ValueError: If organelle_type isn't one GetOrganelle actually supports.
    """
    if organelle_type not in VALID_ORGANELLE_TYPES:
        raise ValueError(f"Unsupported organelle type: {organelle_type!r}. Must be one of {sorted(VALID_ORGANELLE_TYPES)}.")

    return {
        "sample_id": sample_id,
        "organelle_type": organelle_type,
        "seeded": seed_file is not None,
        "seed_file": str(seed_file) if seed_file is not None else None,
        "assembly_output_dir": str(assembly_output_dir),
    }

def main() -> int:
    """
    Resolves and records one sample's GetOrganelle configuration.

    Returns:
        int: The return code for the process.
    """
    args = parse_args()

    logger = make_logger(
        name = "getorganelle_config",
        log_file_path = args.log_file,
        global_log_file_path = args.global_log_file,
    )
    context = f"{{{args.sample_id}}}"

    try:
        config = resolve_config(
            sample_id = args.sample_id,
            organelle_type = args.organelle_type,
            seed_file = args.seed_file,
            assembly_output_dir = args.assembly_output_dir,
        )
    except ValueError as e:
        logger.error(f"{context} Could not resolve GetOrganelle configuration: {e}")
        return 1

    args.output_dir.mkdir(parents = True, exist_ok = True)
    config_path = args.output_dir / f"{args.sample_id}.getorganelle_config.json"
    with config_path.open("w", encoding = "utf-8") as handle:
        json.dump(config, handle, indent = 4)

    seed_desc = config["seed_file"] if config["seeded"] else "none (unseeded, fallback database)"
    logger.info(f"{context} Resolved GetOrganelle config: organelle_type={config['organelle_type']}, seed={seed_desc}. Written to {config_path}.")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
