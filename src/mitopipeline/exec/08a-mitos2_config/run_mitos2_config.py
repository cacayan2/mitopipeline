"""
run_mitos2_config.py

Execution layer that resolves one sample's MITOS2 configuration
deterministically -- the selected FASTA to annotate, genetic code,
taxonomy reference (refdir/refseqver), and the output directory the
annotation exec stage will write to -- and records it as a
JSON manifest that stage reads next.

Doesn't run MITOS2 itself; that's the annotation stage's job. This
step exists so an invalid option (an unsupported genetic code, a
missing taxonomy reference directory) is caught and logged before a
real MITOS2 invocation starts.
"""

# Imports
import argparse
import json
from pathlib import Path
from mitopipeline.logging.logger_factory import make_logger

def parse_args() -> argparse.Namespace:
    """
    Parses the command line arguments for this script.

    Returns:
        argparse.Namespace: The parsed command line arguments.
    """
    parser = argparse.ArgumentParser(description = "Resolves one sample's MITOS2 configuration.")

    parser.add_argument("--sample-id", type = str, required = True, help = "The sample id.")
    parser.add_argument("--input-fasta", type = Path, required = True, help = "The selected mitogenome FASTA to annotate (contig selection's own output).")
    parser.add_argument("--genetic-code", type = int, required = True, help = "The genetic code to use.")
    parser.add_argument("--refdir", type = Path, required = True, help = "The base directory containing the MITOS2 reference data.")
    parser.add_argument("--refseqver", type = str, required = True, help = "The reference data version, a subdirectory of refdir.")
    parser.add_argument("--circular", action = "store_true", default = True, help = "Treat the sequence as circular (MITOS2's own default).")
    parser.add_argument("--linear", dest = "circular", action = "store_false", help = "Treat the sequence as linear.")
    parser.add_argument("--noplots", action = "store_true", default = False, help = "Do not create MITOS2's own plots.")
    parser.add_argument("--best", action = "store_true", default = False, help = "Annotate only the best copy of each feature.")
    parser.add_argument("--ncbicode", action = "store_true", default = False, help = "Use NCBI start/stop codons instead of MITOS2's learned ones.")
    parser.add_argument("--zip-output", action = "store_true", default = False, help = "Create a zip of MITOS2's own output.")
    parser.add_argument("--annotation-output-dir", type = Path, required = True, help = "The directory the annotation stage will write this sample's output to.")
    parser.add_argument("--output-dir", type = Path, required = True, help = "The directory to write this config step's own output to.")
    parser.add_argument("--log-file", type = Path, required = True, help = "The path to the log file.")
    parser.add_argument("--global-log-file", type = Path, required = True, help = "The path to the global log file.")

    return parser.parse_args()

def resolve_config(sample_id: str, input_fasta: Path, genetic_code: int, refdir: Path, refseqver: str,
                    circular: bool, tool_options: dict, annotation_output_dir: Path) -> dict:
    """
    Resolves one sample's deterministic MITOS2 configuration. Pure
    function (no I/O beyond what's already in its arguments) -- the
    same inputs always resolve to the same configuration.

    Args:
        sample_id (str): The sample id.
        input_fasta (Path): The selected mitogenome FASTA to annotate.
        genetic_code (int): The requested genetic code.
        refdir (Path): The base taxonomy reference directory.
        refseqver (str): The reference data version.
        circular (bool): Whether to treat the sequence as circular.
        tool_options (dict): Boolean options (noplots, best, ncbicode, zip_output).
        annotation_output_dir (Path): Where the annotation stage will write this sample's output.

    Returns:
        dict: The resolved configuration, JSON-serializable.

    Raises:
        FileNotFoundError: If input_fasta doesn't exist.
        ValueError: If genetic_code isn't a positive integer, or refseqver is empty.
    """
    if not input_fasta.is_file():
        raise FileNotFoundError(f"Selected mitogenome FASTA does not exist: {input_fasta}.")

    if not isinstance(genetic_code, int) or genetic_code <= 0:
        raise ValueError(f"Invalid genetic code: {genetic_code!r}.")

    if not refseqver or not refseqver.strip():
        raise ValueError("No reference data version (refseqver) given.")

    reference_path = refdir / refseqver
    if not reference_path.is_dir():
        raise ValueError(f"No such reference data directory: {reference_path}.")

    return {
        "sample_id": sample_id,
        "input_fasta": str(input_fasta),
        "genetic_code": genetic_code,
        "refdir": str(refdir),
        "refseqver": refseqver,
        "circular": circular,
        "tool_options": dict(tool_options),
        "annotation_output_dir": str(annotation_output_dir),
    }

def main() -> int:
    """
    Resolves and records one sample's MITOS2 configuration.

    Returns:
        int: The return code for the process.
    """
    args = parse_args()

    logger = make_logger(
        name = "mitos2_config",
        log_file_path = args.log_file,
        global_log_file_path = args.global_log_file,
    )
    context = f"{{{args.sample_id}}}"

    try:
        config = resolve_config(
            sample_id = args.sample_id,
            input_fasta = args.input_fasta,
            genetic_code = args.genetic_code,
            refdir = args.refdir,
            refseqver = args.refseqver,
            circular = args.circular,
            tool_options = {
                "noplots": args.noplots, "best": args.best,
                "ncbicode": args.ncbicode, "zip_output": args.zip_output,
            },
            annotation_output_dir = args.annotation_output_dir,
        )
    except (ValueError, FileNotFoundError) as e:
        logger.error(f"{context} Could not resolve MITOS2 configuration: {e}")
        return 1

    args.output_dir.mkdir(parents = True, exist_ok = True)
    config_path = args.output_dir / f"{args.sample_id}.mitos2_config.json"
    with config_path.open("w", encoding = "utf-8") as handle:
        json.dump(config, handle, indent = 4)

    logger.info(f"{context} Resolved MITOS2 config: genetic_code={config['genetic_code']}, refseqver={config['refseqver']}. Written to {config_path}.")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
