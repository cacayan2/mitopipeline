"""
manifest_parser.py

Parses and validates a pipeline sample manifest.
"""

# Imports
from __future__ import annotations
import logging
from pathlib import Path
import pandas as pd
from mitopipeline.models.sample import Sample

REQUIRED_COLUMNS = ["sample_id", "r1", "r2"]

def parse_sample_manifest(
        manifest_path: Path,
        logger: logging.Logger | None = None,
        input_directory: Path | None = None
) -> list[Sample]:
    """
    Parses a TSV manifest while retaining every nonblank column.

    Args:
        manifest_path (Path): The path to the manifest file.
        logger (logging.Logger | None): Logger for logging messages.
        input_directory (Path | None): The input directory where the fastq files live.

    Returns:
        list[Sample]: A list of Sample objects parsed from the manifest.
    """
    # First, we normalize input_directory.
    if input_directory is not None:
        input_directory = Path(input_directory).expanduser().resolve()

    # First, we normalize the path to the manifest and check if it exists.
    manifest_path = Path(manifest_path).resolve()
    if not manifest_path.exists():
        if logger is not None: logger.error(f"Sample manifest file does not exist: {manifest_path}")
        raise FileNotFoundError(f"Sample manifest file does not exist: {manifest_path}")

    # Then we verify if the manifest is nonempty.
    if manifest_path.stat().st_size == 0:
        if logger is not None: logger.error(f"Sample manifest file is empty: {manifest_path}")
        raise ValueError(f"Sample manifest file is empty: {manifest_path}")

    # Now we read the manifest into a pandas DataFrame.
    if logger is not None: logger.info(f"Reading sample manifest from: {manifest_path}.")
    table = pd.read_csv(manifest_path, sep = "\t", dtype = str, keep_default_na = False)

    # Normalizing the column names by stripping whitespace from each column name.
    table.columns = [str(column).strip() for column in table.columns]

    # We verify that all required columns are present in the manifest.
    if table.empty:
        if logger is not None: logger.error(f"Sample manifest file contains no columns: {manifest_path}")
        raise ValueError(f"Sample manifest file contains no columns: {manifest_path}")

    # Here we verify that all required columns are present in the manifest. We subtract the set of required columns
    # from the set of columns in the manifest - if any remain, those are missing.
    missing_columns = set(REQUIRED_COLUMNS) - set(table.columns)
    if missing_columns:
        if logger is not None: logger.error(f"Sample manifest file is missing required columns: {missing_columns}")
        raise ValueError(f"Sample manifest file is missing required columns: {missing_columns}")

    # Now we start construction of the sample object from each row of the manifest. First we obtain the sample id's
    # and verify that there aren't any duplicates.
    sample_ids = table["sample_id"].astype(str).str.strip()
    duplicate_ids = sorted(sample_ids[sample_ids.duplicated(keep = False)].unique())
    if duplicate_ids:
        if logger is not None: logger.error(f"Sample manifest file contains duplicate sample_id's: {duplicate_ids}")
        raise ValueError(f"Sample manifest file contains duplicate sample_id's: {duplicate_ids}")

    # Now we create a list of sample objects from the manifest and define the parent directory of the manifest.
    samples: list[Sample] = []
    manifest_dir = manifest_path.parent

    # We iterate through each row and add R1 and R2, validating along the way.
    for index, row in table.iterrows():
        line_number = index + 2  # +2 because pandas is zero-indexed and we have a header row
        metadata = _row_metadata(row)

        # We verify that sample_id is present and nonblank.
        sample_id = metadata.get("sample_id", "")
        if not sample_id:
            if logger is not None: logger.error(f"Missing sample_id at line {line_number} in manifest.")
            raise ValueError(f"Missing sample_id at line {line_number} in manifest.")

        # We verify that r1 and r2 are present and nonblank.
        r1_value = metadata.get("r1", "")
        r2_value = metadata.get("r2", "")
        if not r1_value:
            if logger is not None: logger.error(f"Missing r1 at line {line_number} in manifest.")
            raise ValueError(f"Missing r1 at line {line_number} in manifest.")
        if not r2_value:
            if logger is not None: logger.error(f"Missing r2 at line {line_number} in manifest.")
            raise ValueError(f"Missing r2 at line {line_number} in manifest.")

        # Now we resolve the manifest paths to absolute paths.
        r1_path = _resolve_manifest_path(value = r1_value, manifest_dir = manifest_dir, input_directory = input_directory)
        r2_path = _resolve_manifest_path(value = r2_value, manifest_dir = manifest_dir, input_directory = input_directory)

        # Verifying if the paths lead to a file.
        if not r1_path.is_file():
            if logger is not None: logger.error(f"Input file not found: {r1_path}.")
            raise FileNotFoundError(f"Input file not found: {r1_path}.")
        if not r2_path.is_file():
            if logger is not None: logger.error(f"Input file not found: {r2_path}.")
            raise FileNotFoundError(f"Input file not found: {r2_path}.")

        # Now adding the r1 and r2 to the metadata data table.
        metadata["r1"] = str(r1_path)
        metadata["r2"] = str(r2_path)

        # Now we construct a Sample object and add to the samples list.
        sample = Sample(sample_id = sample_id,
                        r1 = r1_path,
                        r2 = r2_path,
                        genus = _optional_value(metadata, "genus"),
                        species = _optional_value(metadata, "species"),
                        source = _optional_value(metadata, "source"),
                        metadata = metadata)
        samples.append(sample)

        # We log that we have validated a sample.
        if logger is not None:
            logger.info("[%s] Validated paired FASTQ input paths and retained %d metdata fields.", sample_id, len(sample.metadata)) 
            logger.debug("[%s] Manifest metadata: %s", sample_id, dict(sample.metadata))

    # Log that we logged all samples from the manifest.
    if logger is not None:
        logger.info("Validated %d samples from %s.", len(samples), manifest_path)

    return samples        

def _row_metadata(row: pd.Series) -> dict[str, str]:
    """
    Extracts nonblank metadata from a row of the manifest, excluding required columns.

    Args:
        row (pd.Series): A row from the manifest DataFrame.

    Returns:
        dict[str, str]: A dictionary of metadata fields and their values.
    """
    # We start with an empty dictionary to store the metadata.
    metadata: dict[str, str] = {}

    # This for loop does the following:
    # 1. Iterates through each column and value in the row.
    # 2. Strips whitespace from the column name and value.
    # 3. If either the column name or value is empty, it skips to the next iteration.
    # 4. Otherwise, it adds the column name and value to the metadata dictionary.
    for column, value in row.items():
        key = str(column).strip()
        text = str(value).strip()
        if not key or not text: continue
        metadata[key] = text

    # Finally, we return the metadata dictionary.
    return metadata

def _resolve_manifest_path(value: str, manifest_dir: Path, input_directory: Path | None = None) -> Path:
    """
    Resolves a path from the manifest, handling both absolute and relative paths.

    Args:
        value (str): The path value from the manifest.
        manifest_dir (Path): The directory of the manifest file.
        input_directory (Path | None = None): The input directory.

    Returns:
        Path: The resolved path.
    """
    # We first strip whitespace from the value.
    path = Path(str(value).strip()).expanduser()

    # If the path is absolute, we resolve it.
    if path.is_absolute(): return path.resolve()

    # Logic for dealing with a relative path.
    ## 1. We resolve the relative path and return it if it is a dictionary.
    ## 2. If not, we resolve the relative path to the input_directory and return it. 
    manifest_relative = (manifest_dir / path).resolve()
    if manifest_relative.is_file(): return manifest_relative
    if input_directory is not None:
        input_relative = (input_directory / path).resolve()
        if input_relative.is_file(): return input_relative
    return manifest_relative

def _optional_value(metadata: dict[str, str], column_name: str) -> str | None:
    """
    Retrieves an optional value from the metadata dictionary.

    Args:
        metadata (dict[str, str]): The metadata dictionary.
        column_name (str): The name of the column to retrieve.

    Returns:
        str | None: The value of the column if present and nonblank, otherwise None.
    """
    # We attempt to get the value from the metadata dictionary.
    value = metadata.get(column_name)

    # Returns the value if it is not none.
    return value if value else None