"""
reference_database.py

Builds and validates a local nucleotide BLAST reference database of
complete mitochondrial genomes, scoped to a taxonomic group broad
enough to include real alternative/confound species -- so a correct
match is genuine evidence, not a foregone conclusion -- and pulled
exhaustively within that scope rather than an arbitrary small sample
(the actual root cause of a historical bug: BLAST hits coming back as
distantly related organisms, because the old default query had no
taxonomic targeting at all and capped at 250 records across every
RefSeq mitochondrial genome in all of life).
"""

# Imports
import json
import shutil
import subprocess
import tempfile
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from datetime import datetime
from logging import Logger
from pathlib import Path

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# The taxonomic scope broad enough to include real confound species
# (e.g. marine fish, when the target is freshwater) rather than being
# narrowed to the sample's own lineage -- narrowing to the answer
# being tested would make a "match" prove nothing.
# A single named constant, never hard-coded past this point.
DEFAULT_TAXON_SCOPE = "Actinopterygii"

# A high ceiling, not a real limit -- NCBI's own esearch supports this
# easily, well above the real record count for any one well-studied
# clade (~3900 for Actinopterygii, measured 2026-10-19). The old
# v0.1.0 default of 250 was the actual bug: an arbitrary small
# subsample of a much larger, untargeted pool (17716 records with no
# taxon filter at all).
MAX_RECORDS_CEILING = 10_000

def build_ncbi_query(taxon_scope: str = DEFAULT_TAXON_SCOPE) -> str:
    """
    Builds the NCBI Entrez query for complete RefSeq mitochondrial
    genomes within one taxonomic scope.

    Args:
        taxon_scope (str): The NCBI organism-level scope to search within. Defaults to DEFAULT_TAXON_SCOPE.

    Returns:
        str: The Entrez query string.
    """
    return f'{taxon_scope}[Organism] AND mitochondrion[filter] AND refseq[filter] AND "complete genome"[Title]'

def run_command(command: list[str], logger: Logger, cwd: Path | None = None) -> subprocess.CompletedProcess:
    """
    Runs a command and captures its output.

    Args:
        command (list[str]): The command to run.
        logger (Logger): The logger to use.
        cwd (Path | None): The working directory to run the command in.

    Returns:
        subprocess.CompletedProcess: The completed process.
    """
    logger.debug(f"Command: {command}")
    result = subprocess.run(command, cwd = cwd, capture_output = True, text = True)
    logger.debug(f"Return code: {result.returncode}")
    if result.stderr:
        logger.debug(f"stderr: {result.stderr}")
    return result

def fetch_url_text(url: str, logger: Logger, retries: int = 3, delay_seconds: float = 1.0) -> str:
    """
    Fetches text from a URL with simple retry handling.

    Args:
        url (str): The URL to fetch.
        logger (Logger): The logger to use.
        retries (int): The number of attempts before giving up.
        delay_seconds (float): The base delay between retries, scaled by attempt number.

    Returns:
        str: The decoded response text.

    Raises:
        Exception: Whatever the final attempt's failure was, if every retry is exhausted.
    """
    for attempt in range(1, retries + 1):
        try:
            with urllib.request.urlopen(url, timeout = 120) as response:
                return response.read().decode("utf-8")
        except Exception as error:
            logger.warning(f"URL fetch failed on attempt {attempt}/{retries}: {error}")
            if attempt == retries:
                raise
            time.sleep(delay_seconds * attempt)
    raise RuntimeError("Unreachable URL retry state.")

def ncbi_esearch(query: str, email: str, api_key: str | None, max_records: int, logger: Logger) -> list[str]:
    """
    Searches NCBI nucleotide and returns matching record IDs.

    Args:
        query (str): The Entrez query to run.
        email (str): The contact email NCBI requires for Entrez requests.
        api_key (str | None): An optional NCBI API key, for a higher rate limit.
        max_records (int): The maximum number of IDs to retrieve -- a safety ceiling, not meant to be the limiting factor.
        logger (Logger): The logger to use.

    Returns:
        list[str]: The matching record IDs.
    """
    parameters = {"db": "nuccore", "term": query, "retmode": "xml", "retmax": str(max_records),
                  "email": email, "tool": "mitopipeline"}
    if api_key:
        parameters["api_key"] = api_key
    url = f"{EUTILS_BASE}/esearch.fcgi?" + urllib.parse.urlencode(parameters)
    text = fetch_url_text(url, logger = logger)
    root = ET.fromstring(text)
    ids = [element.text for element in root.findall(".//IdList/Id") if element.text]
    count_text = root.findtext(".//Count")
    logger.info(f"NCBI Entrez query matched {count_text or 'unknown'} records; retrieving {len(ids)}.")
    if count_text and int(count_text) > max_records:
        logger.warning(f"Query matched more records ({count_text}) than max_records ({max_records}) allows; "
                        f"the reference set will be incomplete. Raise max_records to cover the full scope.")
    return ids

def ncbi_efetch_fasta(ids: list[str], email: str, api_key: str | None, batch_size: int, logger: Logger) -> str:
    """
    Fetches nucleotide FASTA records from NCBI, in batches.

    Args:
        ids (list[str]): The record IDs to fetch.
        email (str): The contact email NCBI requires for Entrez requests.
        api_key (str | None): An optional NCBI API key, for a higher rate limit.
        batch_size (int): How many records to fetch per request.
        logger (Logger): The logger to use.

    Returns:
        str: The concatenated FASTA text for every fetched record.
    """
    fasta_chunks: list[str] = []
    for start in range(0, len(ids), batch_size):
        batch = ids[start:start + batch_size]
        parameters = {"db": "nuccore", "id": ",".join(batch), "rettype": "fasta", "retmode": "text",
                      "email": email, "tool": "mitopipeline"}
        if api_key:
            parameters["api_key"] = api_key
        url = f"{EUTILS_BASE}/efetch.fcgi?" + urllib.parse.urlencode(parameters)
        logger.info(f"Fetching NCBI FASTA records {start + 1}-{start + len(batch)} of {len(ids)}.")
        fasta_chunks.append(fetch_url_text(url, logger = logger))
        # Being polite to NCBI when no API key is configured.
        if not api_key:
            time.sleep(0.4)
    return "\n".join(chunk.strip() for chunk in fasta_chunks if chunk.strip()) + "\n"

def count_fasta_records(fasta_path: Path) -> int:
    """
    Counts FASTA records in a file.

    Args:
        fasta_path (Path): The FASTA file to count records in.

    Returns:
        int: The number of FASTA records (lines starting with ">").
    """
    count = 0
    with fasta_path.open("r", encoding = "utf-8", errors = "replace") as handle:
        for line in handle:
            if line.startswith(">"):
                count += 1
    return count

def download_reference_fasta(reference_fasta: Path, taxon_scope: str, email: str, api_key: str | None,
                              max_records: int, batch_size: int, overwrite: bool, logger: Logger) -> int:
    """
    Downloads every complete RefSeq mitochondrial genome within
    taxon_scope from NCBI Entrez and writes them to reference_fasta.

    Args:
        reference_fasta (Path): The path to write the downloaded FASTA to.
        taxon_scope (str): The NCBI organism-level scope to search within.
        email (str): The contact email NCBI requires for Entrez requests.
        api_key (str | None): An optional NCBI API key.
        max_records (int): The ceiling on records to retrieve -- a safety bound; a warning is logged if the real match count exceeds it.
        batch_size (int): How many records to fetch per efetch request.
        overwrite (bool): Whether to re-download even if reference_fasta already has records.
        logger (Logger): The logger to use.

    Returns:
        int: 0 on success, 1 on failure.
    """
    if reference_fasta.exists() and not overwrite:
        record_count = count_fasta_records(reference_fasta)
        if record_count > 0:
            logger.info(f"Existing reference FASTA found with {record_count} records: {reference_fasta}")
            return 0

    if not email:
        logger.error("NCBI Entrez download requires an email address.")
        return 1

    reference_fasta.parent.mkdir(parents = True, exist_ok = True)
    query = build_ncbi_query(taxon_scope)

    ids = ncbi_esearch(query = query, email = email, api_key = api_key, max_records = max_records, logger = logger)
    if not ids:
        logger.error(f"NCBI query returned no records: {query}")
        return 1

    fasta_text = ncbi_efetch_fasta(ids = ids, email = email, api_key = api_key, batch_size = batch_size, logger = logger)
    if not fasta_text.strip():
        logger.error("NCBI efetch returned an empty FASTA response.")
        return 1

    temporary_fasta = reference_fasta.with_suffix(reference_fasta.suffix + ".tmp")
    temporary_fasta.write_text(fasta_text, encoding = "utf-8")

    record_count = count_fasta_records(temporary_fasta)
    if record_count == 0:
        logger.error("Downloaded reference FASTA contains zero records.")
        temporary_fasta.unlink(missing_ok = True)
        return 1

    shutil.move(str(temporary_fasta), str(reference_fasta))
    logger.info(f"Wrote {record_count} {taxon_scope} mitochondrial references to {reference_fasta}.")
    return 0

def validate_database(database_prefix: Path, logger: Logger) -> bool:
    """
    Validates a local BLAST database using blastdbcmd.

    Args:
        database_prefix (Path): The database prefix to validate.
        logger (Logger): The logger to use.

    Returns:
        bool: True if the database is valid.
    """
    result = run_command(["blastdbcmd", "-db", str(database_prefix), "-dbtype", "nucl", "-info"], logger = logger)
    return result.returncode == 0

def build_blast_database(reference_fasta: Path, database_prefix: Path, database_title: str,
                          parse_seqids: bool, logger: Logger) -> int:
    """
    Builds a local nucleotide BLAST database from a reference FASTA.

    Args:
        reference_fasta (Path): The FASTA file to build the database from.
        database_prefix (Path): The path prefix for the built database files.
        database_title (str): A human-readable title recorded in the database.
        parse_seqids (bool): Whether to preserve original sequence IDs.
        logger (Logger): The logger to use.

    Returns:
        int: 0 on success, the underlying makeblastdb return code otherwise.
    """
    database_prefix.parent.mkdir(parents = True, exist_ok = True)
    with tempfile.TemporaryDirectory(dir = database_prefix.parent) as temporary_directory:
        temporary_prefix = Path(temporary_directory) / database_prefix.name
        command = ["makeblastdb", "-in", str(reference_fasta), "-dbtype", "nucl",
                   "-out", str(temporary_prefix), "-title", database_title]
        if parse_seqids:
            command.append("-parse_seqids")

        result = run_command(command, logger = logger)
        if result.returncode != 0:
            logger.error(f"makeblastdb failed: {result.stderr.strip()}")
            return result.returncode

        if not validate_database(temporary_prefix, logger):
            logger.error("Newly built BLAST database failed validation.")
            return 1

        for source_path in sorted(Path(temporary_directory).glob(f"{temporary_prefix.name}*")):
            suffix = source_path.name[len(temporary_prefix.name):]
            shutil.move(str(source_path), str(Path(f"{database_prefix}{suffix}")))

    return 0

def write_metadata(metadata_file: Path, taxon_scope: str, database_prefix: Path,
                    reference_fasta: Path, query: str) -> None:
    """
    Writes database setup metadata for provenance.

    Args:
        metadata_file (Path): The path to write metadata to.
        taxon_scope (str): The taxonomic scope the database was built from.
        database_prefix (Path): The built database's path prefix.
        reference_fasta (Path): The reference FASTA the database was built from.
        query (str): The exact NCBI Entrez query used.

    Returns:
        None
    """
    metadata = {
        "taxon_scope": taxon_scope,
        "ncbi_query": query,
        "database_prefix": str(database_prefix.resolve()),
        "reference_fasta": str(reference_fasta.resolve()),
        "reference_fasta_record_count": count_fasta_records(reference_fasta),
        "created_at": datetime.now().isoformat(),
    }
    metadata_file.parent.mkdir(parents = True, exist_ok = True)
    metadata_file.write_text(json.dumps(metadata, indent = 2), encoding = "utf-8")

def setup_reference_database(reference_fasta: Path, database_prefix: Path, database_title: str,
                              metadata_file: Path, email: str, logger: Logger,
                              taxon_scope: str = DEFAULT_TAXON_SCOPE, api_key: str | None = None,
                              max_records: int = MAX_RECORDS_CEILING, batch_size: int = 200,
                              parse_seqids: bool = True, overwrite: bool = False) -> int:
    """
    Downloads (if needed) and builds a local BLAST reference database
    of complete mitochondrial genomes for one taxonomic scope.

    Args:
        reference_fasta (Path): The path to write/read the reference FASTA.
        database_prefix (Path): The path prefix for the built database.
        database_title (str): A human-readable title recorded in the database.
        metadata_file (Path): The path to write setup metadata to.
        email (str): The contact email NCBI requires for Entrez requests.
        logger (Logger): The logger to use.
        taxon_scope (str): The NCBI organism-level scope to search within. Defaults to DEFAULT_TAXON_SCOPE.
        api_key (str | None): An optional NCBI API key.
        max_records (int): The ceiling on records to retrieve. Defaults to MAX_RECORDS_CEILING.
        batch_size (int): How many records to fetch per efetch request.
        parse_seqids (bool): Whether to preserve original sequence IDs.
        overwrite (bool): Whether to rebuild even if a valid database already exists.

    Returns:
        int: 0 on success, non-zero otherwise.
    """
    if not overwrite and validate_database(database_prefix, logger):
        logger.info(f"Existing BLAST database is already valid: {database_prefix}")
        return 0

    return_code = download_reference_fasta(
        reference_fasta = reference_fasta, taxon_scope = taxon_scope, email = email, api_key = api_key,
        max_records = max_records, batch_size = batch_size, overwrite = overwrite, logger = logger,
    )
    if return_code != 0:
        return return_code

    return_code = build_blast_database(
        reference_fasta = reference_fasta, database_prefix = database_prefix,
        database_title = database_title, parse_seqids = parse_seqids, logger = logger,
    )
    if return_code != 0:
        return return_code

    write_metadata(metadata_file = metadata_file, taxon_scope = taxon_scope, database_prefix = database_prefix,
                    reference_fasta = reference_fasta, query = build_ncbi_query(taxon_scope))
    return 0
