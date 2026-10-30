"""
genbank_feature_table.py

Builds an NCBI GenBank 5-column feature table (.tbl) for one sample's
assembled mitochondrial genome from MITOS2's own real GFF annotation
-- the one artifact MITOS2 itself never produces natively (its own
output is GFF3/BED/FASTA/protein FASTA only, confirmed directly
against a real annotation run: no .tbl anywhere in its output
directory). This is the genuine missing piece between "MITOS2
annotated this genome" and "this genome is ready for table2asn."

Vertebrate mitochondrial genes are intron-less, so each gene's own
GFF coordinates ARE its CDS/tRNA/rRNA coordinates -- no splicing to
compute. Strand is encoded the way NCBI's 5-column format requires:
by coordinate order (a line's first column greater than its second
means minus strand), never a separate qualifier. An origin-crossing
feature (MITOS2 reports end > genome_length for these, never
start > end -- confirmed against circular_map.py's own real handling
of the same convention) is written as two location lines, the real
NCBI join convention for a circular-origin-spanning feature; only the
first line carries the feature key.
"""

# Imports
from __future__ import annotations
from logging import Logger
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from mitopipeline.utils.annotation.models import GeneFeature
from mitopipeline.utils.annotation.mitos_features import parse_mitos_gff

# Real NCBI gene symbol + product name for each of the 13 real
# vertebrate mitochondrial protein-coding genes, keyed by MITOS2's
# own lowercase gene_id (confirmed against a real GFF output). The
# same 13 names are what validation/gates.py's
# CORE_ANIMAL_MT_PROTEIN_GENES already checks for completeness --
# this table just attaches the real product description each needs
# for GenBank, which that set has no reason to carry.
PROTEIN_CODING_PRODUCTS: dict[str, tuple[str, str]] = {
    "atp6": ("ATP6", "ATP synthase F0 subunit 6"),
    "atp8": ("ATP8", "ATP synthase F0 subunit 8"),
    "cob": ("CYTB", "cytochrome b"),
    "cox1": ("COX1", "cytochrome c oxidase subunit I"),
    "cox2": ("COX2", "cytochrome c oxidase subunit II"),
    "cox3": ("COX3", "cytochrome c oxidase subunit III"),
    "nad1": ("ND1", "NADH dehydrogenase subunit 1"),
    "nad2": ("ND2", "NADH dehydrogenase subunit 2"),
    "nad3": ("ND3", "NADH dehydrogenase subunit 3"),
    "nad4": ("ND4", "NADH dehydrogenase subunit 4"),
    "nad4l": ("ND4L", "NADH dehydrogenase subunit 4L"),
    "nad5": ("ND5", "NADH dehydrogenase subunit 5"),
    "nad6": ("ND6", "NADH dehydrogenase subunit 6"),
}

# MITOS2's own rRNA gene_id convention, confirmed against a real GFF
# ("rrnS"/"rrnL"), matched case-insensitively.
RRNA_PRODUCTS: dict[str, tuple[str, str]] = {
    "rrns": ("s-rRNA", "12S ribosomal RNA"),
    "rrnl": ("l-rRNA", "16S ribosomal RNA"),
}

# MITOS2's own tRNA gene_id convention is "trn" plus a single-letter
# amino acid code, with a trailing "1"/"2" disambiguating the two
# real tRNAs for Leu/Ser (e.g. "trnG", "trnL1", "trnS2" -- confirmed
# against a real GFF output).
TRNA_AMINO_ACID_CODES: dict[str, str] = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys",
    "Q": "Gln", "E": "Glu", "G": "Gly", "H": "His", "I": "Ile",
    "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro",
    "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val",
}

# NCBI translation table 2 -- the real vertebrate mitochondrial
# genetic code. Every CDS here is translated under it, never the
# nuclear default (table 1): vertebrate mitochondria reassign AGA/
# AGG to stop and ATA to Met, among other real differences that a
# default-table translation would get wrong.
VERTEBRATE_MITOCHONDRIAL_GENETIC_CODE = 2


def trna_product(gene_id: str) -> tuple[str, str]:
    """
    Resolves a MITOS2 tRNA gene_id to its real (gene symbol, product)
    pair, e.g. "trnL1" -> ("trnL1", "tRNA-Leu"). MITOS2's own "tRNA"
    (not "ncRNA_gene") row embeds the real anticodon in the Name
    qualifier -- e.g. "trnG(tcc)" -- confirmed against a real GFF
    output; that's real, useful provenance but not a valid NCBI gene
    symbol, so it's stripped from the returned symbol here rather
    than written verbatim into a /gene qualifier.

    Args:
        gene_id (str): MITOS2's own tRNA gene_id (e.g. "trnG", "trnL1", "trnG(tcc)").

    Returns:
        tuple[str, str]: (gene symbol, product string).

    Raises:
        ValueError: If gene_id isn't a recognizable "trn"-prefixed amino acid code.
    """
    symbol = gene_id.split("(", 1)[0]
    if not symbol.lower().startswith("trn") or len(symbol) < 4:
        raise ValueError(f"Not a recognizable tRNA gene_id: {gene_id!r}.")
    code = symbol[3:].rstrip("12")
    letter = code[:1].upper()
    amino_acid = TRNA_AMINO_ACID_CODES.get(letter)
    if amino_acid is None:
        raise ValueError(f"Unrecognized tRNA amino-acid code in {gene_id!r}.")
    return symbol, f"tRNA-{amino_acid}"


def _location_segments(start: int, end: int, strand: str | None, genome_length: int) -> list[tuple[int, int]]:
    """
    Resolves one feature's real NCBI 5-column location column pairs,
    in the order they must be written. A feature confined to the
    linear sequence (end <= genome_length) is one segment; an
    origin-crossing feature (end > genome_length -- MITOS2's own
    convention, never start > end) is two, in real transcription
    order for whichever strand this feature is on. Minus strand is
    encoded by writing each segment's larger coordinate first, which
    is what tells table2asn/GenBank this is the complementary strand
    -- never a separate qualifier.

    Args:
        start (int): The feature's 1-based start, as MITOS2 reports it.
        end (int): The feature's 1-based end, as MITOS2 reports it (may exceed genome_length for an origin-crossing feature).
        strand (str | None): "+" or "-", or None (treated as "+" -- MITOS2 always reports a real strand for a kept feature).
        genome_length (int): The linear assembled sequence's own real length.

    Returns:
        list[tuple[int, int]]: [(col1, col2), ...] -- one tuple per location line, in the order they must appear.
    """
    minus = strand == "-"
    if end <= genome_length:
        return [(end, start)] if minus else [(start, end)]

    wrapped_end = end - genome_length
    if wrapped_end < 1 or wrapped_end > genome_length:
        raise ValueError(f"Invalid circular coordinates {start}-{end} for a {genome_length}bp genome.")

    # Plus-strand transcription order is start -> genome_length, then
    # 1 -> wrapped_end. Minus-strand transcription runs the opposite
    # direction through the same two spans, each one itself reversed.
    if minus:
        return [(wrapped_end, 1), (genome_length, start)]
    return [(start, genome_length), (1, wrapped_end)]


def _write_feature(lines: list[str], feature_key: str, start: int, end: int, strand: str | None,
                    genome_length: int, qualifiers: list[tuple[str, str]]) -> None:
    """Appends one feature's real location line(s) plus its qualifier lines to lines, in real 5-column .tbl format."""
    segments = _location_segments(start, end, strand, genome_length)
    for index, (col1, col2) in enumerate(segments):
        key = feature_key if index == 0 else ""
        lines.append(f"{col1}\t{col2}\t{key}")
    for qualifier_key, qualifier_value in qualifiers:
        lines.append(f"\t\t\t{qualifier_key}\t{qualifier_value}")


def build_feature_table(sample_id: str, seqid: str, genes: list[GeneFeature], trnas: list[GeneFeature],
                         rrnas: list[GeneFeature], genome_length: int, organism: str,
                         logger: Logger | None = None) -> tuple[str, list[str]]:
    """
    Builds one sample's real NCBI 5-column feature table text from
    its already-parsed MITOS2 gene/tRNA/rRNA features (the shared
    parser output). Protein-coding genes get a "gene" line plus a
    matching "CDS" line (transl_table=2, codon_start=1 -- vertebrate
    mitochondrial genes are intron-less and MITOS2's own coordinates
    already mark the real start/stop, so no frame inference is
    needed). tRNA/rRNA genes get a single "gene" plus matching
    "tRNA"/"rRNA" line pair. Every gene not in PROTEIN_CODING_PRODUCTS/
    RRNA_PRODUCTS or TRNA_AMINO_ACID_CODES is skipped with a logged
    warning rather than guessed at.

    Args:
        sample_id (str): The sample id (only used for logging/warnings).
        seqid (str): The real sequence id the ">Feature" header names.
        genes (list[GeneFeature]): MITOS2's own "gene" (protein-coding) features.
        trnas (list[GeneFeature]): MITOS2's own "tRNA" features.
        rrnas (list[GeneFeature]): MITOS2's own "rRNA" features.
        genome_length (int): The assembled sequence's own real length.
        organism (str): The real organism name for the source feature.
        logger (Logger | None): Logger to use. Defaults to None.

    Returns:
        tuple[str, list[str]]: (the real .tbl text, a list of real warnings for any feature skipped or ambiguous).
    """
    warnings: list[str] = []
    lines = [f">Feature {seqid}"]
    _write_feature(lines, "source", 1, genome_length, "+", genome_length,
                    [("organism", organism), ("mol_type", "genomic DNA"), ("organelle", "mitochondrion")])

    tagged = [(g, "gene") for g in genes] + [(g, "rRNA") for g in rrnas] + [(g, "tRNA") for g in trnas]
    tagged.sort(key = lambda pair: pair[0].start)

    for gene, category in tagged:
        key = gene.gene_id.lower()
        if category == "gene":
            product_entry = PROTEIN_CODING_PRODUCTS.get(key)
            if product_entry is None:
                warnings.append(f"{sample_id}: unrecognized protein-coding gene {gene.gene_id!r}, skipped.")
                if logger is not None: logger.warning(warnings[-1])
                continue
            symbol, product = product_entry
            _write_feature(lines, "gene", gene.start, gene.end, gene.strand, genome_length, [("gene", symbol)])
            _write_feature(lines, "CDS", gene.start, gene.end, gene.strand, genome_length,
                            [("gene", symbol), ("product", product), ("codon_start", "1"),
                             ("transl_table", str(VERTEBRATE_MITOCHONDRIAL_GENETIC_CODE))])
        elif category == "rRNA":
            product_entry = RRNA_PRODUCTS.get(key)
            if product_entry is None:
                warnings.append(f"{sample_id}: unrecognized rRNA gene {gene.gene_id!r}, skipped.")
                if logger is not None: logger.warning(warnings[-1])
                continue
            symbol, product = product_entry
            _write_feature(lines, "gene", gene.start, gene.end, gene.strand, genome_length, [("gene", symbol)])
            _write_feature(lines, "rRNA", gene.start, gene.end, gene.strand, genome_length,
                            [("gene", symbol), ("product", product)])
        else:
            try:
                symbol, product = trna_product(gene.gene_id)
            except ValueError:
                warnings.append(f"{sample_id}: unrecognized tRNA gene {gene.gene_id!r}, skipped.")
                if logger is not None: logger.warning(warnings[-1])
                continue
            _write_feature(lines, "gene", gene.start, gene.end, gene.strand, genome_length, [("gene", symbol)])
            _write_feature(lines, "tRNA", gene.start, gene.end, gene.strand, genome_length,
                            [("gene", symbol), ("product", product)])

    return "\n".join(lines) + "\n", warnings


def _extract_cds_nucleotides(sequence: str, start: int, end: int, strand: str | None, genome_length: int) -> str:
    """Extracts one CDS's real nucleotide sequence (5'->3', coding strand), handling an origin-crossing span."""
    minus = strand == "-"
    if end <= genome_length:
        bases = sequence[start - 1:end]
    else:
        wrapped_end = end - genome_length
        bases = sequence[start - 1:genome_length] + sequence[0:wrapped_end]
    if minus:
        return str(Seq(bases).reverse_complement())
    return bases


def check_cds_translations(genes: list[GeneFeature], genome_fasta: str | Path, genome_length: int,
                            logger: Logger | None = None) -> list[str]:
    """
    Translates every real protein-coding CDS under the vertebrate
    mitochondrial genetic code (table 2) and flags a genuine problem
    -- a premature internal stop codon -- as a warning string, never
    a fabricated pass. Does NOT flag a non-stop terminal codon as a
    problem: an incomplete stop codon completed by post-transcriptional
    polyadenylation is real, well-documented vertebrate mitochondrial
    biology (e.g. a terminal "T"/"TA" that becomes TAA), not an
    annotation defect -- table2asn's own real validation is the
    authoritative check for anything subtler than this.

    Args:
        genes (list[GeneFeature]): The real protein-coding gene features to check.
        genome_fasta (str | Path): The assembled genome FASTA MITOS2 annotated.
        genome_length (int): The assembled sequence's own real length.
        logger (Logger | None): Logger to use. Defaults to None.

    Returns:
        list[str]: One warning per CDS with a real premature internal stop codon. Empty when every CDS translates cleanly.
    """
    records = list(SeqIO.parse(Path(genome_fasta), "fasta"))
    if not records:
        raise ValueError(f"Genome FASTA contains no records: {genome_fasta}.")
    sequence = str(records[0].seq).upper()

    warnings: list[str] = []
    for gene in genes:
        key = gene.gene_id.lower()
        if key not in PROTEIN_CODING_PRODUCTS:
            continue
        nucleotides = _extract_cds_nucleotides(sequence, gene.start, gene.end, gene.strand, genome_length)
        usable_length = len(nucleotides) - (len(nucleotides) % 3)
        protein = str(Seq(nucleotides[:usable_length]).translate(
            table = VERTEBRATE_MITOCHONDRIAL_GENETIC_CODE, to_stop = False))
        if "*" in protein[:-1]:
            message = (f"{gene.gene_id}: premature internal stop codon under the vertebrate mitochondrial "
                       f"genetic code (table {VERTEBRATE_MITOCHONDRIAL_GENETIC_CODE}) -- real annotation "
                       f"problem, not a formatting issue.")
            warnings.append(message)
            if logger is not None: logger.warning(message)
    return warnings


def build_sample_feature_table(sample_id: str, gff_path: str | Path, genome_fasta: str | Path, organism: str,
                                output_tbl_path: str | Path, logger: Logger | None = None) -> list[str]:
    """
    Builds and writes one sample's complete, real GenBank feature
    table from its real MITOS2 GFF and assembled genome FASTA, and
    checks every real CDS's translation under the vertebrate
    mitochondrial genetic code.

    Args:
        sample_id (str): The sample id.
        gff_path (str | Path): MITOS2's own real result GFF.
        genome_fasta (str | Path): The assembled genome FASTA MITOS2 annotated.
        organism (str): The real organism name for the source feature.
        output_tbl_path (str | Path): Where to write the real .tbl file.
        logger (Logger | None): Logger to use. Defaults to None.

    Returns:
        list[str]: Every real warning found (an unrecognized feature, or a real premature internal stop codon).
    """
    records = list(SeqIO.parse(Path(genome_fasta), "fasta"))
    if not records:
        raise ValueError(f"Genome FASTA contains no records: {genome_fasta}.")
    genome_length = len(str(records[0].seq))
    seqid = records[0].id

    genes = parse_mitos_gff(gff_path, logger = logger, feature_type = "gene")
    trnas = parse_mitos_gff(gff_path, logger = logger, feature_type = "tRNA")
    rrnas = parse_mitos_gff(gff_path, logger = logger, feature_type = "rRNA")

    table_text, warnings = build_feature_table(sample_id, seqid, genes, trnas, rrnas, genome_length, organism, logger)
    warnings += check_cds_translations(genes, genome_fasta, genome_length, logger)

    output_tbl_path = Path(output_tbl_path)
    output_tbl_path.parent.mkdir(parents = True, exist_ok = True)
    output_tbl_path.write_text(table_text, encoding = "utf-8")
    if logger is not None:
        logger.info(f"{{{sample_id}}} Wrote GenBank feature table to {output_tbl_path} "
                    f"({len(genes)} CDS, {len(trnas)} tRNA, {len(rrnas)} rRNA, {len(warnings)} warning(s)).")
    return warnings
