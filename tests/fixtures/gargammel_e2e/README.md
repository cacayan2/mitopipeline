# Gargammel-simulated end-to-end fixture set

Synthetic, realistic paired-end Illumina reads for a full pipeline
end-to-end smoke test (issue #403), requested by Dr. Wheeler: a small,
downloadable test dataset that can run raw QC -> trimming -> post-QC ->
GetOrganelle assembly and produce real, expected outputs without
needing an actual museum specimen.

Three samples, one per assembly outcome, so the smoke test (and any
later issue that validates assembly-status handling, e.g. #405/#422's
rescue-trigger logic) has a real fixture for each branch. All three
have now been run end to end through the pipeline's own exec scripts
(not just raw GetOrganelle) and confirmed to produce exactly the
outcome below:

| sample_id             | source                         | target coverage | confirmed outcome |
|------------------------|---------------------------------|------------------|--------------------|
| `gargammel_circular`   | full NC_002333.2 (16596bp)      | ~150x bulk + 60 origin-spanning pairs | `circular`, exact 16596bp assembly |
| `gargammel_fragmented` | full NC_002333.2 (16596bp)      | ~20x             | `fragmented` (1 scaffold, no complete candidate) |
| `gargammel_incomplete` | first 6000bp of NC_002333.2 only | ~2x              | `missing` (GetOrganelle: "Too few seed reads found!") |

Total size: 2.0MB (well under GitHub's recommended per-file and
repository size guidance).

## Source sequence

[NC_002333.2](https://www.ncbi.nlm.nih.gov/nuccore/NC_002333.2) —
*Danio rerio* (zebrafish) mitochondrion, complete genome, 16,596bp.
Fetched directly from NCBI via the EUtils API (`efetch`,
`db=nuccore`, `rettype=fasta`). Kept here for reference/reproducibility
as `source_reference_NC_002333.2.fasta`. Zebrafish was chosen as a
real, well-characterized, freely available reference — not a stand-in
for this project's actual target taxa (Actinopterygii more broadly);
the point of this fixture is exercising the PIPELINE's wiring, not
producing a scientifically meaningful zebrafish result.

## How these were generated

[Gargammel](https://github.com/grenaud/gargammel) v1.1.4 (bioconda),
which simulates realistic Illumina paired-end reads from a target
("endogenous") FASTA by fragmenting it and running it through ART's
read-simulation profiles. No contamination/bacterial component was
simulated (`--comp 0,0,1`, the default -- 100% endogenous), since the
point here is assembly-outcome variation, not contamination handling.

```bash
# Directory layout gargammel expects per sample:
#   <workdir>/endo/<target>.fasta   (bact/ and cont/ present but empty)

# fragmented -- full genome, low coverage
gargammel -c 20 -l 200 -rl 100 --comp 0,0,1 -o fragmented/sim fragmented/

# incomplete -- truncated genome (first 6000bp only), low coverage
gargammel -c 2 -l 200 -rl 100 --comp 0,0,1 -o incomplete/sim incomplete/
```

`-l 200` explicitly sets gargammel's fragment length above the
requested 100bp read length (`-rl 100`). Without it, gargammel
defaults to ~35bp fragments, so each "100bp read" is mostly Illumina
adapter sequence padding -- fastp's own adapter trimming then strips
nearly the whole read as adapter and discards it as too-short before
GetOrganelle ever sees it. All three samples here had to be generated
(or regenerated) with this flag set explicitly.

### The circular sample needed a second step

A linear reference FASTA can never produce a read pair spanning a
circular genome's origin (position 16596 wrapping to position 1) --
gargammel fragments the linear sequence as given, so no fragment ever
crosses that junction. Without such reads, GetOrganelle cannot resolve
a closed, circular path no matter how high the coverage, and reports a
fragmented/incomplete graph instead. (A tried workaround --
physically circularizing the input FASTA before handing it to
gargammel -- made gargammel's `fragSim` step hang indefinitely and was
abandoned.)

The fix actually used: run gargammel's normal bulk simulation on the
plain linear reference at high coverage, then separately generate a
small number of synthetic "origin-spanning" read pairs by sampling
directly from the doubled sequence (`(genome + genome)[start:start +
insert_size]` for insert sizes/positions that straddle the wrap point),
and concatenate the two read sets:

```bash
gargammel -c 150 -l 200 -rl 100 --comp 0,0,1 -o circular/sim circular/
# + 60 manually-generated origin-spanning read pairs, concatenated in
```

This combination (12,507 read pairs total) is what's actually in
`circular/gargammel_circular_R1/R2.fastq.gz` -- confirmed to produce
`Result status of animal_mt: circular genome` with an assembled length
of exactly 16,596bp, matching the source genome.

Illumina HiSeq 2500 profile (gargammel's default), 100bp paired-end
reads -- comfortably above this pipeline's default fastp
`length_required` (50bp), so these aren't filtered out by trimming.

Read pair counts, directly measured: circular 12,507 pairs, fragmented
1,659 pairs, incomplete 60 pairs.

## Using this fixture set

`manifest.tsv` in this directory is a ready-to-use pipeline manifest
(`sample_id`, `r1`, `r2`, plus descriptive metadata columns the parser
preserves but doesn't require) with paths relative to this directory,
so pointing the pipeline's `manifest` config key at it needs no other
setup:

```yaml
manifest: tests/fixtures/gargammel_e2e/manifest.tsv
```

## Verified end to end

All three samples have been run through the pipeline's real exec
scripts in sequence (qc_raw -> trimming -> qc_trimmed ->
getorganelle_config -> getorganelle_assembly, not just raw
`get_organelle_from_reads.py`), and each produced the outcome in the
table above. This run also caught and fixed two real wiring bugs in
`GetOrganelleRunner` that had never been exercised by a real
invocation before: `build_command()` was pre-creating the output
directory GetOrganelle itself requires to not exist yet, and the
assembly classification was being (re-)computed after
`postprocess_outputs()` had already deleted the native files it looks
for. Both are fixed in `src/mitopipeline/api/getorganelle.py`.

`gargammel_incomplete`'s `missing` outcome is GetOrganelle's own
genuine failure ("Too few seed reads found!"), not a pipeline crash --
exactly the recorded-failure case #403 asks for.
