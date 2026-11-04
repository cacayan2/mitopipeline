# MitoPipeline (mitogenome-batch-processor)
MitoPipeline is a reproducible Snakemake/Python workflow for assembly and analysis of mitogenomes from paired-end and whole-genome sequencing reads. This pipeline exists in a broader context of ecological surveying of freshwater fish, and was designed to create a database of mitogenomes for freshwater fish (~140 samples). However, the pipeline's design focuses on reproducibility, error handling, reporting, and portability to future datasets and organism types and could be readily adapted for assemblies of different types of organisms.

# Table of Contents
- [Pipeline Overview](#pipeline-overview)
- [Dependencies](#dependencies)
- [Installation](#installation)
- [Setup](#setup)
    + [Config](#config)
    + [Manifest](#manifest)
- [Running the Pipeline](#running-the-pipeline)
- [Background](#background)
- [v0.2.0 Changelist](#v020-changelist)

# Pipeline Overview



# Dependencies

# Installation

# Setup

## Config

## Manifest

# Running the Pipeline

# Background
The goal of this project is the development of a pipeline which will batch process and assemble raw reads from samples of fish tissue into completed mitogenomes that comprise a reference database for freshwater fish species in Illinois. This reference database would ideally then be published to NCBI GenBank and compiled into a microbiology resource announcement. 

# v0.2.0 Changelist
> **v0.2.0:** This version represents a continuation of the work started in v0.1.0 and is a complete overhaul of the original project design. While the core functionality and deliverables are preserved, when v0.1.0 was executed on sequencing data: ~40% of the data assembled, with less than 10% reaching the end of the pipeline before crashing. Changelist:
> - **Complete Metadata Collection for BioProject, Biosample, SRA, and GenBank**
>   + The pipeline is able to aggregate metadata in order to create the required submittables for NCBI data archival.
>   + Hints for how to add missing metadata are generated at the end of a pipeline run.
> - **Improved Reporting Capability**
>   + More aesthetically pleasing reports that are easier and color-coded to identify problems with the job.
>   + A more dynamic methods section with robust metadata collection protocols.
>   + Added capability to run fastQC on reads recruited by GetOrganelle to help diagnose troublesome data.
>   + Reports also include specs from the machine that ran the pipeline.
> - **Better Repository Organization**
>   + Files are organized in folders by step number for easier maintenance and extension.
>   + All files and folders less ambiguously and generous in-code commenting/docstrings.
>   + All artifacts required for publication are contained in a folder called "aggregated_data" in the output folder.
> - **Seed Selection**
>   + Fully implemented a canonical fallback sequence for biomarkers based on coverage to use for read recruitment in assembly.
>   + Supports assembly rescue capability.
> - **BLAST Phylogeny Inference**
>   + A local BLAST database is constructed to create phylogenies. 
>   + This database contains entries for *Actinopterygii* only to save space (can be changed in config).
> - **Improved Rescue Capability**
>   + <ins>Original Design:</ins> Say a hypothetical assembly would fail. v0.1.0 would take reads recruited by GetOrganelle and run them through MiniMap2 to obtain additional sequences from the original pool of sequences. GetOrganelle would then run on the originally recruited sequences and the new ones obtained by Minimap2 - this would not protect against (but instead reinforced) gaps in the genome. 
>   + <ins>New Design:</ins> Three different strategies are employed for rescue, done sequentially until a full genome is obtained or all fail. The first is the original design - this will help boost coverage of existing reads and help cut out noise. The second changes which seed is utilized, ensuring the seed isn't necessarily the issue with assembly. The third reads the full pool of reads again with the original seed used for the assembly but widens the k-mer list, raises target coverage, and allows more assembly rounds - assists with targeting complex or low-coverage regions.
> - **Improved Testing Capability**
>   + Test were written to assist with detecting edge cases and failures to detect regression when coding new sections of the pipeline and to prevent debugging on actual pipeline runs (which would take hours per run).
>   + As per Dr. Burns' advice, fixtures for quick testing were produced using Gargammel. These are small enough to upload to GitHub and can be used for a quick test of the pipeline before running on real data. 
>   + Test fixtures are simulated Illumina NGS sequencing reads of a completed *Danio rerio* mitogenome. Three different simulated datasets were produced: a complete genome, a genome with gaps in the middle, and one that is so broken that no assembly would be produced, ensuring the pipeline runs despite problematic datasets. 
> - **Robust Error Handling and Logging**
>   + Visbility for crashes and stages in the pipeline was difficult to track in v0.1.0. Logging provides both centralized (in a master log file) and decentralized (in log files designated to samples/steps) formatting allowing for thorough debugging of problematic reads.
>   + Errors and missing data are handled more intelligently preventing silent failures from crashing downstream steps.
