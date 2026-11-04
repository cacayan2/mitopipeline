# MitoPipeline (mitogenome-batch-processor)
MitoPipeline is a reproducible Snakemake/Python workflow for assembly and analysis of mitogenomes from paired-end and whole-genome sequencing reads. This pipeline exists in a broader context of ecological surveying of freshwater fish, and was designed to create a database of mitogenomes for freshwater fish (~140 samples). However, the pipeline's design focuses on reproducibility, error handling, reporting, and portability to future datasets and organism types and could be readily adapted for assemblies of different types of organisms.

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
>   + All artifacts required for publication are contained in a folder called "aggregated_data"