#! /bin/bash

# Construct aligned data matrix and constraint
# file in Nexus format for Braconidae analysis
# with MrBayes

# Merge DNA sequences and taxonomic info
R --no-save < merge_braconidae_data.R

# Align sequences
../code/align_seqs.sh expanded_braconidae.fasta

# Delete sites with 90% or more gaps
Rscript --vanilla ../code/trim_gaps.R expanded_braconidae_aligned.fasta expanded_braconidae_aligned_trimmed.fasta 0.9

# Make the NEXUS data and constraint files
R --no-save < generate_braconidae_nexus_files.R

