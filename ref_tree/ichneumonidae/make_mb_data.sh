#! /bin/bash

# Construct aligned data matrix and constraint
# file in Nexus format for Ichneumonidae analysis
# with MrBayes

# Update taxonomic info
R --no-save < update_chesters_ichneumonidae_taxonomy.R

# Merge DNA sequences and taxonomic info
R --no-save < merge_ichneumonidae_data.R

# Align sequences
../code/align_seqs.sh expanded_ichneumonidae.fasta

# Delete sites with 90% or more gaps
Rscript --vanilla ../code/trim_gaps.R expanded_ichneumonidae_aligned.fasta expanded_ichneumonidae_aligned_trimmed.fasta 0.9

# Make the NEXUS data and constraint files
R --no-save < generate_ichneumonidae_nexus_files.R

