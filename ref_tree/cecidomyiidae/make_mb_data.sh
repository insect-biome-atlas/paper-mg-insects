#! /bin/bash
# Construct aligned data matrix in Nexus format for
# Cecidomyiidae analysis with MrBayes

# The script starts from data extracted from Sundh et al (2024)

# Update the Chesters taxonomy using NCBI taxonomy data
Rscript --vanilla update_chesters_cecidomyiidae_taxonomy.R

# Generate sikora data
R --no-save < generate_sikora_data.R

# Merge Cecidomyiidae data
R --no-save < merge_cecidomyiidae_data.R

# Align sequences
../code/align_seqs.sh expanded_cecidomyiidae.fasta

# Delete sites with 90% or more gaps
Rscript --vanilla ../code/trim_gaps.R expanded_cecidomyiidae_aligned.fasta expanded_cecidomyiidae_aligned_trimmed.fasta 0.9

# Make the NEXUS data and constraint files
R --no-save < generate_cecidomyiidae_nexus_files.R

