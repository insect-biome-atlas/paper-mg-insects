#! /bin/bash
# Construct aligned data matrix in Nexus format for
# Entognatha and outgroup analysis with MrBayes

# Merge taxonomies and sequence data
R --no-save < merge_entognatha_outgroup_data.R

# Align sequences
../code/align_seqs.sh entognatha_outgroup.fasta

# Delete sites with 90% or more gaps
Rscript --vanilla ../code/trim_gaps.R entognatha_outgroup_aligned.fasta entognatha_outgroup_aligned_trimmed.fasta 0.9

# Make the NEXUS data and constraint files
R --no-save < generate_entognatha_outgroup_nexus_files.R

