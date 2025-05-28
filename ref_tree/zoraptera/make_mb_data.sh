#! /bin/bash

# Construct aligned data matrix in Nexus format for
# Zoraptera analysis with MrBayes

# Fetch coding CO1 sequences and label with TipLabel
R --no-save < fetch_zoraptera_seqs.R

# Align sequences
../code/align_seqs.sh zoraptera.fasta

# Delete sites with 90% or more gaps
Rscript --vanilla ../code/trim_gaps.R zoraptera_aligned.fasta zoraptera_aligned_trimmed.fasta 0.9

# Make the NEXUS data and constraint files
R --no-save < generate_zoraptera_nexus_files.R

