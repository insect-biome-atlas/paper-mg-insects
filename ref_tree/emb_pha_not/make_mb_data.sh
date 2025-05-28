#! /bin/bash

# Construct aligned data matrix in Nexus format for
# analysis of Embioptera, Phasmatodea and Notoptera
# (Grylloblattodea + Mantophasmatodea) with MrBayes

# Fetch coding CO1 sequences and label with TipLabel
R --no-save < fetch_emb_pha_not_seqs.R

# Align sequences
../code/align_seqs.sh emb_pha_not.fasta

# Delete sites with 90% or more gaps
Rscript --vanilla ../code/trim_gaps.R emb_pha_not_aligned.fasta emb_pha_not_aligned_trimmed.fasta 0.9

# Make the NEXUS data and constraint files
R --no-save < generate_emb_pha_not_nexus_files.R

