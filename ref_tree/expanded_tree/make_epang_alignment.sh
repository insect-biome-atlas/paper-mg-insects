#! /bin/bash

# Construct aligned and trimmed data matrix for use
# by epang

# Align sequences
../code/align_seqs.sh chesters_expanded.fasta

# Delete sites with 90% or more gaps
Rscript --vanilla ../code/trim_gaps.R chesters_expanded_aligned.fasta chesters_expanded_aligned_trimmed.fasta 0.1

