#! /bin/bash

# The script assumes that the alignments directory exists

# The first command line argument ($1) should be the original unaligned fasta file

# Create root name
Name=$(echo $1 | cut -d "." -f 1)

# Read coding DNA sequences and translate to aa sequences
Rscript --vanilla ../code/generate_root_aa_seqs.R $1 ${Name}_taxonomy.tsv ./alignments/${Name}_aa.fasta

# Align aa sequences using mafft --auto
mafft --auto ./alignments/${Name}_aa.fasta > ./alignments/${Name}_aa_aligned.fasta

# Convert them back to nucleotide sequences

# pal2nal will correctly report a few errors in the handling of ambiguities. This should have negligible effect on
# the analysis and is ignored here.

pal2nal.pl ./alignments/${Name}_aa_aligned.fasta $1 -codontable 5 -output fasta > ${Name}_aligned.fasta

