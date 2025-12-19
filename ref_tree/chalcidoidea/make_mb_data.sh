#! /bin/bash
# Construct aligned data matrices in Nexus format
# for Chalcidoidea analyses with MrBayes

# Merge taxonomies and sequence data
R --no-save < merge_chalcidoidea_data.R

# Align sequences
../code/align_root_seqs.sh expanded_chalcidoidea1.fasta
../code/align_root_seqs.sh expanded_chalcidoidea2.fasta
../code/align_root_seqs.sh expanded_chalcidoidea3.fasta

# Delete sites with 90% or more gaps
Rscript --vanilla ../code/trim_gaps.R expanded_chalcidoidea1_aligned.fasta expanded_chalcidoidea1_aligned_trimmed.fasta 0.9
Rscript --vanilla ../code/trim_gaps.R expanded_chalcidoidea2_aligned.fasta expanded_chalcidoidea2_aligned_trimmed.fasta 0.9
Rscript --vanilla ../code/trim_gaps.R expanded_chalcidoidea3_aligned.fasta expanded_chalcidoidea3_aligned_trimmed.fasta 0.9

# Make the NEXUS data and constraint files
R --no-save < generate_chalcidoidea_nexus_files.R

