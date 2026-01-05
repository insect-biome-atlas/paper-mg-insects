#! /bin/bash
# Construct aligned data matrices in Nexus format
# for Chalcidoidea analyses with MrBayes

# Merge taxonomies and sequence data
R --no-save < merge_chalcidoidea_data.R

# Align sequences
../code/align_seqs.sh expanded_chalcidoidea1.fasta
../code/align_seqs.sh expanded_chalcidoidea2.fasta
../code/align_seqs.sh expanded_chalcidoidea3.fasta

# Make the NEXUS data and constraint files
R --no-save < generate_chalcidoidea_nexus_files.R

