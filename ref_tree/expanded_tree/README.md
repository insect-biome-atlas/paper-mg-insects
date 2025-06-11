# expanded_tree
Folder containing scripts for generating an expanded and corrected version of the Sundh et al (2024) version of the Chesters (2017) tree.

The script `generate_expanded_tree.R` fetches information from each of the expanded analyses of insect clades, and then replaces the old version of each clade in the corrected Sundh et al (2024) version of the Chesters (2017) Insecta tree (see the folder `chesters_2017`) with the new version. Finally, the Entognatha representatives and outgroups are added in (based on the data assembled in Sundh et al, 2024).

 The expanded tree is written to `chesters_expanded.nwk`. The corresponding taxonomy is written to `chesters_expanded_taxonomy.tsv` and the unaligned sequences to `chesters_expanded.fasta`.

