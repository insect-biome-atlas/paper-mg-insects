# ichneumonidae

The taxonomic info for the Ichneumonidae in the Sundh et al (2024) tree was first complemented with subfamily, tribe and subtribe ranks from the NCBI taxonomy, and this information was then corrected or updated to reflect current classification. This was accomplished using the script `update_chesters_ichneumonidae_taxonomy.R`, generating the file `chesters_ichneumonidae_taxonomy_updated.tsv`.

CO1 sequences for missing ichneumonid taxa that would be valuable to have in the analysis because they were likely to be encountered in Madagascar but were missing from the original Chesters (2017) tree were assembled from GenBank. These data are given in the files `missing_ichneumonidae_CO1.fasta` and `missing_ichneumonidae_taxonomy.tsv`. The files also include a couple of braconid outgroups to root the tree.

The script `merge_ichneumonidae_data` was then used to merge the sequence and taxonomy data together into the files `expanded_ichneumonidae.fasta` and `expanded_ichneumonidae_taxonomy.tsv`.

The sequences were then aligned as amino acid sequences and the alignment converted back to a nucleotide alignment, and sites with 90% or more gaps were trimmed away (see `make_mb_data.sh` for the scripts used to achieve this). The resulting alignment is in `expanded_ichneumonidae_aligned_trimmed.fasta`.

For the phylogenetic analysis, as there is no comprehensive analysis of ichneumonid relationships suitable for our purposes yet, we generated partial constraints for all nodes in the ichneumonid part of the Sundh et al. (2024) version of the Chesters (2017) tree. We then defined hard constraints for the families (Ichneumonidae and Braconidae). This allowed the added ichneumonid taxa to float around in the tree found by Chesters (2017), finding their most likely placement in that tree.

We ran the MrBayes analysis for 10 M generations using a strict clock model with a codon-partitioned GTR+Gamma model (see `mb_runs/run1/run.nex`).

