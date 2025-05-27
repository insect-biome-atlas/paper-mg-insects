# cecidomyiidae

For Cecidomyiidae, we first updated the Sundh et al (2024) taxonomy information with subfamily and supertribe information from NCBI taxonomy (in `update_chesters_taxonomy.R`).

We then retrieved the CO1 sequences from Sikora et al (2014), including outgroups (in `generate_sikora_data.R`).

Finally, we added CO1 sequences for some cecidomyiid taxa that were still poorly represented (files `missing_cecidomyiidae_CO1.fasta` and `missing_cecidomyiidae_taxonomy.tsv`).

We then merged the sequence data and taxonomy from these sources (in `merge_cecidomyiidae_data.R`). In some cases, we replaced sequences in Sundh et al (2024) with longer sequences of the same taxa that we retrieved from Sikora et al (2014) or from GenBank.

The merged sequences were aligned as amino acid sequences and the alignment converted back to a nucleotide alignment, and sites with 90% or more gaps were trimmed away (see `make_mb_data.sh` for the scripts used to achieve this). The resulting alignment is in `expanded_cecidomyiidae_aligned_trimmed.fasta`.

For the phylogenetic analysis, we constrained all clades with > PP 95% in the Bayesian analysis (their Fig. 1C) using soft constraints (in `generate_cecidomyiidae_nexus_files.R`). We then added hard constraints for the families included in the analysis, leaving the three taxa that were _incertae sedis_ in Sikora et al (2014).

We ran the MrBayes analysis for 10 M generations using a strict clock model with a codon-partitioned GTR+Gamma model (see `mb_runs/run1/run.nex`).

