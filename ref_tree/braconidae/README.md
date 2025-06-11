# braconidae

For Braconidae, we retrieved the original subfamily classification from NCBI and updated it using information from the UCE analysis of Jasso-Martínez et al. (2022; Figs. 3-8). The updated classification is in the file `Chesters_Braconidae_taxonomy_updated.xlsx` and exported in ';'-delimited format in `Chesters_Braconidae_taxonomy_updated.csv`.

We then manually added the available CO1 sequences of the  missing subfamilies. The information about these taxa and the corresponding fasta sequences is in the files:
 - `missing_braconidae_subfamilies_CO1.fasta`
 - `missing_braconidae_subfamilies_taxonomy.tsv`

Finally, we added the available CO1 sequences of the genera included in the UCE analysis of Jasso-Martínez et al. (2022) but not present in the Chesters tree. All UCE genera and their subfamily placement is given in the tsv file below. The genera are partly manually coded from the figures and the text, and partly from the tree given in the Insect Phylogeny Database (Chesters 2025), which is also included here.
 - `uce_braconidae_genera_taxonomy.tsv`: Taxonomy of the genera included in the UCE tree (including genera in the main text figures AND in the IPD tree, see next item).
 - `jassomartinez_2022.nwk`: Tree from Insect Phylogeny Database representing the UCE results (note small discrepancies between this tree and the figures in the main text)

We first added standard 658 bp CO1 barcode sequences of the missing genera using an R script. We then manually added genera with non-standard CO1 sequences. 
 - `add_missing_braconid_genera.R`: Script checking for the availability of standard barcode sequences for missing braconid genera, generating the two files below.
 - `missing_braconidae_genera_taxonomy.tsv`: Taxonomy and sequence accession data for genera with standard barcode sequences.
 - `missing_braconidae_genera_check_needed.tsv`: Taxonomy for genera for which a manual check is needed (CO1 sequences available, but not standard barcode).
 - `missing_braconidae_genera_checked_taxonomy.tsv`: Manually retrieved taxonomy and sequence accession data for genera with nonstandard barcode sequences.

We also manually retrieved full CO1 sequence data from mitochondrial genomes from two ichneumonid representatives suitable as outgroups. The taxonomy and sequence accession info are given in the file
 - `braconidae_outgroups_taxonomy.tsv`

The info in the taxonomy files is used to fetch the actual sequence data from GenBank. The taxonomy and sequence files are then merged.
 - `fetch_braconidae_sequences.R`: Script for fetching Braconidae sequences from GenBank.
 - `missing_braconidae_data.R`: Script merging sequence and taxonomy files, producing the `expanded_braconidae.fasta` and `expanded_braconidae_taxonomy.tsv` files.

The sequences were then aligned as amino acid sequences and the alignment converted back to a nucleotide alignment, and sites with 90% or more gaps were trimmed away (see `make_mb_data.sh` for the scripts used to achieve this). The resulting alignment is in `expanded_braconidae_aligned_trimmed.fasta`.

For the phylogenetic analysis, we used all higher clades (subfamily or above) with more than 95% bootstrap support in the UCE analysis as partial constraints. We added hard constraints for the families. The data and constraint files are generated in the script `generate_braconidae_nexus_files.R`, and the output files are `braconidae.nex` and `braconidae_constraints.nex`. They are placed in the `mb_runs` folder.

We ran the MrBayes analysis for 100 M generations using a strict clock model with a codon-partitioned GTR+Gamma model (see `mb_runs/run2/run.nex`).

