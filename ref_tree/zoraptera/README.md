# zoraptera

There is only one Zoraptera sequence in the Chesters (2017) tree, and the same species is represented by a mitochondrial genome in GenBank. Therefore, we assembled a Zoraptera dataset from scratch, downloading all available CO1 sequences (max one per species). As outgroups, we chose three diverse representatives of Dermaptera.

The taxonomic and ncbi info for the Zoraptera and outgroup data is summarized in `missing_zoraptera_ncbi.tsv`, and this info is used to fetch the corresponding sequences in FASTA format to the file `missing_zoraptera_CO1.fasta` using the script `fetch_zoraptera_seqs.R`. The same script also reformats the sequence file so that the sequences are labeled with appropriate tip labels; the output is written to the file `zoraptera.fasta`. Finally, the script augments the taxonomy information, producing the file `zoraptera_taxonomy.tsv`.

The sequences were then aligned as amino acid sequences and the alignment converted back to a nucleotide alignment, and sites with 90% or more gaps were trimmed away (see `make_mb_data.sh` for the scripts used to achieve this). The resulting alignment is in `zoraptera_aligned_trimmed.fasta`.

For the phylogenetic analysis, we generated hard constraints for the orders (Dermaptera and Zoraptera). The generated constraints are in the file `zoraptera_constraints.nex` in the `mb_runs` folder. 

We ran the MrBayes analysis for 10 M generations using a strict clock model with a codon-partitioned GTR+Gamma model (see `mb_runs/run1/run.nex`).

