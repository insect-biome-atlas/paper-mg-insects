# emb_pha_not

Embioptera, Phasmatodea and Notoptera (Mantophasmatodea + Grylloblattodea) form a clade in Misof et al. (2014). They are sparsely represented in Chesters (2017): one sequence each for Mantophasmatodea and Grylloblattodea, and two for Embioptera. Phasmatodea is more richly represented (38 sequences).

We collected sequences for Embioptera and Notoptera from GenBank with the aim to replace the sequences in Chesters (2017), using some Phasmatodea sequences as an internal outgroup or reference group (sister group of Embioptera). This allows us to replace the Embioptera and Notoptera parts of the Chesters (2017) tree in two steps, using the sister group in both cases as the outgroup.

There is only one useful CO1 sequence of Mantophasmatodea in GenBank; there has been a lot of CO1 sequencing in the group but of CO1 parts that do not overlap at all with the barcoding region. There are more sequences of the other groups. We downloaded all CO1 sequences from specimens identified to species, if they had any overlap with the barcoding region. We added up to five different unnamed species of each genus, if available. If there was a choice, we preferred complete CO1 sequences from mitochondrial genomes, otherwise as complete as possible given that there was a significant overlap with the barcoding region.

For Phasmatodea, we chose three CO1 sequences from complete mitochondrial genomes, matching the Misof et al (2014) representatives as closely as possible.

The NCBI record information was assembled in `emb_pha_not_ncbi_recs.tsv`. The script `fetch_emb_pha_not_seqs.R` fetches the fasta sequences into the file `emb_pha_not_CO1.fasta`. It then relabels the sequences with appropriate tip labels and place them in the file `emb_pha_not.fasta`. Finally, it completes the taxonomic info and places it in the file `emb_pha_not_taxonomy.tsv`.

The sequences were then aligned as amino acid sequences and the alignment converted back to a nucleotide alignment, and sites with 90% or more gaps were trimmed away (see `make_mb_data.sh` for the scripts used to achieve this). The resulting alignment is in `emb_pha_not_aligned_trimmed.fasta`.

For the phylogenetic analysis, we generated hard constraints for the orders (Embioptera, Phasmatodea, Mantophasmatodea and Grylloblattodea), as well as for the higher clades Embioptera+Phasmatodea and Notoptera (Mantophasmatodea+Grylloblattodea). The generated constraints are in the file `emb_pha_not_constraints.nex` in the `mb_runs` folder. 

We ran the MrBayes analysis for 10 M generations using a strict clock model with a codon-partitioned GTR+Gamma model (see `mb_runs/run1/run.nex`).

