# entognatha_outgroup
Data and scripts for generating a tree for Entognatha and outgroups, to replace the outgroups in the Chesters (2017) tree.

## Details
For Entognatha and outgroups, we used the data from Sundh et al (2024). Briefly, it consists of data assembled from Bellini et al (2023) on Collembola, and missing data on other Entognatha groups and on suitable outgroups from GenBank. As far as possible, in this process, data from mitochondrial genomes was favored, as the Bellini et al (2023) paper is based on mitogenomes. The data are in the files `collembola_CO1.fasta`, `missing_entognatha_CO1.fasta` and `root_mt_genomes.fasta`.

The corresponding taxonomy files are `collembola_taxonomy.tsv`, `missing_entognatha_taxonomy.tsv` and `root_taxonomy.tsv`. They contain data on the start and stop of the CO1 sequences (in frame).

The preferred tree in Bellini et al (2023; Figure 1) is given in the file collembola_tree.nwk.

During the reanalysis, it was discovered that the sequence attributed to Fujientomon (Protura) in NCBI is likely to be a fungal sequence. The sequence was removed and replaced with the CO1 sequence from the mitochondrial genome of Sinentomon erythranum, which belongs to the same order.

The order classification of Andinentulus and Yamatentomon was also corrected from Sinentomata to Acerentomata.

A bacterial sequence was added to help in the capture of bacterial contaminants, namely a sequence from the genus Wolbachia (a common insect symbiont in the alphaproteobacteria, the bacterial relatives of mitochondria).

Finally, a plant sequence was added (from the mitochondrial genome of Oryza), representing grass pollen that is often found in Malaise trap samples. Note that the location of the CO1 sequence is incorrectly specified in the GenBank record. We have corrected it in the `root_taxonomy.tsv` file.

We merged the sequence data and taxonomy from these corrected sources (in `merge_entognatha_outgroup_data.R`).

The merged sequences were aligned as amino acid sequences and the alignment converted back to a nucleotide alignment, and sites with 90% or more gaps were trimmed away (see `make_mb_data.sh` for the scripts used to achieve this). The resulting alignment is in `entognatha_outgroup_aligned_trimmed.fasta`. Note that this alignment is slightly trimmed compared to the one used in Sundh et al (2024).

For the conversion from aa to nucleotide alignment with pal2nal, we used translation table 5 throughout, which generates a number of warnings for taxa with other translation tables. However, the resulting nucleotide alignments are not affected by this, so we did not attempt to remove these warnings.

For the phylogenetic analysis, we generated partial constraints for all of the clades and taxa included in the Bellini et al (2023) analysis. The generated constraints are in the file `collembola_constraints.nex`. We also generated hard constraints for families and for well-established relationships among outgroup taxa. This is controlled in the script `generate_entognatha_outgroup_nexus_files.R`.

We ran the MrBayes analysis for 10 M generations using a strict clock model with a codon-partitioned GTR+Gamma model (see `mb_runs/run1/run.nex`).

## References
Bellini BC, Zhang F, Cavalcante de Souza PG, Clicia dos Santos-Costa R, da Silva Medeiros G, Godeiro NN (2023) The Evolution of Collembola Higher Taxa (Arthropoda, Hexapoda) Based on Mitogenome Data. Diversity 15: 7. https://doi.org/10.3390/d15010007.

For additional references, see the README in the parent folder.
