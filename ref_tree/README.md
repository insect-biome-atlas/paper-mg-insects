# ref-tree
This repo contains all data and scripts used to construct the reference tree used for taxonomic annotation

The folders are as follows:

## `chesters_2017`
Contains the original tree from Sundh et al (2024), derived from Chesters (2017). It also contains the backbone tree from Misof et al (2014), and the scripts used to correct the original Insecta tree by removing an incorrectly placed sequence and correcting the order-level backbone so that it agrees with the tree in Misof et al (2014). The output tree has branch lengths corresponding to time (Myr) according to the divergence time estimates in Misof et al (2014) (median of medians).

## `ncbi`
Contains detailed classification for arthropod genera from the NCBI taxonomy, including subtribe, tribe, supertribe and subfamily where available. Also contains a script for regenerating these data from a dump of the NCBI taxonomy.

## `ichneumonidae`
Contains data and scripts to generate an expanded tree for Ichneumonidae.

## `braconidae`
Contains data and scripts to generate an expanded and more accurate tree for Braconidae.

## `cecidomyiidae`
Contains data and scripts to generate a considerably expanded and more accurate tree for Cecidomyiidae.

## `zoraptera`
Contains data and scripts to generate an expanded tree for Zoraptera from scratch.

## `emb_pha_not`
Contains data and scripts to generate an expanded tree for the clade Embioptera + Phasmatodea + Notoptera. Notoptera = Grylloblattodea + Mantophasmatodea.

## `entognatha_outgroup`
Contains data and scripts to generate a tree from scratch for Entognatha and outgroups. The data and scripts are the same as in Sundh et al (2024).

## `code`
Contains bash and R scripts used to generate the new trees

## `expanded_tree`
Contains a script used to replace the old subtrees in the corrected Insecta tree from the folder `chesters_2017` with the new trees generated for the target taxa above. The script then adds in the root part of the tree, comprising Entognatha representatives and outgroups.

## References
Chesters D (2017) Construction of a Species-Level Tree of Life for the Insects and Utility in Taxonomic Profiling. Systematic Biology 66: 426–39. https://doi.org/10.1093/sysbio/syw099.

Misof B, Liu S, Meusemann K, Peters RS, Donath A, Mayer C, Frandsen PB, et al. (2014) Phylogenomics Resolves the Timing and Pattern of Insect Evolution. Science 346: 763–67. https://doi.org/10.1126/science.1257570.

Sundh et al (2024) HAPP: High-Accuracy Pipeline for Processing deep metabarcoding data
Sundh J, Granqvist E, Iwaszkiewicz-Eggebrecht E, Manoharan L, van Dijk LJA, Goodsell R, et al. (2024) High-Accuracy Pipeline for Processing deep metabarcoding data. bioRxiv 2024.12.20.629441; doi: https://doi.org/10.1101/2024.12.20.629441

## More info
There is a README file in each of the folders listed above, giving more information about the contents, how it was derived, and relevant references.
