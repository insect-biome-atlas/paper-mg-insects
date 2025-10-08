# chesters_2017

## Background
This directory contains files and scripts related to the Chesters (2017) tree, as modified in Sundh et al (2024).

As pointed out by us previously (Sundh et al. 2024), the original tree
published by Chesters did not retain the monohyly of Insecta with
respect to the outgroups. After replacing the outgroup part of the
tree, and further quality filtering of the data as decribed in
Sundh et al (2024), we discovered additional problems with the backbone
of the Chesters tree.

First, we discovered that Orthoptera and Diptera were not monophyletic. This
is due to the fact that Lilaea_fuliginosa, a tabanid, is placed in Orthoptera.
It is an obvious error and this species should be removed from the tree. When this
is done, Orthoptera and Diptera both become monophyletic.

Second, even after this fix is introduced, it is evident that many of the higher
clades in the preferred backbone tree (Fig. 3 in the paper) are not
retained in the published tree. This is true for the order Hemiptera, as well as
for a number of higher interorder clades, especially in the Exopterygota.
The text also conflicts with the figure in that the text suggests monophyly of
Mantodea and Isoptera, while the depicted tree (Fig. 3) groups Blattodea and
Isoptera. The name "Condylognatha" is applied incorrectly in the text of the paper;
it refers to (Thysanoptera + Hemiptera) and not to (Mantodea + Isoptera).

The script `backbone_info.R` reads in the modified Chesters tree from Sundh
et al (2024), removes Lilaea_fuliginosa, and then systematically tests the
monohyly of the orders and higher clades in the preferred backbone tree shown in
Chesters (2017; Fig. 3).

The results are summarized below.
- Neuropterida is monophyletic.
- Psocodea is monophyletic. Psocoptera without Phthiraptera is not (as expected).
- Hemiptera is not monophyletic. However, Homoptera, Heteroptera (ex Fulgoromorpha)
and Fulgoromorpha are all monophyletic. If hemipterans are extracted from the tree, the relationship is (Homoptera,(Fulgoromorpha,Heteroptera)),
which makes sense.

With respect to fixing these problems, we made the following observations:
- The extracted Hemiptera subtree from Chesters can be pasted in by matching the root
with the root of the Hemiptera in the Misof tree.
- The Grylloblattodea + Mantophasmatodea can be pasted in by similar root matching.
- All others can be pasted in by enforcing the "stalk" length with respect to other insect
groups in the Chesters tree to remain the same after the clade has been pasted in.

The orders with less than 10 representatives in Chesters 2017:
- Grylloblattodea 1, 7 spp in BOLD ("Grylloblattidae", part of "Notoptera")
- Mantophasmatodea 1, 2 spp in BOLD ("Mantophasmatidae", part of "Notoptera")
- Zoraptera 1, 8 spp in BOLD
- Embioptera 2, 72 spp in BOLD
These should be expanded to facilitate correct placement in the backbone tree, and increase the chances of picking up barcodes of these taxa.

## Taxonomic annotation errors
In analyzing the original taxonomic annotations from the Sundh et al (2024) paper, we discovered the following mismatches to the current NCBI taxonomy at the family level. These mismatches must be corrected to match the family names in the trait data files, which follow the current NCBI taxonomy.
- Xylophagaidae should be corrected to Xylophagidae
- Pemphigidae is now treated as a subfamily within Aphididae
- Kerriidae is now called Tachardiidae

Also, it was discovered that the Cynipoidea families were not annotated according to the revised classification of the superfamily (see Hearn et al. 2024).
- Cecinothofagus is now placed in Paraulacidae
- Diplolepis and Liebelia are now placed in Diplolepididae

We created a script correcting these family annotations in the taxonomic annotation data file from the Sundh et al (2024) paper.

We also discovered that the following families are missing from the NCBI classification, but we have trait data for them so this causes no problem in the analysis.
- Projapygidae 
- Anajapygidae

Note also that the NCBI classification treats Diplura as an order, but it is now elevated to a Class. We are using the updated classification, as in Sundh et al. (2024).

## Data and scripts
The files `chesters_new_outgroups_taxonomy.tsv`, `chesters_new_outgroups.fasta` and `chesters_new_outgroups.nwk` contain the taxonomy, sequences and tree from Sundh et al (2024). Note, however, that the family annotation has been corrected for the four cases mentioned above.

The script `correct_taxonomic_annotations.R` will correct the annotations (see above) in the `chesters_new_outgroups_taxonomy.tsv` file to create the updated file `chesters_new_outgroups_taxonomy_updated.tsv`.

The files `misof_2014_fig1_insecta_orders.nwk` and `misof_2014_fig1_insecta_orders_dated.nwk` contain order-level backbone trees with and without branch lengths in terms of time units (Myr). These trees were hand-coded from Misof et al (2014) Fig. 1 and the median age estimates given in Fig. 2, as specified in the Supplementary Material (see file `misof_2014_node_age_medians.tsv`).

The script `backbone_info.R` analyzes the Sundh et al (2024) tree in light of the Misof et al (2014) tree as specified above.

The script `generate_corrected_tree.R` will generate an Insecta tree based on the Misof et al (2014) backbone, and order-level trees based on Sundh et al (2024). The resulting tree is written to the file `s24insecta_m14backbone.nwk` and the taxonomy info to the file `s24insecta_m14backbone_taxonomy.tsv`.

The script `extract_clades.R` extracts the info relating to the families Ichneumonidae, Braconidae and Cecidomyiidae. The info is written to the directories `../ichneumonidae/`, `../braconidae/` and `../cecidomyiidae`. It uses the Sundh et al (2024) data but the result should be the same if it were to be applied to the corrected Insecta tree.

