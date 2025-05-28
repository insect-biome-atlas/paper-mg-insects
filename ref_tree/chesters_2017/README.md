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
- All others can be pasted in by the "stalk" length with respect to other insect
groups in the Chesters tree. Alternatively, we use the 

The orders with less than 10 representatives in Chesters 2017:
- Grylloblattodea 1, 7 spp in BOLD ("Grylloblattidae", part of "Notoptera")
- Mantophasmatodea 1, 2 spp in BOLD ("Mantophasmatidae", part of "Notoptera")
- Zoraptera 1, 8 spp in BOLD
- Embioptera 2, 72 spp in BOLD

## Data and scripts
The files `chesters_new_outgroups_taxonomy.tsv`, `chesters_new_outgroups.fasta` and `chesters_new_outgroups.nwk` contain the taxonomy, sequences and tree from Sundh et al (2024).

The files `misof_2014_fig1_insecta_orders.nwk` and `misof_2014_fig1_insecta_orders_dated.nwk` contain order-level backbone trees with and without branch lengths in terms of time units (Myr). These trees were hand-coded from Misof et al (2014) Fig. 1 and the median age estimates given in Fig. 2, as specified in the Supplementary Material.

The script `backbone_info.R` analyzes the Sundh et al (2024) tree in light of the Misof et al (2014) tree as specified above.

The script `generate_corrected_tree.R` will generate an Insecta tree based on the Misof et al (2014) backbone, and order-level trees based on Sundh et al (2024). The resulting tree is written to the file `s24insecta_m14backbone.nwk` and the taxonomy info to the file `s24insecta_m14backbone_taxonomy.tsv`.

The script `extract_clades.R` extracts the info relating to the families Ichneumonidae, Braconidae and Cecidomyiidae. The info is written to the directories `../ichneumonidae/`, `../braconidae/` and `../cecidomyiidae`. It uses the Sundh et al (2024) data but the result should be the same if it were to be applied to the corrected Insecta tree.

