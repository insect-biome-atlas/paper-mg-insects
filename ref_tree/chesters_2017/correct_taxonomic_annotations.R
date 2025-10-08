# Correct the taxonomic annotations in the modified Chesters tree
# from Sundh et al. (2024)

D <- read.delim("chesters_new_outgroups_taxonomy.tsv")

# Correct family names
D$Family[D$Family=="Xylophagaidae"] <- "Xylophagidae"
D$Family[D$Family=="Pemphigidae"] <- "Aphididae"
D$Family[D$Family=="Kerriidae"] <- "Tachardiidae"

# Update family names for Cynipoidea
D$Family[D$Genus=="Diplolepis"] <- "Diplolepididae"
D$Family[D$Genus=="Liebelia"] <- "Diplolepididae"
D$Family[D$Genus=="Cecinothofagus"] <- "Paraulacidae"

write.table(D, "chesters_new_outgroups_taxonomy_updated.tsv", row.names=FALSE, sep="\t")

