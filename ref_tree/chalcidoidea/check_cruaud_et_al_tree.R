# Check monophyly for duplicated genera in Cruaud et al tree

library(ape)

tree <- read.tree("IQ_COMBINED.tre")

D <- read.delim("cruaud_taxonomy.tsv")
dup_genus <- unique(D$Genus[duplicated(D$Genus) & !D$Genus=="Genus"])

for (genus in dup_genus) {

    tip_labels <- D$TipLabel[D$Genus==genus]
    if (!is.monophyletic(tree,tip_labels))
        cat ("Genus '",genus,"' is not monophyletic\n", sep="")
}

