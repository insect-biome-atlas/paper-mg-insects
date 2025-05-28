library(ape)

# Function to extract clean family data from Chesters
extract_clean_fam_data <- function(taxonomy, tree, seqs, family, path=paste0("../",tolower(family),"/")) {

    tips <- taxonomy$TipLabel[taxonomy$Family==family]
    if (length(tips)>2)
        tree <- extract.clade(tree,getMRCA(tree,tips))
    seqs <- seqs[tips]

    start1 <- no_stop_codons(seqs, 1)
    start2 <- no_stop_codons(seqs, 2)
    start3 <- no_stop_codons(seqs, 3)
    keep <- (start1 | start2 | start3)

    start1_tips <- names(seqs)[start1]
    start2_tips <- names(seqs)[start2]
    start3_tips <- names(seqs)[start3]
    drop_tips <- names(seqs)[!keep]
    keep_tips <- names(seqs)[keep]

    if (length(tips)>2)
        tree <- drop.tip(tree, drop_tips)
    seqs <- seqs[keep]
    meta <- taxonomy[match(keep_tips, taxonomy$TipLabel),]
    meta$Start <- 3
    meta$Start[match(start2_tips,meta$TipLabel)] <- 2
    meta$Start[match(start1_tips,meta$TipLabel)] <- 1
    # Make sure length is divisible by 3
    meta$Stop <- 0
    for (i in 1:nrow(meta)) {
        len <- length(seqs[[i]])
        meta$Stop[i] <- len - (meta$Start[i] - 1)
    }

    write.FASTA(seqs,paste0(path,"chesters_",tolower(family),".fasta"))
    write.tree(tree,paste0(path,"chesters_",tolower(family),".nwk"))
    write.table(meta,row.names=FALSE,sep="\t",paste0(path,"chesters_",tolower(family),"_taxonomy.tsv"))
}

# Read in seq functions
source("../code/seq_fxns.R")

# Read data
taxonomy <- read.delim("chesters_new_outgroups_taxonomy.tsv")
tree <- read.tree("chesters_new_outgroups.nwk")
seqs <- read.FASTA("chesters_new_outgroups.fasta")

# Extract clean data
extract_clean_data(taxonomy, tree, seqs, "Ichneumonidae", "Family")
extract_clean_data(taxonomy, tree, seqs, "Braconidae", "Family")
extract_clean_data(taxonomy, tree, seqs, "Cecidomyiidae", "Family")

