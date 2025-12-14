library(ape)

# Function to extract family data from Chesters
# We also add start and stop (=1 and sequence length)
# positions for translation, for use in downstream
# analyses.
extract_data <- function(taxonomy, tree, seqs, family) {

    tips <- taxonomy$TipLabel[taxonomy$Family %in% family]
    if (length(tips)<=2) {
        cat("ERROR: There should be more than 3 tips; found only",length(tips),"tips\n")
        return (list())
    }
    tree <- extract.clade(tree,getMRCA(tree,tips))
    seqs <- seqs[tips]
    meta <- taxonomy[match(tips, taxonomy$TipLabel),]
    meta$Start <- 1
    meta$Stop <- 0
    for (i in 1:length(seqs)) {
        meta$Stop[i] <- length(seqs[[i]])
    }

    list(tree=tree, seqs=seqs, taxonomy=meta)
}

# Read in seq functions
source("../code/seq_fxns.R")

# Read data
taxonomy <- read.delim("chesters_new_outgroups_taxonomy_updated.tsv")
tree <- read.tree("chesters_new_outgroups_updated.nwk")
seqs <- read.FASTA("chesters_new_outgroups_updated.fasta")

# Extract and write clade-specific data
write.tsv <- function(x, file_name) { write.table(x, file=file_name, row.names=FALSE, sep="\t") }

res <- extract_data(taxonomy, tree, seqs, "Ichneumonidae")
write.FASTA(res$seqs,  "../ichneumonidae/chesters_ichneumonidae.fasta")
write.tree(res$tree,   "../ichneumonidae/chesters_ichneumonidae.nwk")
write.tsv(res$taxonomy,"../ichneumonidae/chesters_ichneumonidae_taxonomy.tsv")

res <- extract_data(taxonomy, tree, seqs, "Braconidae")
write.FASTA(res$seqs,  "../braconidae/chesters_braconidae.fasta")
write.tree(res$tree,   "../braconidae/chesters_braconidae.nwk")
write.tsv(res$taxonomy,"../braconidae/chesters_braconidae_taxonomy.tsv")

res <- extract_data(taxonomy, tree, seqs, "Cecidomyiidae")
write.FASTA(res$seqs,  "../cecidomyiidae/chesters_cecidomyiidae.fasta")
write.tree(res$tree,   "../cecidomyiidae/chesters_cecidomyiidae.nwk")
write.tsv(res$taxonomy,"../cecidomyiidae/chesters_cecidomyiidae_taxonomy.tsv")

# Define chalcidoid families
# Include Mymarommatidae outgroup
chalcidoids <- c("Pteromalidae",
                "Aphelinidae",
                "Perilampidae",
                "Eulophidae",
                "Encyrtidae",
                "Eupelmidae",
                "Trichogrammatidae",
                "Chalcididae",
                "Ormyridae",
                "Megastigmidae",
                "Eurytomidae",
                "Torymidae",
                "Leucospidae",
                "Mymaridae",
                "Mymarommatidae",
                "Eucharitidae",
                "Agaonidae",
                "Tetracampidae",
                "Rotoitidae")
res <- extract_data(taxonomy, tree, seqs, chalcidoids)
write.FASTA(res$seqs,  "../chalcidoidea/chesters_chalcidoidea.fasta")
write.tree(res$tree,   "../chalcidoidea/chesters_chalcidoidea.nwk")
write.tsv(res$taxonomy,"../chalcidoidea/chesters_chalcidoidea_taxonomy.tsv")

