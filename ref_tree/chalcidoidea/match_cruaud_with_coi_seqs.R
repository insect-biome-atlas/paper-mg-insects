# Script for matching the cruaud et al tree with suitable coi
# reference sequences in GenBank

# Read in libraries and functions needed
source("gb_fxns.R")
source("search_chalcidoid_coi_seqs_ncbi.R")
source("select_ref_seq_fxn.R")

# Read in Cruaud et al tree
tree <- read.tree("IQ_COMBINED.tre")

# Read in base taxonomy for tips in the Cruaud et al tree
D <- read.delim("cruaud_taxonomy.tsv")

# Restrict to ingroups
D <- D[D$IN.OUT=="INGROUP",]
D$gb_accn <- ""

# Create data frame for results
# E <- data.frame()
E <- read.delim("cruaud_coi_extension_taxonomy_step2_part1.tsv")

# Create fresh file for sequences
seq_file <- "cruaud_coi_extension_step2.fasta"
cat(file=seq_file,"") # Make sure we start from scratch

# Mark the taxa that already have suitable sequences
for (i in 1:nrow(D)) {

    if (D$TipLabel[i] %in% E$TipLabel) {
        F <- E[E$TipLabel==D$TipLabel[i],]
        D$gb_accn[i] <- select_best_ref_seq(F)
    }
}

# 1. Get all species-level sequences
for (i in 1:nrow(D)) {

    next

    res <- data.frame()

    if (D$Species_epithet[i] != "sp") {
        cat("Fetching common species:",D$Species[i],"\n")
        res <- fetch_coi_sequences(D$Species[i], seq_file)
    } else {
        sp_indet <- c("sp","n.a.","sp01","sp1","sp2","nsp018","nsp035","nsp045")
        AHE_species <- D$AHE_species[i][length(D$AHE_species[i])] # Get rid of prefixes like aff, cf, nr and the like
        UCE_species <- D$UCE_species[i][length(D$UCE_species[i])] # Ditto
        if (!(AHE_species %in% sp_indet)  || !(UCE_species %in% sp_indet)) {
            species <- character()
            if (!(AHE_species %in% sp_indet))
                species <- paste(D$AHE_genus[i],AHE_species)
            if (!(UCE_species %in% sp_indet))
                species <- c(species, paste(D$UCE_genus[i],UCE_species))
            cat("Fetching unique species:",species,"\n")
            res <- fetch_coi_sequences(species, seq_file)
        }
    }

    if (nrow(res)>0) {
        res$TipLabel <- D$TipLabel[i]
        E <- rbind(E,res)
        res <- data.frame()
    }
}

# 2. Get all singleton genus sequences
for (i in 280:nrow(D)) {

    res <- data.frame()

    cat("Processing genus",genus,"\n")
    genus <- D$Genus[i]
    if (sum(D$Genus==genus)==1 && D$gb_accn[i] == "") {
        cat("Fetching missing singleton genus",genus,"\n")
        res <- fetch_coi_sequences(genus, seq_file)
    }
    if (nrow(res)>0) {
        res$TipLabel <- D$TipLabel[i]
        E <- rbind(E,res)
        res <- data.frame()
    }
}

# Mark the taxa that now have suitable sequences
for (i in 1:nrow(D)) {

    if (D$gb_accn[i]=="" && D$TipLabel[i] %in% E$TipLabel) {
        F <- E[E$TipLabel==D$TipLabel[i],]
        D$gb_accn[i] <- select_best_ref_seq(F)
    }
}

# 3. Fill in missing genera if monophyletic
multi_genera <- unique(D$Genus[duplicated(D$Genus)])
for (genus in multi_genera[15:length(multi_genera)]) {

    res <- data.frame()

    cat("Processing multigenus",genus,"\n")
    idx <- which(D$Genus==genus)
    if (sum(D$gb_accn[idx]!="")==0 && is.monophyletic(tree,D$TipLabel[idx])) {
        cat("Fetching missing monophyletic multigenus",genus,"\n")
        res <- fetch_coi_sequences(genus, seq_file)
    }
    if (nrow(res)>0) {
        res$TipLabel <- D$TipLabel[idx[1]]
        E <- rbind(E,res)
        res <- data.frame()
    }
}

write.table(D,"cruaud_coi_taxonomy_step2.tsv")
write.table(E,"cruaud_coi_extension_taxonomy_step2.tsv", sep="\t",row.names=FALSE)
