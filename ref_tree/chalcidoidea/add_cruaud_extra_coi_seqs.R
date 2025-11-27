# Script for adding some coi sequences that were missed
# in the first round

# Read in libraries and functions needed
source("gb_fxns.R")
source("search_chalcidoid_coi_seqs_ncbi.R")
source("select_ref_seq_fxn.R")

# Read in Cruaud et al tree
tree <- read.tree("IQ_COMBINED.tre")

# Read in metadata on sequences that have been found so far
D <- read.delim("cruaud_coi_extension_min_bc_overlap_300_seq_meta.tsv")

# Read in data frame for previous results
P <- read.delim("cruaud_coi_extension_taxonomy_step2.tsv")

# Create data fram for new results
E <- data.frame()

# Create fresh file for sequences
seq_file <- "cruaud_coi_extension_step3.fasta"
cat(file=seq_file,"") # Make sure we start from scratch

# 0. Eplore all sequences of taxa with >=50 hits and no previously found COI sequences
# Identified manually based on previous results, so we now the taxa and that the relevant
# search is a genus-level search
taxa <- c("Tetrapus","Dirphys","Otitesella")
for (taxon in taxa) {

    idx <- which(D$Genus==taxon)

    cat("Fetching all sequences for TipLabel:",D$TipLabel[idx],"and taxon", taxon, "\n")
    res <- fetch_all_coi_sequences(taxon, seq_file, P)

    cat("Retrieved",nrow(res),"sequences\n")

    if (nrow(res)>0) {
        res$TipLabel <- D$TipLabel[idx]
        E <- rbind(E,res)
        res <- data.frame()
    }
}

# Mark taxa that now have suitable sequences
for (i in 1:nrow(D)) {

    if (D$gb_accn[i]=="" && D$TipLabel[i] %in% E$TipLabel) {
        F <- E[E$TipLabel==D$TipLabel[i],]
        D$gb_accn[i] <- select_best_ref_seq(F,300)
    }
}

# 1. Get all species-level sequences matching "cytochrome oxidase" [but not COI or COX]
for (i in 1:nrow(D)) {

    if (D$gb_accn[i]!="")
        next

    res <- data.frame()

    if (D$Species_epithet[i] != "sp") {
        cat("Fetching cytochrome oxidase matches of common species:",D$Species[i],"\n")
        res <- fetch_cytochrome_oxidase_sequences(D$Species[i], seq_file, P)
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
            cat("Fetching cytochrome oxidase matches of unique species:",species,"\n")
            res <- fetch_cytochrome_oxidase_sequences(species, seq_file, P)
        }
    }

    cat("Retrieved",nrow(res),"sequences\n")

    if (nrow(res)>0) {
        res$TipLabel <- D$TipLabel[i]
        E <- rbind(E,res)
        res <- data.frame()
    }
}

# 2. Get all singleton genus sequences
for (i in 1:nrow(D)) {

    if (D$gb_accn[i]!="")
        next

    res <- data.frame()

    genus <- D$Genus[i]
    if (sum(D$Genus==genus)==1 && D$gb_accn[i] == "") {
        cat("Fetching missing singleton genus",genus,"\n")
        res <- fetch_cytochrome_oxidase_sequences(genus, seq_file, P)
    }

    cat("Retrieved",nrow(res),"sequences\n")

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
for (genus in multi_genera[1:length(multi_genera)]) {

    res <- data.frame()

    cat("Processing multigenus",genus,"\n")
    idx <- which(D$Genus==genus)
    if (sum(D$gb_accn[idx]!="")==0 && is.monophyletic(tree,D$TipLabel[idx])) {
        cat("Fetching missing monophyletic multigenus",genus,"\n")
        res <- fetch_cytochrome_oxidase_sequences(genus, seq_file, P)
    }

    cat("Retrieved",nrow(res),"sequences\n")

    if (nrow(res)>0) {
        res$TipLabel <- D$TipLabel[idx[1]]
        E <- rbind(E,res)
        res <- data.frame()
    }
}

# Mark taxa that now have suitable sequences
for (i in 1:nrow(D)) {

    if (D$gb_accn[i]=="" && D$TipLabel[i] %in% E$TipLabel) {
        F <- E[E$TipLabel==D$TipLabel[i],]
        D$gb_accn[i] <- select_best_ref_seq(F,300)
    }
}

write.table(D,"cruaud_coi_extension_min_bc_overlap_300_seq_meta_step3.tsv")
write.table(E,"cruaud_coi_extension_taxonomy_step3.tsv", sep="\t",row.names=FALSE)

