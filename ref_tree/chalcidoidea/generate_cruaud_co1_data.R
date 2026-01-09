# Script for combining cruaud taxonomy data with gb_accn annotations
# with extensive sequence metadata, and then generating the
# corresponding FASTA file with sequences

source("../code/seq_fxns.R")
source("../code/fetch_fasta_fxn.R")

# Combine taxonomy and sequence metadata
# ======================================
D <- read.delim("cruaud_taxonomy_coi_min_bc_cov_300_steps1-4.tsv")
E <- read.delim("cruaud_ncbi_coi_survey_steps1-4.tsv")

# Remove a couple of matches that are untenable
# These are matched by NCBI to "Encarsia" and "Ormyrus", respectively,
# and there are several tips in each of these genera already.
D$gb_accn[match("Dirphys",D$Genus)] <- ""
D$gb_accn[match("Ormyrulus",D$Genus)] <- ""

# Remove an erroneous sequence (a Diptera sequence)
D$gb_accn[match("Cyrtogaster",D$Genus)] <- ""

# Merge with sequence metadata
D <- merge(D, E, all.x=TRUE)

# Capture species names when Cruaud taxonomy is inconclusive
for (i in 1:nrow(D)) {
    if (D$Species_epithet[i]=="sp" && !(D$gb_species[i] %in% c("","sp")))
        D$Species[i] <- paste(D$gb_genus[i],D$gb_species[i],sep=" ")
}

# Manually add Perilampus sequences
# John Heraty indicates that Perilampus hyalinus sequences
# are on BOLD but not published to NCBI yet, and
# that Perilampus chrysopae can be added to one
# of the tips in the other Perilampus clade.
D$bold_accn <- ""
X1 <- read.delim("perilampus_hyalinus_bold.tsv")
idx <- which(D$TipLabel==X1$TipLabel)
for (i in 1:ncol(X1)) {
    D[idx,which(colnames(D)==colnames(X1)[i])] <- X1[,i]
}
D$bold_accn[idx] <- X1$aln_query
X2 <- read.delim("perilampus_chrysopae_ncbi.tsv")
idx <- which(D$TipLabel==X2$TipLabel)
for (i in 1:ncol(X2)) {
    D[idx,which(colnames(D)==colnames(X2)[i])] <- X2[,i]
}

# Add metadata on barcode coverage
D$aln_chalc_bc_start <- D$aln_t_start
D$aln_chalc_bc_end <- D$aln_t_end
D$aln_chalc_bc_cov <- D$aln_t_end

idx <- which(!is.na(D$aln_query))
D$aln_chalc_bc_end[idx] <- 412 - (652 - D$aln_t_end[idx])
for (i in 1:length(idx)) { D$aln_chalc_bc_start[idx[i]] <- max(1, D$aln_t_start[idx[i]] - 240) }
D$aln_chalc_bc_cov[idx] <- D$aln_chalc_bc_end[idx] - D$aln_chalc_bc_start[idx] + 1

# Write the table for future reference
write.table(D,"cruaud_taxonomy_coi_min_bc_cov_300_seq_meta.tsv", sep="\t", row.names=FALSE)


# Assemble prefetched sequences
# =============================

core_seqs <- read.FASTA("cruaud_ncbi_coi_survey_steps1-4.fasta")
perilampus_hyalinus_seq <- read.FASTA("perilampus_hyalinus_bold.fasta")
perilampus_chrysopae_seq <- read.FASTA("perilampus_chrysopae_ncbi.fasta")
all_seqs <- c(core_seqs,perilampus_hyalinus_seq,perilampus_chrysopae_seq)

# Treat either gb_accn or bold_accn (the latter only for Perilampus hyalinus)
# as "GenBank" identifiers
D$GenBank <- paste0(D$gb_accn,D$bold_accn) 
D <- D[D$GenBank!="",]
D$Cruaud_TipLabel <- D$TipLabel
D$TipLabel <- paste(D$gb_genus,D$gb_species,sep="_")
D$TipLabel[D$Cruaud_TipLabel=="PERI_Perilampus_MERG00218_0101"] <- "Perilampus_chrysopae"
D$TipLabel[D$Cruaud_TipLabel=="PERI_Perilampus_MERG00219_0101"] <- "Perilampus_hyalinus"

# We need to filter out the relevant sequences based on
# the GenBank identifiers from all the prefetched seqs
a <- strsplit(names(all_seqs), split=" ")
accn <- character()
for (i in 1:length(a)) accn[i] <- a[[i]][1]
seqs <- all_seqs[accn %in% D$GenBank]
seqs <- seqs[!duplicated(names(seqs))]

# Rename sequences with tip labels
seqs <- rename_seqs(seqs, D)

# Extract the barcode part and add info about
# the barcode sequence
D$Start <- 0
D$Stop <- 0
idx <- match(names(seqs),D$TipLabel)
for (i in 1:length(seqs)) {

    # The following will reverse the sequence if the alignment
    # info (aln_q_neg) is TRUE, so that it is the complement of
    # the sequence that we need. Thus, we only need to complement
    # it afterwards.
    if (D$aln_q_neg[idx[i]]==FALSE)
        seqs[[i]] <- seqs[[i]][D$aln_q_start[idx[i]]:D$aln_q_end[idx[i]]]
    else {
        seqs[[i]] <- seqs[[i]][D$aln_q_end[idx[i]]:D$aln_q_start[idx[i]]]
        seqs[[i]] <- as.raw(complement(seqs[[i]]))
    }

#    cat("For tip label:",D$TipLabel[idx[i]],"length is",length(seqs[[i]]),"\n")

    # The following code should ensure that the sequence start and
    # stop are set in frame
    start <- 1 + ((3 - ((D$aln_t_start[idx[i]] + 1) %% 3)) %% 3)
    stop  <- length(seqs[[i]]) - ((D$aln_t_end[idx[i]] + 2) %% 3)

#    cat("start:",start," -- stop:",stop,"\n")

    seqs[[i]] <- seqs[[i]][start:stop]

    # Set Start and Stop after extraction of coding part
    D$Start[idx[i]] <- 1
    D$Stop[idx[i]] <- length(seqs[[i]])
}

# Write final sequence file version
write.FASTA(seqs,"cruaud_CO1.fasta")


# Generate cruaud CO1 seqs tree
# =============================

tree <- read.tree("IQ_COMBINED.tre")
tip_labels <- D$Cruaud_TipLabel
tree <- keep.tip(tree, tip_labels)
tree$tip.label <- D$TipLabel[match(tree$tip.label,D$Cruaud_TipLabel)]
tree <- root(tree,D$TipLabel[D$Familycode=="MYMO"],r=TRUE)
write.tree(tree, "cruaud_CO1.nwk")


# Update taxonomy
# ===============

# Read in manually curated annotations
T <- read.delim("cruaud_taxonomy_coi_min_bc_cov_300_seq_meta_updated.csv", sep=";")

# Start with base assumptions
D$Family <- D$OldFamily
D$Subfamily <- D$OldSubfamily
D$Tribe <- D$OldTribe
D$Clade <- D$Family

# Update with manually updated annotations
idx <- match(D$Cruaud_TipLabel,T$TipLabel)
for (i in 1:nrow(D)) {
    if (is.na(idx[i]))
        next
    if (T$TentativeNewFamily[idx[i]]!="") D$Family[i] <- T$TentativeNewFamily[idx[i]]
    if (T$TentativeNewSubfamily[idx[i]]!="") D$Subfamily[i] <- D$TentativeNewSubfamily[idx[i]]
    if (T$TentativeNewTribe[idx[i]]!="") D$Tribe[i] <- D$TentativeNewTribe[idx[i]]
    if (T$Clade[idx[i]]!="") D$Clade[i] <- T$Clade[idx[i]]
}


# Write relevant parts of the resulting taxonomy
D <- D[,c("TipLabel","Family","Subfamily","Tribe","Genus","Species","Clade","Start","Stop")]
write.table(D,"cruaud_CO1_taxonomy.tsv",row.names=FALSE,sep="\t")

