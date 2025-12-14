# Script for combining cruaud taxonomy data with gb_accn annotations
# with extensive sequence metadata, and then generating the
# corresponding FASTA file with sequences

source("../code/seq_fxns.R")

D <- read.delim("cruaud_taxonomy_coi_min_bc_cov_300_steps1-4.tsv")
E <- read.delim("cruaud_ncbi_coi_survey_steps1-4.tsv")

# Remove a couple of matches that are untenable
# These are matched by NCBI to "Encarsia" and "Ormyrus", respectively,
# and there are several tips in each of these genera already.
D$gb_accn[match("Dirphys",D$Genus)] <- ""
D$gb_accn[match("Ormyrulus",D$Genus)] <- ""

# Merge with sequence metadata
D <- merge(D, E, all.x=TRUE)

# Manually add Perilampus sequences
# See mail from John Heraty (this repo) where
# he indicates that Perilampus hyalinus sequences
# are on BOLD but not published to NCBI yet, and
# that Perilampus chrysopae can be added to one
# of the tips in the other Perilampus clade.
D$bold_accn <- ""
X <- read.delim("perilampus_hyalinus_bold.tsv")
idx <- which(D$TipLabel==X$TipLabel)
for (i in 1:ncol(X)) {
    D[idx,which(colnames(D)==colnames(X)[i])] <- X[,i]
}
X <- read.delim("perilampus_chrysopae_ncbi.tsv")
idx <- which(D$TipLabel==X$TipLabel)
for (i in 1:ncol(X)) {
    D[idx,which(colnames(D)==colnames(X)[i])] <- X[,i]
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


# Assemble sequence file
core_seqs <- read.FASTA("cruaud_ncbi_coi_survey_steps1-4.fasta")
perilampus_hyalinus_seq <- read.FASTA("perilampus_hyalinus_bold.fasta")
perilampus_chrysopae_seq <- read.FASTA("perilampus_chrysopae_ncbi.fasta")
all_seqs <- c(core_seqs,perilampus_hyalinus_seq,perilampus_chrysopae_seq)

# Treat either gb_accn or bold_accn (the latter only for Perilampus hyalinus)
# as "GenBank" identifiers
D$GenBank <- paste0(D$gb_accn,D$bold_accn) 
D <- D[D$GenBank!="",]
D$Cruaud_TipLabel <- D$TipLabel
D$TipLabel <- paste(D$gb_genus,D$gb_species,sep=" ")
D$TipLabel[D$Cruaud_TipLabel=="PERI_Perilampus_MERG00218_0101"] <- "Perilampus_chrysopae"
D$TipLabel[D$Cruaud_TipLabel=="PERI_Perilampus_MERG00219_0101"] <- "Perilampus_hyalinus"

# We need to filter out the relevant sequences based on
# the GenBank identifiers
a <- strsplit(names(cruaud_seqs))
accn <- character()
for (i in 1:length(a)) accn[i] <- a[[i]][1]
cruaud_seqs <- all_seqs[accn %in% D$GenBank]
cruaud_seqs <- rename_seqs(cruaud_seqs, D)

# Extract the barcode part and add info about
# the barcode sequence
D$Start <- 0
D$Stop <- 0
idx <- match(names(cruaud_seqs),D$TipLabel)
for (i in 1:length(cruaud_seqs)) {

    # The following will reverse the sequence if the alignment
    # info (aln_q_neg) is TRUE, so that it is the complement of
    # the sequence that we need. Thus, we only need to complement
    # it afterwards.
    seqs[[i]] <- seqs[[i]][D$aln_q_start[idx[i]]:D$aln_q_end[idx[i]]]
    if (D$aln_q_neg)
        seqs[[i]] <- complement(seqs[[i]])

    # The following code should ensure that the sequence start and
    # stop are set in frame
    D$Start[idx[i]] <- (3 - ((D$aln_t_start[idx[i]] + 1) %% 3)) %% 3
    D$Stop[idx[i]]  <- D$aln_t_end[idx[i]] - ((D$aln_t_end[idx[i]] + 2) %% 3)

    seqs[[i]] <- seqs[[i]][D$Start[idx[i]]:D$Stop[idx[i]]]
}

# Write the resulting taxonomy
write.table(D,"cruaud_CO1_taxonomy.tsv",row.names=FALSE,sep="\t")

# Write the resulting sequences
write.FASTA(cruaud_seqs, "cruaud_CO1.fasta")

# Generate cruaud CO1 seqs tree
tree <- read.tree("IQ_COMBINED.tre")
tip_labels <- D$Cruad_TipLabel
tree <- drop.tip(D, tip_labels)
tree$tip.labels <- D$TipLabel[match(tree$tip.labels,D$Cruaud_TipLabel)]
write.tree(tree, "cruaud_CO1.nwk")


