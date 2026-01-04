library(ape)
source("../code/seq_fxns.R")


# 1. Merge taxonomy information
# =============================

# Read in taxonomies
T1 <- read.delim("chesters_cecidomyiidae_taxonomy_detailed.tsv")
T2 <- read.delim("missing_cecidomyiidae_taxonomy.tsv")
T3 <- read.delim("sikora_taxonomy.tsv")

# Add tip label and higher ranks to missing taxa taxonomy file (includes outgroups)
T2$TipLabel<-sub(" ","_",T2$Species)
T2$Kingdom <- "Animalia"
T2$Phylum <- "Arthropoda"
T2$Class <- "Insecta"
T2$Order <- "Diptera"

# Ditto for Sikora file
T3$TipLabel<-sub(" ","_",T3$Species)
T3$Kingdom <- "Animalia"
T3$Phylum <- "Arthropoda"
T3$Class <- "Insecta"
T3$Order <- "Diptera"

# Merge into one taxonomy file
cols <- c("TipLabel","Kingdom","Phylum","Class","Order","Family","Subfamily","Supertribe","Genus","Species")

T <- rbind(T1[,cols],T2[,cols],T3[,cols])
T <- T[!duplicated(T$TipLabel),]    # Need to remove info for duplicated tip labels (the Chesters seqs are removed, the taxonomy info is the same)

# Add info on Clade, matching life history data file
T$Clade <- T$Family
for (i in 1:nrow(T)) {
    if (!is.na(T$Subfamily[i]) && T$Subfamily[i]!="") {
        T$Clade[i] <- T$Subfamily[i]
    }
}

cols <- c("TipLabel","Kingdom","Phylum","Class","Order","Family","Subfamily","Supertribe","Clade","Genus","Species")
T <- T[,cols]
write.table(T,"expanded_cecidomyiidae_taxonomy.tsv",sep="\t",row.names=FALSE)


# 2. Merge sequences
# ==================

# Read in sequences
chesters_seqs <- read.FASTA("chesters_cecidomyiidae.fasta")
extra_seqs <- read.FASTA("missing_cecidomyiidae_CO1.fasta")
sikora_seqs <- read.FASTA("sikora_CO1.fasta")

# Use more meaningful names
chesters_meta <- T1
extra_meta <- T2
sikora_meta <- T3

# Relabel missing sequences from NCBI strings to tip labels
extra_seqs <- rename_seqs(extra_seqs, extra_meta)
sikora_seqs <- rename_seqs(sikora_seqs, sikora_meta)

# Remove Chesters sequences that are also in Sikora et al, but more complete there
chesters_seqs <- chesters_seqs[!(names(chesters_seqs) %in% sikora_meta$TipLabel)]
chesters_meta <- chesters_meta[!(chesters_meta$TipLabel %in% sikora_meta$TipLabel),]

# Remove Chesters sequences that are also in missing sequences, but more complete there
chesters_seqs <- chesters_seqs[!(names(chesters_seqs) %in% extra_meta$TipLabel)]
chesters_meta <- chesters_meta[!(chesters_meta$TipLabel %in% extra_meta$TipLabel),]

# Extract the coding part
chesters_seqs <- extract_coding(chesters_seqs, chesters_meta)
extra_seqs <- extract_coding(extra_seqs, extra_meta)
sikora_seqs <- extract_coding(sikora_seqs, sikora_meta)

# Write the resulting sequences
out_file <- "expanded_cecidomyiidae.fasta"
write.FASTA(chesters_seqs, out_file)
write.FASTA(extra_seqs, out_file, append=TRUE)
write.FASTA(sikora_seqs, out_file, append=TRUE)

