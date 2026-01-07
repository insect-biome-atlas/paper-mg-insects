library(ape)
source("../code/seq_fxns.R")


# 1. Merge taxonomy information
# =============================

# Read in taxonomic info
data_path <- ""
T1 <- read.delim(paste0(data_path,"collembola_taxonomy.tsv"))
T2 <- read.delim(paste0(data_path,"missing_entognatha_taxonomy.tsv"))
T3 <- read.delim(paste0(data_path,"root_taxonomy.tsv"))

T1$Kingdom <- "Animalia"
T1$Phylum <- "Arthropoda"
T1$TranslationTable <- 5

T2$Kingdom <- "Animalia"
T2$Phylum <- "Arthropoda"
T2$TranslationTable <- 5

# Remove Sinentomon from the Collembola files, as it is given in more detail in the missing Entognatha files.
T1 <- T1[T1$Genus!="Sinentomon",]

# Synonymize Paronellidae with Entomobryidae (see comment in the traits source file)
T1$Family[T1$Family=="Paronellidae"] <- "Entomobryidae"

cols <- c("TipLabel","Kingdom","Phylum","Class","Order","Family","Genus","Species","TranslationTable")

# Merge taxon info files
T <- rbind(T1[,cols],T2[,cols],T3[,cols])

# Add info on Clade, matching life history data file

# For now, we set Clade equal to family; it might change in future versions
T$Clade <- T$Family

# Reorganize taxonomic info
cols <- c("TipLabel","Kingdom","Phylum","Class","Order","Family","Clade","Genus","Species","TranslationTable")
T <- T[,cols]

# Write merged taxonomy file
write.table(T,paste0(data_path,"entognatha_outgroup_taxonomy.tsv"),sep="\t",row.names=FALSE)


# 2. Merge sequences
# ==================


# Generate root sequences
# -----------------------

seqs <- read.FASTA(paste0(data_path,"root_mt_genomes.fasta"))
meta <- read.delim(paste0(data_path,"root_taxonomy.tsv"))

# Relabel sequences with species names
root_seqs <- rename_seqs(seqs, meta)

# Extract the coding part
root_seqs <- extract_coding(root_seqs, meta)


# Generate extra entognatha sequences
# -----------------------------------

seqs <- read.FASTA(paste0(data_path,"missing_entognatha_CO1.fasta"))
meta <- read.delim(paste0(data_path,"missing_entognatha_taxonomy.tsv"))

# Relabel sequences with species names
ento_seqs <- rename_seqs(seqs, meta)

# Extract the coding part
ento_seqs <- extract_coding(ento_seqs, meta)


# Read in Collembola sequences
# ----------------------------

seqs <- read.FASTA(paste0(data_path,"collembola_CO1.fasta"))
meta <- read.delim(paste0(data_path,"collembola_taxonomy.tsv"))

# Sequences should already be named correctly with tip labels

# Remove the Sinentomon sequence (also in missing Entognatha sequence file)
seqs <- seqs[names(seqs)!="Sinentomon_erythranum"]

# Extract the coding part
seqs <- extract_coding(seqs, meta)


# Write the resulting sequences
# -----------------------------

out_file <- paste0(data_path,"entognatha_outgroup.fasta")
write.FASTA(seqs, out_file)
write.FASTA(ento_seqs, out_file, append=TRUE)
write.FASTA(root_seqs, out_file, append=TRUE)

