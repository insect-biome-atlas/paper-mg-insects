library(ape)

source("../code/seq_fxns.R")

# Read in outgroups
seqs <- read.FASTA("braconidae_outgroups_CO1.fasta")
seqs_meta <- read.delim("braconidae_outgroups_taxonomy.tsv")
seqs_meta <- add_tiplabel(seqs_meta)
seqs <- rename_seqs(seqs, seqs_meta)
seqs <- extract_coding(seqs, seqs_meta)

merged_tax <- seqs_meta[,c("TipLabel","Family","Subfamily","Genus","Species")]
merged_seqs <- seqs

# Add Chesters data
# The taxonomy already contains TipLabel, and the seqs are labeled
# with TipLabel and shortened to the coding part
seqs <- read.FASTA("chesters_braconidae.fasta")
seqs_meta <- read.delim("chesters_braconidae_taxonomy.tsv")
temp <- read.delim("chesters_braconidae_taxonomy_updated.csv",sep=";")
seqs_meta$Subfamily <- temp$Subfamily_updated[match(seqs_meta$TipLabel, temp$TipLabel)]

merged_tax <- rbind(merged_tax, seqs_meta[,c("TipLabel","Family","Subfamily","Genus","Species")])
merged_seqs <- c(merged_seqs,seqs)

# Add missing subfamilies
seqs <- read.FASTA("missing_braconidae_subfamilies_CO1.fasta")
seqs_meta <- read.delim("missing_braconidae_subfamilies_taxonomy.tsv")
seqs_meta <- add_tiplabel(seqs_meta)
seqs <- rename_seqs(seqs, seqs_meta)
seqs <- extract_coding(seqs, seqs_meta)

merged_tax <- rbind(merged_tax, seqs_meta[,c("TipLabel","Family","Subfamily","Genus","Species")])
merged_seqs <- c(merged_seqs,seqs)

# Add missing automatically added genera
seqs <- read.FASTA("missing_braconidae_genera_CO1.fasta")
seqs_meta <- read.delim("missing_braconidae_genera_taxonomy.tsv")
seqs_meta <- add_tiplabel(seqs_meta)
seqs_meta$Family <- "Braconidae"
seqs <- rename_seqs(seqs, seqs_meta)
seqs <- extract_coding(seqs, seqs_meta)

merged_tax <- rbind(merged_tax, seqs_meta[,c("TipLabel","Family","Subfamily","Genus","Species")])
merged_seqs <- c(merged_seqs,seqs)

# Add missing manually added genera
seqs <- read.FASTA("missing_braconidae_genera_checked_CO1.fasta")
seqs_meta <- read.delim("missing_braconidae_genera_checked_taxonomy.tsv")
seqs_meta <- add_tiplabel(seqs_meta)
seqs_meta$Family <- "Braconidae"
seqs <- rename_seqs(seqs, seqs_meta)
seqs <- extract_coding(seqs, seqs_meta)

merged_tax <- rbind(merged_tax, seqs_meta[,c("TipLabel","Family","Subfamily","Genus","Species")])
merged_seqs <- c(merged_seqs, seqs)

# Expand merged_tax
merged_tax$Kingdom <- "Animalia"
merged_tax$Phylum <- "Arthropoda"
merged_tax$Class <- "Insecta"
merged_tax$Order <- "Hymenoptera"
merged_tax$Clade <- merged_tax$Subfamily
merged_tax <- merged_tax[,c("TipLabel","Kingdom","Phylum","Class","Order","Family","Subfamily","Clade","Genus","Species")]

# Write merged data
write.FASTA(merged_seqs, file="expanded_braconidae.fasta")
write.table(merged_tax, file="expanded_braconidae_taxonomy.tsv", row.names=FALSE, sep="\t")

