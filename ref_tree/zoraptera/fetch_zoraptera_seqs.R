# Dependent on ape
library(ape)

# Read in needed functions
source("../code/fetch_fasta_fxn.R")
source("../code/seq_fxns.R")

# Set data path
meta <- read.delim("missing_zoraptera_ncbi.tsv")
cat("",file="missing_zoraptera_CO1.fasta")
fetch_fasta(meta,"missing_zoraptera_CO1.fasta")

# Read in the fetched sequences
seqs <- read.FASTA("missing_zoraptera_CO1.fasta")

# Relabel missing sequences from NCBI strings to tip labels
seqs <- rename_seqs(seqs,meta)

# Extract the coding part
seqs <- extract_coding(seqs,meta)

# Write the resulting sequences
out_file <- "zoraptera.fasta"
write.FASTA(seqs, out_file)

# Update taxonomy file
meta$Kingdom <- "Animalia"
meta$Phylum <- "Arthropoda"
meta$Class <- "Insecta"
meta$Clade <- meta$Family
meta$TipLabel <- sub(" ","_", meta$Species)
meta <- meta[,c("GenBank","Order","Family","Genus","Species","Start","Stop","TipLabel","Kingdom","Phylum","Class","Order","Clade")]

# Write the final taxonomy file
write.table(meta,"zoraptera_taxonomy.tsv")
