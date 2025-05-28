# Dependent on ape
library(ape)

# Read in needed functions
source("../code/fetch_fasta_fxn.R")
source("../code/seq_fxns.R")

# Fetch sequence data
recs <- read.delim("emb_pha_not_ncbi_recs.tsv")
cat("",file="emb_pha_not_CO1.fasta")    # start from scratch
fetch_fasta(recs,"emb_pha_not_CO1.fasta")   # collect fasta recs

# Finalize taxonomy data
meta <- recs
meta$TipLabel <- gsub(" ", "_", meta$Species)
meta$Kingdom <- "Animalia"
meta$Phylum <- "Arthropoda"
meta$Class <- "Insecta"
meta$Clade <- meta$Family
write.table(meta, "emb_pha_not_taxonomy.tsv", sep="\t", row.names=FALSE)

# Read in sequences
seqs <- read.FASTA("emb_pha_not_CO1.fasta")

# Relabel missing sequences from NCBI strings to tip labels
seqs <- rename_seqs(seqs,meta)

# Extract the coding part
seqs <- extract_coding(seqs,meta)

# Write the resulting sequences
out_file <- "emb_pha_not.fasta"
write.FASTA(seqs, out_file)

