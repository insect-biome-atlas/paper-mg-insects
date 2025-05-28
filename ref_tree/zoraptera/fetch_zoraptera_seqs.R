# Dependent on ape
library(ape)

# Read in needed functions
source("../code/fetch_fasta_fxn.R")
source("../code/seq_fxns.R")

# Set data path
meta <- read.delim("zoraptera_taxonomy.tsv")
cat("",file="zoraptera_CO1.fasta")
fetch_fasta(meta,"zoraptera_CO1.fasta")

# Read in sequences
seqs <- read.FASTA("zoraptera_CO1.fasta")

# Relabel missing sequences from NCBI strings to tip labels
seqs <- rename_seqs(seqs,meta)

# Extract the coding part
seqs <- extract_coding(seqs,meta)

# Write the resulting sequences
out_file <- "zoraptera.fasta"
write.FASTA(seqs, out_file)
