source("../code/fetch_fasta_fxn.R")

# Define function that fetches sequences into a file
# deleting the content of the file first
fetch_seqs <- function(D, seq_file) {
    cat("",file=seq_file)
    fetch_fasta(D, seq_file)
}

# Fetch outgroup sequences
seq_file <- "braconidae_outgroups_CO1.fasta"
T <- read.delim("braconidae_outgroups_taxonomy.tsv")
seqs <- fetch_seqs(T, seq_file)

# Fetch missing subfamily sequences
seq_file <- "missing_braconidae_subfamilies_CO1.fasta"
T <- read.delim("missing_braconidae_subfamilies_taxonomy.tsv")
seqs <- fetch_seqs(T, seq_file)

# Fetch missing genera with automatically generated metadata
seq_file <- "missing_braconidae_genera_CO1.fasta"
T <- read.delim("missing_braconidae_genera_taxonomy.tsv")
seqs <- fetch_seqs(T, seq_file)

# Fetch missing genera with manually retrieved metadata
seq_file <- "missing_braconidae_genera_checked_CO1.fasta"
T <- read.delim("missing_braconidae_genera_checked_taxonomy.tsv")
seqs <- fetch_seqs(T, seq_file)

