# Read in function that converts fasta to nexus
source("../code/seq_fxns.R")

# Read in constraint functions
source("../code/constraint_fxns.R")

# Generate nexus data file
fasta2nexus("zoraptera_aligned_trimmed.fasta","mb_runs/zoraptera.nex")

# Set name of constraint file
out_file <- "mb_runs/zoraptera_constraints.nex"

# Generate hard constraints for orders
# --------------------------------------

# Read in metadata
taxa <- read.delim("zoraptera_taxonomy.tsv")

# Print block header to output file
cat("NEXUS\n\nbegin mrbayes;\n", file=out_file)

# Output constraint partitions for orders
for (order in unique(taxa$Order)) {
    ingroup <- taxa$TipLabel[taxa$Order==order]
    add_hard_constraint(order, ingroup, out_file)
}

# Print tail to output file
cat("end;\n", file=out_file, append=TRUE)

