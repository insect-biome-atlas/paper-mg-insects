# Read in function that converts fasta to nexus
source("../code/seq_fxns.R")

# Read in constraint functions
source("../code/constraint_fxns.R")

# Generate nexus data file
fasta2nexus("emb_pha_not_aligned_trimmed.fasta","mb_runs/emb_pha_not.nex")

# Set name of constraint file
out_file <- "mb_runs/emb_pha_not_constraints.nex"

# Generate hard constraints for orders
# --------------------------------------

# Read in metadata
taxa <- read.delim("emb_pha_not_taxonomy.tsv")

# Print block header to output file
cat("\nbegin mrbayes;\n", file=out_file, append=TRUE)

# Output constraint partitions for orders
for (order in unique(taxa$Order)) {
    ingroup <- taxa$TipLabel[taxa$Order==order]
    add_hard_constraint(order, ingroup, out_file)
}

# Add constraints for higher-level clades in Misof et al (2014)
emb_pha <- taxa$TipLabel[taxa$Order %in% c("Embioptera","Phasmatodea")]
add_hard_constraint("emb_pha", emb_pha, out_file)
notoptera <- taxa$TipLabel[taxa$Order %in% c("Grylloblattodea","Mantophasmatodea")]
add_hard_constraint("notoptera", notoptera, out_file)

# Print tail to output file
cat("end;\n", file=out_file, append=TRUE)

