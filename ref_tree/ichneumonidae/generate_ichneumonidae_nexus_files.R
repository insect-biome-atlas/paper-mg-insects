# Read in function that converts fasta to nexus
source("../code/seq_fxns.R")

# Read in functions that generate hard constraints and
# partial constraints from an input tree
source("../code/constraint_fxns.R")

# Generate nexus data file
fasta2nexus("expanded_ichneumonidae_aligned_trimmed.fasta","mb_runs/ichneumonidae.nex")

# Set name of constraint file
out_file <- "mb_runs/ichneumonidae_constraints.nex"

# Generate partial constraints from chesters tree
gen_mb_con_file("chesters_ichneumonidae.nwk", out_file)

# Generate hard constraints for families
# --------------------------------------

# Read in metadata
taxa <- read.delim("expanded_ichneumonidae_taxonomy.tsv")

# We append the constraints to the same file
# as a separate mrbayes block

# Print block header to output file
cat("\nbegin mrbayes;\n", file=out_file, append=TRUE)

# Output constraint partitions for families
for (fam in unique(taxa$Family)) {
    ingroup <- taxa$TipLabel[taxa$Family==fam]
    add_hard_constraint(fam, ingroup, out_file)
}

# Print tail to output file
cat("end;\n", file=out_file, append=TRUE)

