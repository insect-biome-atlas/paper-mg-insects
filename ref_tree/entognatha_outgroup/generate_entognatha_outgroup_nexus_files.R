# Read in function that converts fasta to nexus
source("../code/seq_fxns.R")

# Read in functions that generate hard constraints and
# partial constraints from an input tree
source("../code/constraint_fxns.R")

# Generate nexus data file
fasta2nexus("entognatha_outgroup_aligned_trimmed.fasta","mb_runs/entognatha_outgroup.nex")

# Generate collembola constraints
gen_mb_con_file("collembola_tree.nwk", "mb_runs/collembola_constraints.nex")

# Generate higher constraints
out_file <- "mb_runs/higher_constraints.nex"

# Read in metadata
root_taxa <- read.delim("root_taxonomy.tsv")
collembola_taxa <- read.delim("collembola_taxonomy.tsv")
entognatha_taxa <- read.delim("missing_entognatha_taxonomy.tsv")

# Print header to output file
cat("#NEXUS\n\nbegin mrbayes;\n",file=out_file)

# Assemble all tip labels and higher classification info
taxa <- data.frame(TipLabel=c(root_taxa$TipLabel, collembola_taxa$TipLabel, entognatha_taxa$TipLabel))
n <- nrow(taxa) - nrow(root_taxa)
taxa$Kingdom <- c(root_taxa$Kingdom, rep("Animalia",n))
taxa$Phylum <- c(root_taxa$Phylum, rep("Arthropoda",n))
taxa$Class <- c(root_taxa$Class, collembola_taxa$Class, entognatha_taxa$Class)
taxa$Order <- c(root_taxa$Order, collembola_taxa$Order, entognatha_taxa$Order)

# Output informative constraint partitions for higher taxa
for (kingdom in unique(taxa$Kingdom)) {
    ingroup <- taxa$TipLabel[taxa$Kingdom==kingdom]
    add_hard_constraint(kingdom, ingroup, out_file)
}

for (phylum in unique(taxa$Phylum)) {
    ingroup <- taxa$TipLabel[taxa$Phylum==phylum]
    add_hard_constraint(phylum, ingroup, out_file)
}


for (class in unique(taxa$Class)) {
    ingroup <- taxa$TipLabel[taxa$Class==class]
    add_hard_constraint(class, ingroup, out_file)
}

for (order in unique(taxa$Order)) {
    ingroup <- taxa$TipLabel[taxa$Order==order]
    add_hard_constraint(order, ingroup, out_file)
}

# Add special constraints
# -----------------------

# Opisthokonta
ingroup <- taxa$TipLabel[taxa$Kingdom!="Diphoda"]
add_hard_constraint("Opisthokonta", ingroup, out_file)

# Protostomia (Deuterostomia only includes Chordata, so not needed)
ingroup <- taxa$TipLabel[taxa$Kingdom=="Animalia" & !(taxa$Phylum=="Chordata")]
add_hard_constraint("Protostomia", ingroup, out_file)

# Lophotrochozoa
ingroup <- taxa$TipLabel[taxa$Phylum=="Mollusca" | taxa$Phylum=="Annelida"]
add_hard_constraint("Lophotrochozoa", ingroup, out_file)

# Ecdysozoa
ingroup <- taxa$TipLabel[taxa$Phylum=="Arthropoda" | taxa$Phylum=="Tardigrada" | taxa$Phylum=="Nematoda"]
add_hard_constraint("Ecdysozoa", ingroup, out_file)

# Panarthropoda
ingroup <- taxa$TipLabel[taxa$Phylum=="Arthropoda" | taxa$Phylum=="Tardigrada"]
add_hard_constraint("Panarthropoda", ingroup, out_file)

# Mandibulata
ingroup <- taxa$TipLabel[taxa$Phylum=="Arthropoda" & !(taxa$Class=="Arachnida")]
add_hard_constraint("Mandibulata", ingroup, out_file)

# Pancrustacea
ingroup <- taxa$TipLabel[taxa$Phylum=="Arthropoda" & !(taxa$Class=="Arachnida") & !(taxa$Class=="Diplopoda")]
add_hard_constraint("Pancrustacea", ingroup, out_file)

# Hexapoda
ingroup <- taxa$TipLabel[taxa$Class=="Insecta" | taxa$Class=="Collembola" | taxa$Class=="Diplura" | taxa$Class=="Protura"]
add_hard_constraint("Hexapoda", ingroup, out_file)

# Dicondylia
ingroup <- taxa$TipLabel[taxa$Order=="Zygentoma" | taxa$Order=="Ephemeroptera" | taxa$Order=="Diptera"]
add_hard_constraint("Dicondylia", ingroup, out_file)

# Pterygota
ingroup <- taxa$TipLabel[taxa$Order=="Ephemeroptera" | taxa$Order=="Diptera"]
add_hard_constraint("Pterygota", ingroup, out_file)

# Print tail to output file
cat ("end;\n", file=out_file, append=TRUE)

