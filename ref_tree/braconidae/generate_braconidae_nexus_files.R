# Read in function that converts fasta to nexus
source("../code/seq_fxns.R")

# Read in functions that generate hard constraints and
# partial constraints from an input tree
source("../code/constraint_fxns.R")

# Generate nexus data file
fasta2nexus("expanded_braconidae_aligned_trimmed.fasta","mb_runs/braconidae.nex")

# Set name of constraint file
constraints_file <- "mb_runs/braconidae_constraints.nex"

# Print block header to constraint file
cat("#NEXUS\n\nbegin mrbayes;\n", file=constraints_file)

# Read in taxonomy for expanded Braconidae data
D <- read.delim("expanded_braconidae_taxonomy.tsv")

# Generate taxon sets for families 
fams <- unique(D$Family)
for (fam in fams) {
    taxon_set <- D$TipLabel[D$Family == fam]
    add_taxset(fam, taxon_set, constraints_file)
}

# Generate taxon sets for subfamilies 
subfams <- unique(D$Subfamily[D$Family=="Braconidae" & !grepl("?", D$Subfamily, fixed=TRUE)])
for (subfam in subfams) {
    taxon_set <- D$TipLabel[D$Subfamily == subfam]
    add_taxset(subfam, taxon_set, constraints_file)
}

# Generate hard family constraints
for (fam in fams) {
    add_hard_constraint(fam, fam, constraints_file, taxset=TRUE)
}

# Generate partial constraints for subfamilies
# We do not place Ontsira because it is uncertain if the included sequence matches
# the species analyzed by Jasso-Martinez et al, and our preliminary analyses did
# not clarify its placement
placed_taxa <- D$TipLabel[D$Family=="Braconidae" & !grepl("Doryctinae_Ontsira",D$Subfamily, fixed=TRUE)]
for (subfam in subfams) {
    ingroup <- D$TipLabel[D$Subfamily==subfam]
    outgroup <- placed_taxa[!(placed_taxa %in% ingroup)]
    add_soft_constraint(subfam, ingroup, outgroup, constraints_file)
}

# Manual coding of higher clades in UCE tree
clades <- list()
clades[[1]] <- c("Alysiinae", "Opiinae")
clades[[2]] <- c(clades[[1]], "Exothecinae")
clades[[3]] <- c("Braconinae", "Telengaiinae")
clades[[4]] <- c(clades[[2]], clades[[3]])
clades[[5]] <- c("Rogadinae", "Hormiinae")
clades[[6]] <- c(clades[[5]], "Rhysipolinae")
clades[[7]] <- c(clades[[4]], clades[[6]], "Mesostoinae_Xenosternum_Avga", "Doryctinae_s_str")
clades[[8]] <- c(clades[[7]], "Doryctinae_South_America")
clades[[9]] <- c(clades[[8]], "Rhyssalinae")                # cyclostomes_s_str
clades[[10]] <- c("Mesostoinae_s_str", "Aphidiinae")        # aphidioid complex   
clades[[11]] <- c(clades[[10]], "Masoninae")
clades[[12]] <- c(clades[[9]], clades[[11]])                # cyclostomes_s_lat

clades[[13]] <- c("Miracinae", "Cardiochilinae")
clades[[14]] <- c(clades[[13]], "Microgastrinae")
clades[[15]] <- c(clades[[14]], "Khoikhoiinae")
clades[[16]] <- c(clades[[15]], "Mendesellinae")
clades[[17]] <- c(clades[[16]], "Cheloninae")               # microgastroid complex
clades[[18]] <- c(clades[[17]], "Proteropinae")             # proteropinae + microgastroid complex
clades[[19]] <- c(clades[[18]], "Proteropinae")
clades[[20]] <- c("Agathidinae", "Sigalphinae")
clades[[21]] <- c(clades[[20]], "Ichneutinae")              # sigalphoid complex
clades[[22]] <- c(clades[[21]], clades[[18]])
clades[[23]] <- c("Cenocoeliinae", "Euphorinae")            # euphoroid complex
clades[[24]] <- c(clades[[23]], clades[[22]])

clades[[25]] <- c("Macrocentrinae", "Charmontinae")
clades[[26]] <- c("Homolobinae", "Microtypinae")
clades[[27]] <- c(clades[[26]], "Orgilinae")
clades[[28]] <- c(clades[[25]], clades[[27]]) 
clades[[29]] <- c(clades[[28]], "Helconinae")
clades[[30]] <- c(clades[[29]], "Brachistinae_s_str")
clades[[31]] <- c(clades[[30]], "Acampsohelconinae")        # helconoid complex
clades[[32]] <- c(clades[[24]], clades[[31]])               # non_cyclostomes_s_str

clades[[33]] <- c(clades[[32]], "Trachypetinae")
clades[[34]] <- c(clades[[33]], "Meteorideinae")            # non_cyclostomes_s_lat


# Generate partial constraints for well supported (>50%) higher clades in UCE tree
for (i in 1:length(clades)) {
    ingroup <- clades[[i]]
    outgroup <- subfams[!(subfams %in% ingroup) & subfams!="Doryctinae_Ontsira"]
    add_soft_constraint(paste0("clade",i), ingroup, outgroup, constraints_file)
}

# Print tail to output file
cat("end;\n", file=constraints_file, append=TRUE)

