# Script for extracting taxonomy information from the tree in the
# chalcidoid "bush of life" analysis (Cruaud et al 2024, Cladistics)

# Libraries needed
library(ape)

# Read in the combined tree from Cruaud et al
tree <- read.tree("IQ_COMBINED.tre")

# Create a data frame for the taxonomy information
# contained in the tip labels
n <- length(tree$tip.label)
D <- data.frame(list(TipLabel=character(n),
                     Familycode=character(n),
                     Subfamily=character(n),
                     Tribe=character(n),
                     Genus=character(n),
                     Mergecode=character(n),
                     Binarycode=character(n),
                     Comment=character(n)
                     ))

# Split the tip labels into components
a <- strsplit(tree$tip.label,split="_")

# Cycle over the parsed strings and extract
# the relevant information
for (i in 1:length(a)) {
    names <- a[[i]]
    D$TipLabel[i] <- tree$tip.label[i]
    D$Familycode[i] <- names[1]
    D$Subfamily[i] <- names[2]
    last_elem <- names[length(names)]
    if (grepl("only",last_elem)) {
        D$Comment[i] <- last_elem
        names <- names[1:(length(names)-1)]
    }
    len <- nchar(names[3])
    if (substr(names[3],start=len-2,stop=len)=="ini")
        D$Tribe[i] <- names[3]
    len <- length(names)
    D$Genus[i] <- names[len-2]
    D$Mergecode[i] <- names[len-1]
    D$Binarycode[i] <- names[len]
    if (D$Subfamily[i]==D$Genus[i])
        D$Subfamily[i] <- ""
}

D$Sample_code_combined <- paste(D$Mergecode,D$Binarycode,sep="_")

E <- read.delim("cruaud_et_al_taxon_sampling.csv",sep=";")
E <- E[grepl("MERG",E$Sample_code_combined),]

D$Species_epithet <- "sp"
for (i in 1:nrow(D)) {
    idx <- match(D$Sample_code_combined[i], E$Sample_code_combined)
    if (E$AHE_species[idx]==E$UCE_species[idx] && !E$AHE_species[idx]=="sp")
        D$Species_epithet[i] <- E$AHE_species[idx]
}
D$Species <- paste(D$Genus, D$Species_epithet)

# Get additional information from the taxon sample table
idx <- match(D$Sample_code_combined,E$Sample_code_combined)
D$OldFamily <- E$Family[idx]
D$OldSubfamily <- E$Subfamily[idx]
D$OldTribe <- E$Tribe[idx]
D$TentativeNewFamily <- E$TentativeNewFamily[idx]
D$TentativeNewSubfamily <- E$TentativeNewSubFamily[idx]
D$TentativeNewTribe <- E$TentativeNewTribe[idx]
D$IN.OUT <- E$IN.OUT[idx]
D$AHE_genus <- E$AHE_genus[idx]
D$AHE_species <- E$AHE_species[idx]
D$UCE_genus <- E$UCE_genus[idx]
D$UCE_species <- E$UCE_species[idx]

# Write resulting table
write.table(D,"cruaud_taxonomy.tsv",sep="\t",row.names=FALSE)

