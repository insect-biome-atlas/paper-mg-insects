# Read in function that converts fasta to nexus
source("../code/seq_fxns.R")

# Read in functions that generate hard constraints and
# partial constraints from an input tree
source("../code/constraint_fxns.R")

# Generate nexus data file
fasta2nexus("expanded_cecidomyiidae_aligned_trimmed.fasta","mb_runs/cecidomyiidae.nex")

# Set name of constraint file
out_file <- "mb_runs/cecidomyiidae_constraints.nex"

# Print header to constraint file
cat("#NEXUS\n\nbegin mrbayes;\n", file=out_file)

# Generate partial constraints from sikora tree
# Their Fig. 1C, all supported clades > PP 95% in Bayesian
# analysis. Note that Pseudomonardia is ignored.
# It is included in the tree but not in the data
# table.
# ---------------------------------------------

sikora_meta <- read.delim("sikora_taxonomy.tsv")
sikora_meta$TipLabel <- sub(" ","_",sikora_meta$Species)
sikora_tips <- sikora_meta$TipLabel

# Build constraints successively
# Note that Pseudomonardia occurs in the tree
lestremiinae <- sikora_meta$TipLabel[sikora_meta$Subfamily=="Lestremiinae"]
add_partial_subset_constraint("lestremiinae",lestremiinae,sikora_tips,out_file)
catotrichinae <- sikora_meta$TipLabel[sikora_meta$Subfamily=="Catotrichinae"]
add_partial_subset_constraint("catotrichinae",catotrichinae,sikora_tips,out_file)
catochini <- sikora_meta$TipLabel[sikora_meta$Genus=="Catocha"]
add_partial_subset_constraint("catochini",catochini,sikora_tips,out_file)
aco_gro <- sikora_meta$TipLabel[sikora_meta$Genus=="Acoenonia" | sikora_meta$Genus=="Groveriella"]
add_partial_subset_constraint("aco_gro",aco_gro,sikora_tips,out_file)
campylomyzini <- sikora_meta$TipLabel[sikora_meta$Genus=="Campylomyza" | sikora_meta$Genus=="Neurolyga"]
add_partial_subset_constraint("campylomyzini",campylomyzini,sikora_tips,out_file)
peromyiini <- sikora_meta$TipLabel[sikora_meta$Genus=="Peromyia"]
add_partial_subset_constraint("peromyiini",peromyiini,sikora_tips,out_file)
bryomyiini <- sikora_meta$TipLabel[sikora_meta$Genus=="Skuhraviana" | sikora_meta$Genus=="Heterogenella" | sikora_meta$Genus=="Bryomyia"]
add_partial_subset_constraint("bryomyiini",bryomyiini,sikora_tips,out_file)
micromyini <- sikora_meta$TipLabel[sikora_meta$Genus=="Polyardis" | sikora_meta$Genus=="Monardia" | sikora_meta$Genus=="unclassified.Micromyiini"]
add_partial_subset_constraint("micromyini",micromyini,sikora_tips,out_file)
aprionini <- sikora_meta$TipLabel[sikora_meta$Genus=="Aprionus" | sikora_meta$Genus=="Tropaprionus"]
add_partial_subset_constraint("aprionini",aprionini,sikora_tips,out_file)
heteropezini <- sikora_meta$TipLabel[sikora_meta$Genus=="Nikandria" | sikora_meta$Genus=="Heteropeza" | sikora_meta$Genus=="Leptosyna"]
add_partial_subset_constraint("heteropezini",heteropezini,sikora_tips,out_file)
diallactiini1 <- sikora_meta$TipLabel[sikora_meta$Genus=="unclassified.Diallactiini" | sikora_meta$Genus=="Gynapteromyia"]
add_partial_subset_constraint("diallactiini1",diallactiini1,sikora_tips,out_file)
asynaptini <- sikora_meta$TipLabel[sikora_meta$Genus=="Camptomyia" | sikora_meta$Genus=="Svenartia"]
add_partial_subset_constraint("asynaptini",asynaptini,sikora_tips,out_file)
porricondylini <- sikora_meta$TipLabel[sikora_meta$Genus=="Zaitzeviola" | sikora_meta$Genus=="Porricondyla" | sikora_meta$Genus=="Coccopsilis" | sikora_meta$Genus=="Claspettomyia"]
add_partial_subset_constraint("porricondylini",porricondylini,sikora_tips,out_file)
dicerurini1 <- sikora_meta$TipLabel[sikora_meta$Genus=="Dicerura"]
add_partial_subset_constraint("dicerurini1",dicerurini1,sikora_tips,out_file)
dicerurini2 <- sikora_meta$TipLabel[sikora_meta$Genus=="Tetraneuromyia"]
add_partial_subset_constraint("dicerurini2",dicerurini2,sikora_tips,out_file)
stomatosematidi <- sikora_meta$TipLabel[sikora_meta$Supertribe=="Stomatosematidi"]
add_partial_subset_constraint("stomatosematidi",stomatosematidi,sikora_tips,out_file)
cecidomyiidi <- sikora_meta$TipLabel[sikora_meta$Supertribe=="Cecidomyiidi"]
add_partial_subset_constraint("cecidomyiidi",cecidomyiidi,sikora_tips,out_file)
lasiopteridi <- sikora_meta$TipLabel[sikora_meta$Supertribe=="Lasiopteridi"]
add_partial_subset_constraint("lasiopteridi",lasiopteridi,sikora_tips,out_file)
cec1 <- c(lasiopteridi,cecidomyiidi)
add_partial_subset_constraint("cec1",cec1,sikora_tips,out_file)
cecidomyiinae <- c(stomatosematidi,cec1)
add_partial_subset_constraint("cecidomyiinae",cecidomyiinae,sikora_tips,out_file)
por1 <- c(porricondylini,asynaptini)
add_partial_subset_constraint("por1",por1,sikora_tips,out_file)
por2 <- c(dicerurini1,por1)
add_partial_subset_constraint("por2",por2,sikora_tips,out_file)
pc1 <- c("Dirhiza_lateritia",cecidomyiinae)
add_partial_subset_constraint("pc1",pc1,sikora_tips,out_file)
pc2 <- c(dicerurini2,pc1)
add_partial_subset_constraint("pc2",pc2,sikora_tips,out_file)
pc3 <- c(por2,pc2)
add_partial_subset_constraint("pc3",pc3,sikora_tips,out_file)
wpc1 <- c(pc3,diallactiini1)
add_partial_subset_constraint("wpc1",wpc1,sikora_tips,out_file)
m1 <- c(micromyini,aprionini)
add_partial_subset_constraint("m1",m1,sikora_tips,out_file)
m2 <- c(peromyiini,bryomyiini)
add_partial_subset_constraint("m2",m2,sikora_tips,out_file)
m3 <- c(m1,m2)
add_partial_subset_constraint("m3",m3,sikora_tips,out_file)
m4 <- c(campylomyzini,m3)
add_partial_subset_constraint("m4",m4,sikora_tips,out_file)
m5 <- c(aco_gro,m4)
add_partial_subset_constraint("m5",m5,sikora_tips,out_file)
winnertziinae <- sikora_meta$TipLabel[sikora_meta$Subfamily=="Winnertziinae"]
porricondylinae <- sikora_meta$TipLabel[sikora_meta$Subfamily=="Porricondylinae"]
mwpc1 <- c(m5,winnertziinae,porricondylinae,cecidomyiinae)
add_partial_subset_constraint("mwpc1",mwpc1,sikora_tips,out_file)
mwpc <- c(mwpc1,catochini)
add_partial_subset_constraint("mwpc",mwpc,sikora_tips,out_file)
cmwpc <- c(mwpc,catotrichinae)
add_partial_subset_constraint("cmwpc",cmwpc,sikora_tips,out_file)

# We ignore relationships in Chesters tree

# Generate hard constraints for families
# --------------------------------------

# Read in metadata
# We do not need the detailed classification of the
# chesters data for the current task
chesters_taxa <- read.delim("chesters_cecidomyiidae_taxonomy.tsv")
extra_taxa <- read.delim("missing_cecidomyiidae_taxonomy.tsv")
extra_taxa$TipLabel <- sub(" ","_",extra_taxa$Species)
sikora_taxa <- read.delim("sikora_taxonomy.tsv")
sikora_taxa$TipLabel <- sub(" ","_",sikora_taxa$Species)

# Assemble all tip labels and higher classification info
taxa <- data.frame(TipLabel=c(chesters_taxa$TipLabel, extra_taxa$TipLabel, sikora_taxa$TipLabel))
taxa$Family <- c(chesters_taxa$Family, extra_taxa$Family, sikora_taxa$Family)

# We append the constraints to the same file
# as a separate mrbayes block

# Print block header to output file
cat("\nbegin mrbayes;\n", file=out_file, append=TRUE)

# Output constraint partitions for families
for (fam in unique(taxa$Family)) {
    if (fam=="")
        next
    ingroup <- taxa$TipLabel[taxa$Family==fam]
    add_hard_constraint(fam, ingroup, out_file)
}

# Print tail to output file
cat("end;\n", file=out_file, append=TRUE)

