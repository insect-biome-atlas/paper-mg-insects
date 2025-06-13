# Script for updating the taxonomy for the Chesters Ichneumonidae sequences

# Read in needed functions
source("../code/taxonomy_fxns.R")

# Read in original taxonomy
T <- read.delim("chesters_ichneumonidae_taxonomy.tsv")

# Update the taxonomy
T$TipLabel <- sub(" ","_",T$Species)

# Get subfamily, tribe and subtribe info from NCBI (supertribes are not used in Ichneumonidae)
ncbi <- read.delim("../ncbi/ncbi_arthropod_genus_classification_detailed.tsv")
T$Subfamily <- ncbi$Subfamily[match(T$Genus,ncbi$Genus)]
T$Tribe <- ncbi$Tribe[match(T$Genus,ncbi$Genus)]
T$Subtribe <- ncbi$Subtribe[match(T$Genus,ncbi$Genus)]

# Corrections in the NCBI classification (see mail from Gavin Broad 20204-03-29)

# Ecphysis now in Claseinae.
T$Subfamily[T$Genus=="Ecphysis"] <- "Claseinae"
T$Tribe[T$Genus=="Ecphysis"] <- ""

# Species in subfamily Phygadeuontinae were also losted as belonging to tribe Phygadeuontini, which I've corrected to just subfamily Phygadeuontinae.
T$Tribe[T$Subfamily=="Phygadeuontinae"] <- ""

# Neotheronia and Theronia are now in the tribe Theroniini.
T$Tribe[T$Genus %in% c("Theronia","Neotheronia")] <- "Theroniini"

# Tribe Hemigasterini should now be Aptesini
T$Tribe[T$Tribe=="Hemigasterini"] <- "Aptesini"

# Tribe Heresiarchini is now a synonym of tribe Ichneumonini.
T$Tribe[T$Tribe=="Heresiarchini"] <- "Ichneumonini"

# Proclitus is in Orthocentrinae (Helictinae is a synonym).
T$Subfamily[T$Genus=="Proclitus"] <- "Orthocentrinae"

# Astrenis is now in Tersilochinae (Phrudinae synonymised)
T$Subfamily[T$Genus=="Atrenis"] <- "Tersilochinae"

# Seleucus is now in Ctenopelmatinae in its own tribe
T$Subfamily[T$Genus=="Seleucus"] <- "Ctenopelmatinae"
T$Tribe[T$Genus=="Seleucus"] <- "Seleucini"

# Echthrus is in Cryptini
T$Tribe[T$Genus=="Echthrus"] <- "Cryptini"

# Erythrodolius is in Sisyrostolinae
T$Subfamily[T$Genus=="Erythrodolius"] <- "Sisyrostolinae"

# Hyperacmus is in Cylloceriinae
T$Subfamily[T$Genus=="Hyperacmus"] <- "Cylloceriinae"

# Exenterini is a synonym of Tryphonini
T$Tribe[T$Tribe=="Exenterini"] <- "Tryphonini"

# Write the table with updated taxonomy
write.table(T,"chesters_ichneumonidae_taxonomy_updated.tsv",sep="\t",row.names=FALSE)

