# Script for correcting and extending raw family traits for MG coded by
# Greg Lamarre (see binary file 'MG_Hexapoda_Family_traits.xlsx')

# Read in help functions
source("trait_fxns.R")


# Read in original data from csv file exported from Excel
# =======================================================

cat("Reading in original data\n")

D <- read.delim("mg_hexapoda_family_traits.csv",sep=";")
D <- D[,c(1:3,5:8)]
column_names <- c("Class","Order","Family","Niche","Habitat","Length_min","Length_max")
colnames(D) <- column_names

# Correct the data in various ways
#=================================

cat("Correcting original data\n")

# Note that Cyphoderidae is not in the NCBI classification but we have
# trait information for this family from Collembola experts, who use
# this family in their classification, so it is not a problem
# The family will be in the expanded Chesters tree.

# Correct niche and habitat information for taxa incorrectly coded by Lamarre
# (typing errors, coding more detailed than specified)
D$Niche[D$Niche=="Parasite "] <- "Parasite"
D$Niche[D$Niche=="Phytophagous"] <- "Phytophage"
D$Niche[D$Niche=="Predators"] <- "Predator"
D$Niche[D$Niche=="Lichen - decaying plant mat"] <- "Saprophage"
D$Niche[D$Niche=="Lichen"] <- "Saprophage"
D$Habitat[D$Habitat=="wood"] <- "Wood"
D$Habitat[D$Habitat=="Homoptera"] <- "Plants"

# Update information for entries marked in red by Lamarre
# Abbreviations:
#  - IoA: insects of Australia
#  - GBIF: GBIF backbone classification
#  - CoL: Insects of Australia 
#  - R20: Ronquist et al (2020) Completing Linnaeus's inventory... (PLoS ONE).

# Blattellidae
# ============
# Subfamily of Ectobiidae in CoL (most genera and spp, listed as separate family also, by mistake?)
# In GBIF listed as synonym of Ectobiidae.
# Separate family status is consistent with phylogeny on Wikipedia. Trait data OK.

# Nitidulidae
# ===========
# Coded for majority of genera/species according to IoA (p649)
idx <- match("Nitidulidae", D$Family)
D$Niche[idx] <- "Saprophage"
D$Habitat[idx] <- "Plants"

# Synneuridae
# ===========
# Should be corrected in the Chesters tree to Canthyloscelidae
idx <- match("Synneuridae", D$Family)
D$Family[idx] <- "Canthyloscelidae"
D$Niche[idx] <- "Saprophage"
D$Habitat[idx] <- "Wood"

# Xylophagaidae
# =============
# Should be corrected to Xylophagidae (see note by Lamarre)
# Trait data OK.
D$Family[D$Family=="Xylophagaidae"] <- "Xylophagidae"

# Hormaphidinae
# =============
# Treated as a separate family by NCBI but as a subfamily of Aphididae
# in GBIF and CoL. For some reason the subfamily name was used in
# Lamarre´s file. The trait data are OK.
idx <- match("Hormaphidinae", D$Family)
D$Family[idx] <- "Hormaphididae"

# Miridae
# =======
# Lamarre expressed some concern but the trait data appear to be OK.
# Some species are clearly predatory but not the majority apparently.

# Formicidae
# ==========
# This family is tricky. Most species appear to be omnivores, with either saprophage
# or predator as reasonable codings. In the choice between these two, saprophage
# appears the most commonly used niche, as suggested by R20 for Sweden.
# Soil appears to be the appropriate coding for habitat, as suggested by R20 for Sweden.
idx <- match("Formicidae", D$Family)
D$Niche[idx] <- "Saprophage"
D$Habitat[idx] <- "Soil"

# Torymidae
# =========
# The majority of species appear to be phytophage-parasitoids
# and not phytophages
idx <- match("Torymidae", D$Family)
D$Niche[idx] <- "Phytophage-parasitoid"

# Autostichidae
# =============
# This group (the taxa present in Australia) are treated as a subfamily
# of Oecophoridae. The trait data appear to be OK, with the correction
# introduced above and repeated here just in case.
idx <- match("Autostichidae", D$Family)
D$Niche[idx] <- "Saprophage"
D$Habitat[idx] <- "Plants"

# Cossidae
# ========
# Saprophage and phytophage tree borers are difficult to separate and the
# two groups grade into each other. The long development time of cossid
# larvae suggests that the classification as phytophage is correct. This
# is also consistent with the information in IoA. Thus, the trait data
# appear to be correct (and is also constent with that of R20).

# Epipyropidae
# ============
# Coded by Lamarre as "Parasite" and "Homoptera", but this should probably
# be "Parasite" and "Insects", or possibly "Phytophage-parasitoid" and "Plants".
idx <- match("Epipyropidae", D$Family)
D$Niche[idx] <- "Parasite"
D$Habitat[idx] <- "Insects"

# Hepialidae
# ==========
# Coded by Lamarre as "Saprophage" and "Wood", with an indication that their
# Panama database coded them as "Phytopage" and "Wood". However, although there
# are some species with this life history, IoA indicates that most
# Australian species should rather be coded as "Phytophage" and "Soil",
# which is also the coding used by R20.
idx <- match("Hepialidae", D$Family)
D$Niche[idx] <- "Phytophage"
D$Habitat[idx] <- "Soil"

# Meessiidae
# ==========
# Coded by Lamarre as "Lichen" and "Plants", which should clearly be
# "Saprophage" and "Plants".
idx <- match("Meessiidae", D$Family)
D$Niche[idx] <- "Saprophage"
D$Habitat[idx] <- "Plants"

# Oecophoridae
# ============
# The niche coding should be "Saprophage". Lamarre suggested "Wood"
# for habitat, but it appears that "Soil" (among leaf litter) is
# more common that "Plants" or "Wood" (bark).
idx <- match("Oecophoridae", D$Family)
D$Niche[idx] <- "Saprophage"
D$Habitat[idx] <- "Soil"

# Psychidae
# =========
# Lamarre suggests "Predator?" as niche coding but most species
# seem to be phytophagous. This contrasts with Sweden, were many
# species seem to be saprophagous on occur on the ground.
idx <- match("Psychidae", D$Family)
D$Niche[idx] <- "Phytophage"
D$Habitat[idx] <- "Plants"

# Xyloryctidae
# ============
# Lamarre suggests "Lichen" and "Plants" but a better coding
# that is more representative of the family would be "Phytophage"
# (xylophagous) and "Wood".
idx <- match("Xyloryctidae", D$Family)
D$Niche[idx] <- "Phytophage"
D$Habitat[idx] <- "Wood"

# Anostostomatidae
# ================
# Lamarre suggests "Saprophage" and "Soil". The group is omnivorous
# but "Saprophage" appears to be a reasonable coding. No update needed.


# Additional corrections
# ======================

# Pemphigidae
# ===========
# Now included in Aphididae, which is already coded,
# so just delete this entry
D <- D[!D$Family=="Pemphigidae",]

# Subfamilies
# ===========
# Subfamilies of Cecidomyiidae, Braconidae and Ichneumonidae
# coded by Lamarre but more accurately coded in collaboration
# with experts on these groups and in relation to the available
# reference phylogeny for these groups (see separate directories)
D <- D[!grepl("inae",D$Family),]

# Cynipoidea
# ==========
# The superfamily has been split (see Hearn et al. 2024), and the
# trait data need to be updated to reflect the reference tree.
D <- rbind(D,list(
    Class = "Arthropoda",
    Order = "Hymenoptera",
    Family = "Paraulacidae",
    Niche = "Phytophage-parasitoid",
    Habitat = "Plants",
    Length_min = 0.0,
    Length_max = 0.0))
D <- rbind(D,list(
    Class = "Arthropoda",
    Order = "Hymenoptera",
    Family = "Diplolepididae",
    Niche = "Phytophage",
    Habitat = "Plants",
    Length_min = 0.0,
    Length_max = 0.0))

# Check that coding is OK
check_family(D)
check_coding(D,"Family")


# Replace Entognatha with updated information
# ===========================================
cat("Replacing Entognatha data\n")
D <- D[D$Class!="Collembola",]
E <- read.delim("entognatha_taxonomy_traits.csv",sep=";")
E <- E[1:44,1:7]
colnames(E) <- c("Class","Order","Family","Niche","Habitat","Soil_depth","Size_range")
E$Niche[E$Niche=="Saprophagous"]<-"Saprophage"
E$Length_min <- 0.0
E$Length_max <- 0.0
a <- strsplit(E$Size_range, split="–")
for (i in 1:nrow(E)) {
    E$Length_min[i] <- a[[i]][1]
    if (grepl(" ",a[[i]][2]))
        E$Length_max[i] <- strsplit(a[[i]][2],split=" ")[[1]][1]
    else if (grepl("*",a[[i]][2]))
        E$Length_max[i] <- strsplit(a[[i]][2],split="*")[[1]][1]
    else
        E$Length_max[i] <- a[[i]][2]
}

# Sminthuridae
# ============
# Bruno and Neri code this family as "Water/Soil" for habitat. However, it appears
# that soil is by far the most common feeding habitat.
idx <- match("Sminthuridae",E$Family)
E$Habitat[idx] <- "Soil"

# Dicyrtomidae
# ============
# Bruno and Neri code this family as "Plant/Soil" for habitat. It appears that soil
# is the most common feeding habitat. This also agrees with the coding in R20 for SE.
idx <- match("Dicyrtomidae",E$Family)
E$Habitat[idx] <- "Soil"

# Add the clean data
E <- E[,column_names]
D <- rbind(D,E)

# Focus on the columns needed in the continued processing
# =======================================================
# Clade is the same as family for all taxa treated so far
# but this may not be true for the groups below
D$Clade <- D$Family
D <- D[,c("Order","Family","Clade","Niche","Habitat")]
check_coding(D)


# Replace Ichneumonidae data with updated information
# ===================================================
cat("Extending Ichneumonidae data\n")
D <- D[D$Family!="Ichneumonidae",]
E <- read.delim("ichneumonidae_traits.csv",sep=";")
E <- E[E$Family!="",]
E$Order <- "Hymenoptera"
E$Niche <- E[,which(colnames(E)=="Main.feeding.niche.MG")]
E$Habitat <- E[,which(colnames(E)=="Main.feeding.habitat.MG")]
E <- E[,colnames(D)]
E <- E[!duplicated(E),]
D <- rbind(D,E)
check_coding(D)


# Replace Braconidae data with updated information
# ================================================
cat("Extending Braconidae data\n")
D <- D[D$Family!="Braconidae",]
E <- read.delim("braconidae_traits.csv",sep=";")
E <- E[E$Family!="",]
E$Order <- "Hymenoptera"
E$Niche <- E[,which(colnames(E)=="Main.feeding.niche.MG")]
E$Habitat <- E[,which(colnames(E)=="Main.feeding.habitat.MG")]
E <- E[,colnames(D)]
E <- E[!duplicated(E),]
D <- rbind(D,E)
check_coding(D)


# Replace Cecidomyiidae data with updated information
# ===================================================
cat("Extending Cecidomyiidae data\n")
D <- D[D$Family!="Cecidomyiidae",]
E <- read.delim("cecidomyiidae_taxonomy_traits.csv",sep=";")
E <- E[E$Family!="",]
E$Order <- "Diptera"
E$Niche <- E[,which(colnames(E)=="Feeding.niche")]
E$Habitat <- E[,which(colnames(E)=="Feeding.habitat")]
E <- E[,colnames(D)]
E <- E[!duplicated(E),]
D <- rbind(D,E)
check_coding(D)


# Replace Chalcidoidea data with updated information
# ================================================
cat("Extending Chalcidoidea data\n")
# Hand code Chalcidoidea families in raw data
chalcidoid_fams <- c(
                    "Agaonidae",
                    "Aphelinidae",
                    "Chalcididae",
                    "Encyrtidae",
                    "Eucharitidae",
                    "Eulophidae",
                    "Eupelmidae",
                    "Eurytomidae",
                    "Megastigmidae",
                    "Mymaridae",
                    "Ormyridae",
                    "Perilampidae",
                    "Pteromalidae",
                    "Torymidae",
                    "Trichogrammatidae")
D <- D[!(D$Family %in% chalcidoid_fams),]
E <- read.delim("chalcidoidea_life_history_traits.csv",sep=";")
E <- E[E$Clade!="",]
E$Order <- "Hymenoptera"
E$Niche <- E[,which(colnames(E)=="Feeding.niche")]
E$Habitat <- E[,which(colnames(E)=="Feeding.microhabitat")]
E <- E[,colnames(D)]
E <- E[!duplicated(E),]
D <- rbind(D,E)
check_coding(D)


# Add entries for additional families detected in most recent annotation effort
# =============================================================================


# Remove NAs
# ==========
D$Niche[is.na(D$Niche)] <- ""
D$Habitat[is.na(D$Habitat)] <- ""
D$Niche[D$Niche=="?"] <- ""
D$Habitat[D$Habitat=="?"] <- ""


# Compute info on community and trophic level
D$Community <- unlist(lapply(D$Niche, FUN=community))
D$Trophic_level <- unlist(lapply(D$Niche, FUN=level))

write.table(D, "clade_trait_data_mg.tsv", row.names=FALSE, sep="\t")

