# Script for processing the traits data from Ronquist et al 2020

# Read in help functions
source("trait_fxns.R")

# Read in and process the data
# ============================
D  <- read.delim("ronquist_2020_SE_traits.csv",sep=";")

# Remove empty rows (if still in the data)
D <- D[!is.na(D$Sorting.Number),]

# Split combined taxon names into family and subfamily names
x <- character(nrow(D))
for (i in 1:nrow(D)) { word <- D$Taxon..Dyntaxa.2017.[i]; if (grepl("-",word)) { x[i] <- strsplit(word,"-")[[1]][1] } else x[i] <- word; }
D$Family <- x

for (i in 1:nrow(D)) { word <- D$Taxon..Dyntaxa.2017.[i]; if (grepl("-",word)) { x[i] <- strsplit(word,"-")[[1]][2] } else x[i] <- ""; }
D$Subfamily <- x

# Match family names to NCBI classification

# The below list is based on manual analyses of constituent taxa
# of family names that do not match the NCBI classification

D$NCBI_Family <- D$Family

ncbi_misses <- c("Cyphoderidae",
                 "Bovicolidae",
                 "Gliricolidae",
                 "Goniodidae",
                 "Acanthococcidae",
                 "Cryptococcidae",
                 "Drepanosiphidae",
                 "Lachnidae",
                 "Matsucoccidae",
                 "Pemphigidae",
                 "Rhizoecidae",
                 "Ulopidae",
                 "Plataspididae",
                 "Rhysodidae",
                 "Dasytidae",
                 "Drilidae",
                 "Malachiidae",
                 "Rhipiphoridae",
                 "Dryophthoridae",
                 "Megalopodidae",
                 "Nanophyidae",
                 "Rhynchitidae",
                 "Heptamelidae",
                 "Ismaridae",
                 "Sparasionidae",
                 "Myrmosidae",
                 "Chimabachidae",
                 "Meessidae",
                 "Parametriotidae",
                 "Peleopodidae",
                 "Mycetobiidae",
                 "Borboropsidae",
                 "Chiropteromyzidae",
                 "Coenomyiidae",
                 "Gasterophilidae",
                 "Hypodermatidae",
                 "Opetiidae",
                 "Phaeomyiidae",
                 "Stenomicridae",
                 "Trixoscelididae"
                )

ncbi_names <- c("Paronellidae",
                 "Bovicoliidae",
                 "Gyropidae",
                 "Philopteridae",
                 "Eriococcidae",
                 "Eriococcidae",
                 "Aphididae",
                 "Aphididae",
                 "Margarodidae",
                 "Aphididae",
                 "Pseudococcidae",
                 "Cicadellidae",
                 "Plataspidae",
                 "Carabidae",
                 "Melyridae",
                 "Elateridae",
                 "Melyridae",
                 "Ripiphoridae",
                 "Curculionidae",
                 "Chrysomelidae",
                 "Brentidae",
                 "Attelabidae",
                 "Tenthredinidae",
                 "Diapriidae",
                 "Scelionidae",
                 "Mutillidae",
                 "Lypusidae",
                 "Meessiidae",
                 "Elachistidae",
                 "Depressariidae",
                 "Anisopodidae",
                 "Heleomyzidae",
                 "Heleomyzidae",
                 "Xylophagidae",
                 "Oestridae",
                 "Oestridae",
                 "Platypezidae",
                 "Sciomyzidae",
                 "Aulacigastridae",
                 "Heleomyzidae"
                )

D$NCBI_Family[match(ncbi_misses,D$NCBI_Family)] <- ncbi_names

# # Print data for mismatched ncbi names for manual checking
# a <- unique(ncbi_names)
# for (i in 1:length(a)) {
#    cat ("Family: ", a[i], "\n")
#    print(D[D$NCBI_Family %in% a[i],])
# }

# Correct coding
D$Main.feeding.niche[D$Main.feeding.niche=="Saprophagous"] <- "Saprophage"
D$Main.feeding.niche[D$Main.feeding.niche=="Phytophagous"] <- "Phytophage"
D$Main.feeding.habitat[D$Main.feeding.habitat=="Temporary Habitats"]<-"Temporary habitats"

# Majority rule resolution of conflicts

# Only one species in Rhysodidae, many in Carabidae: use data for Carabidae
D$Main.feeding.niche[match("Carabidae",D$NCBI_Family)] <- "Predator"
D$Main.feeding.habitat[match("Carabidae",D$NCBI_Family)] <- "Soil"

# Slightly more species in Dasytidae (19) than Malachiidae (15): use data for Dasytidae
D$Main.feeding.habitat[match("Melyridae",D$NCBI_Family)] <- "Wood"

# Only two species in Drilidae, many in Elateridae: use data for Elateridae
D$Main.feeding.niche[match("Elateridae",D$NCBI_Family)] <- "Saprophage"

# Only a handful of species in Ismaridae, many in Diapriidae: use data for Diapriidae
D$Main.feeding.niche[match("Diapriidae",D$NCBI_Family)] <- "Saprophage-parasitoid"
D$Main.feeding.habitat[match("Diapriidae",D$NCBI_Family)] <- "Temporary habitats"

# Only 3 species in Sparasionidae, many in Scelionidae: use data for Scelionidae
D$Main.feeding.niche[match("Scelionidae",D$NCBI_Family)] <- "Predator-parasitoid"

# Only 3 species in Mycetobiidae, 7 in Anisopodidae: use data for Anisopodidae
D$Main.feeding.habitat[match("Anisopodidae",D$NCBI_Family)] <- "Temporary habitats"

# Only 1 species in Opetiidae, 38 in Platypezidae: use data for Platypezidae
D$Main.feeding.habitat[match("Platypezidae",D$NCBI_Family)] <- "Fungi"

# Only 2 species in Phaeomyiidae, many in Sciomyzidae: use data for Sciomyzidae
D$Main.feeding.habitat[match("Sciomyzidae",D$NCBI_Family)] <- "Water"

# 2 species in Aulacigastridae, 3 in Stenomicridae: use data for Stenomicridae
D$Main.feeding.habitat[match("Aulacigastridae",D$NCBI_Family)] <- "Water"

# Staphylinidae is dominated by Predator and Soil species (Saprophagous and Temporary habitats also common)
D$Main.feeding.niche[match("Staphylinidae",D$NCBI_Family)] <- "Predator"
D$Main.feeding.habitat[match("Staphylinidae",D$NCBI_Family)] <- "Soil"

# Only keep first NCBI family record
D <- D[!duplicated(D$NCBI_Family),]

# Write updated table (temporary result)
write.table(D, "ronquist_2020_SE_traits_ncbi_taxonomy.tsv",row.names=FALSE,sep="\t")

# Prune down to essential columns
D <- D[,c("Order","NCBI_Family","Main.feeding.niche","Main.feeding.habitat")]
colnames(D) <- c("Order","Family","Niche","Habitat")

# Expand/replace info for key taxa

# Cynipoidea
# ==========
# The superfamily has been split (see Hearn et al. 2024), and the
# trait data need to be updated to reflect the reference tree.
D <- rbind(D,list(
    Order = "Hymenoptera",
    Family = "Diplolepididae",
    Niche = "Phytophage",
    Habitat = "Plants"))

# Eriaporidae is now considered part of Pirenidae, which is not
# included otherwise in the data
D$Family[D$Family=="Eriaporidae"] <- "Pirenidae"

# A new family Athaliidae is now recognized in the Tenthredinoidea
D <- rbind(D,list(
    Order = "Hymenoptera",
    Family = "Athaliidae",
    Niche = "Phytophage",
    Habitat = "Plants"))

# Ectobius now placed in a family separate from Blattellidae (Ectobiidae)
D <- rbind(D,list(
    Order = "Blattodea",
    Family = "Ectobiidae",
    Niche = "Saprophage",
    Habitat = "Temporary habitats"))


# Add Clade with default being family
D$Clade <- D$Family

# Check that coding is OK
check_family(D)
check_coding(D,"Family")


# Replace Entognatha with updated information
# ===========================================
cat("Replacing Entognatha data\n")
entognatha <- c("Collembola","Diplura","Protura")
D <- D[!(D$Order %in% entognatha),]
E <- read.delim("entognatha_taxonomy_traits.csv",sep=";")
E <- E[1:44,2:5]
colnames(E) <- c("Order","Family","Niche","Habitat")
E$Clade <- E$Family
idx <- match("Gulgastruridae",E$Family)
E$Habitat[idx] <- "Soil"    # Given as Soil (cave entrance)
idx <- match("Saprophagous",E$Niche)
E$Niche[idx] <- "Saprophage"    # Correct to state name used elsewhere

# Sminthurididae
# ==============
# Bruno and Neri code this family as "Water/Soil" for habitat. However, it appears
# that soil is by far the most common feeding habitat.
idx <- match("Sminthurididae",E$Family)
E$Habitat[idx] <- "Soil"

# Dicyrtomidae
# ============
# Bruno and Neri code this family as "Plant/Soil" for habitat. It appears that soil
# is the most common feeding habitat. This also agrees with the coding in R20 for SE.
idx <- match("Dicyrtomidae",E$Family)
E$Habitat[idx] <- "Soil"

# Add the clean data
E <- E[,colnames(D)]
D <- rbind(D,E)


# Replace Ichneumonidae data with updated information
# ===================================================
cat("Extending Ichneumonidae data\n")
D <- D[D$Family!="Ichneumonidae",]
E <- read.delim("ichneumonidae_traits.csv",sep=";")
E <- E[E$Family!="",]
E$Order <- "Hymenoptera"
E$Niche <- E[,which(colnames(E)=="Main.feeding.niche.SE")]
E$Habitat <- E[,which(colnames(E)=="Main.feeding.habitat.SE")]
E <- E[,colnames(D)]
E <- E[!duplicated(E),]     # Note that any ambiguous coding of clades will remain
idx <- match("Hybrizontinae",E$Clade)
E$Clade[idx] <- "Hybrizoninae"      # Use NCBI name despite being incorrect...
D <- rbind(D,E)
check_coding(D)


# Replace Braconidae data with updated information
# ================================================
cat("Extending Braconidae data\n")
D <- D[D$Family!="Braconidae",]
E <- read.delim("braconidae_traits.csv",sep=";")
E <- E[E$Family!="",]
E$Order <- "Hymenoptera"
E$Niche <- E[,which(colnames(E)=="Main.feeding.niche.SE")]
E$Habitat <- E[,which(colnames(E)=="Main.feeding.habitat.SE")]
E <- E[,colnames(D)]
E <- E[!duplicated(E),]     # Note that any ambiguous coding of clades will remain
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
E <- E[!duplicated(E),]     # Note that any ambiguous coding of clades will remain
D <- rbind(D,E)
check_coding(D)


# Replace Chalcidoidea data with updated information
# ================================================
cat("Extending Chalcidoidea data\n")
# Hand code the Chalcidoidea fams in R2020
chalcidoid_fams <- c("Aphelinidae",
                "Azotidae",
                "Chalcididae",
                "Encyrtidae",
                "Pirenidae", # Eriaporidae
                "Eulophidae",
                "Eupelmidae",
                "Eurytomidae",
                "Mymaridae",
                "Ormyridae",
                "Perilampidae",
                "Pteromalidae",
                "Signiphoridae",
                "Tetracampidae",
                "Torymidae",
                "Trichogrammatidae")
D <- D[!(D$Family %in% chalcidoid_fams),]
E <- read.delim("chalcidoidea_life_history_traits.csv",sep=";")
E <- E[E$Clade!="",]
E$Order <- "Hymenoptera"
E$Niche <- E[,which(colnames(E)=="Feeding.niche.SE")]
E$Habitat <- E[,which(colnames(E)=="Feeding.microhabitat.SE")]
E <- E[,colnames(D)]
E <- E[!duplicated(E),]     # Note that any ambiguous coding of clades will remain
D <- rbind(D,E)
check_coding(D)


# Add entries for additional families detected in most recent annotation effort
# =============================================================================


# Remove NAs and ?
# ================
D$Niche[is.na(D$Niche)] <- ""
D$Habitat[is.na(D$Habitat)] <- ""
D$Niche[D$Niche=="?"] <- ""
D$Habitat[D$Habitat=="?"] <- ""


# Add community and trophic level info
D$Community <- unlist(lapply(D$Niche, FUN=community))
D$Trophic_level <- unlist(lapply(D$Niche, FUN=level))

write.table(D, "clade_trait_data_se.tsv", row.names=FALSE, sep="\t")

