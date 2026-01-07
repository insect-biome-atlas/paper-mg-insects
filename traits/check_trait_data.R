# Some basic checks of the trait data
# in relation to the reference tree used

# Read in data
D1 <- read.delim("clade_trait_data_se.tsv")
D2 <- read.delim("clade_trait_data_mg.tsv")
E <- read.delim("../ref_tree/expanded_tree/chesters_expanded_taxonomy.tsv")

entognatha <- c("Collembola","Protura","Diplura")

# We are only interested in unique clade info for entognatha or insecta
EE <- E[E$Class %in% entognatha,]
EI <- E[E$Class=="Insecta",]
EE <- EE[!duplicated(EE$Clade),]
EI <- EI[!duplicated(EI$Clade),]

# Check that all clades in tree have trait data for
# the expanded groups. Mismatches are errors
# =================================================

# Entognatha
cat("Mismatches for Entognatha in SE and MG\n")
print(EE[!(EE$Clade %in% D1$Clade),])
print(EE[!(EE$Clade %in% D2$Clade),])

# Ichneumonidae
cat("Mismatches for Ichneumonidae in SE and MG\n")
Eich <- EI[EI$Family=="Ichneumonidae",]
print(Eich[!(Eich$Clade %in% D1$Clade),])
print(Eich[!(Eich$Clade %in% D2$Clade),])

# Braconidae
cat("Mismatches for Braconidae in SE and MG\n")
Ebra <- EI[EI$Family=="Braconidae",]
print(Ebra[!(Ebra$Clade %in% D1$Clade),])
print(Ebra[!(Ebra$Clade %in% D2$Clade),])

# Chalcidoidea
cat("Mismatches for Chalcidoidea in SE and MG\n")
F <- read.delim("../ref_tree/chalcidoidea/expanded_chalcidoidea3_taxonomy.tsv")
cha_fams <- unique(F$Family[F$Family!="Mymarommatidae"])
Echa <- EI[EI$Family %in% cha_fams,]
print(Echa[!(Echa$Clade %in% D1$Clade),])
print(Echa[!(Echa$Clade %in% D2$Clade),])

# Cecidomyiidae
cat("Mismatches for Cecidoymiidae in SE and MG\n")
Ecec <- EI[EI$Family=="Cecidomyiidae",]
print(Ecec[!(Ecec$Clade %in% D1$Clade),])
print(Ecec[!(Ecec$Clade %in% D2$Clade),])

# List clades in trait data missing in tree
# =========================================
cols <- c("Order","Family","Clade")
foo <- function(D) {
    D <- D[!(D$Order %in% EE$Order),]
    D <- D[!(D$Clade %in% EI$Clade),cols]
    cat("There are",nrow(D),"missing entries\n")
    print(D[order(D$Order),])
}
cat("SE trait clades missing in tree\n")
foo(D1)
cat("MG trait clades missing in tree\n")
foo(D2)

# List clades in tree missing in trait data
# =========================================
foo <- function(D) {
    D <- D[!(D$Order %in% EE$Order),]
    EI <- EI[!(EI$Clade %in% D$Clade),cols]
    cat("There are",nrow(EI),"missing entries\n")
    print(EI[order(EI$Order),])
}
cat("Clades in tree missing SE trait data\n")
foo(D1)
cat("Clades in tree missing MG trait data\n")
foo(D2)


