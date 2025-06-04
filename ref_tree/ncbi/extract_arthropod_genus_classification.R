# Script for extracting relevant classification information from
# NCBI taxonomy dump files. We are interested in the Subfamily,
# Supertribe, Tribe and Species group information for the
# genera of arthropods.

# Read in library for handling large files
library(data.table)

# Read in all arthropod taxa, with mandatory ranks
X <- fread("rankedlineage.dmp")
X <- X[,c(1,3,5,7,9,11,13,15,17,21)]
colnames(X) <- c("TaxonID","Scientific.Name","Species","Genus","Family","Order","Class","Phylum","Kingdom","Domain")
D <- data.frame(X[X$Phylum=="Arthropoda",])
rm(X)

# Read in the tree node info
X <- fread("nodes.dmp")
X <- X[,c(1,3,5)]
colnames(X) <- c("TaxonID", "ParentID", "Rank")

# Assemble new data frame with info about ranks between Genus and Family
genera <- unique(D$Genus[D$Genus!=""])
E <- data.frame(list(TaxonID=D$TaxonID[match(genera,D$Scientific.Name)], Genus=genera, Family=D$Family[match(genera,D$Scientific.Name)]))
E$Subfamily <- ""
E$Supertribe <- ""
E$Tribe <- ""
E$Subtribe <- ""
E$Family <- ""

for (i in 1:nrow(E)) {

    if (i%%100==0)
        cat(i,"\n")
    id <- E$TaxonID[i]
    rank <- X$Rank[match(id,X$TaxonID)]
    while (rank != "Family" && id != 6656) {    # TaxonID 6656 is "Arthropoda", catches cases where family is missing

        id <- X$ParentID[match(id,X$TaxonID)]
        rank <- X$Rank[match(id,X$TaxonID)]
        name <- D$Scientific.Name[match(id,D$TaxonID)]
        if (rank=="subtribe") E$Subtribe[i] <- name
        if (rank=="tribe") E$Tribe[i] <- name
        if (rank=="supertribe") E$Supertribe[i] <- name
        if (rank=="subfamily") E$Subfamily[i] <- name
        if (rank=="family") E$Family[i] <- name
    }
}

write.table(E,"ncbi_arhtropod_genus_classification_detailed.tsv")

