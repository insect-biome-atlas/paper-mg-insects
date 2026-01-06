# Script for extracting relevant classification information from
# NCBI taxonomy dump files.

# We need name info with classification to major ranks for
# interpreting ncbi matches.

# We are interested in the Subfamily, Supertribe, Tribe and
# Species group information for the genera of arthropods, for
# taxonomic annotation purposes and for scoring trait data.
# We also need Class and Order information.

# Early experiments showed that the nodes.dmp file does not contain the full
# tree; it contains tribe and subtribe but not supertribe. Go figure...

# We therefore use the fullnamelineage file, which does seem to contain the
# full classification info

# Read in library for handling large files
library(data.table)

# Read in all arthropod taxa, with mandatory ranks
cat("Reading 'ncbi_taxonomy_dump/rankedlineage.dmp'\n")
X <- fread("ncbi_taxonomy_dump/rankedlineage.dmp")
X <- X[,c(1,3,5,7,9,11,13,15,17,21)]
colnames(X) <- c("TaxonID","Scientific.Name","Species","Genus","Family","Order","Class","Phylum","Kingdom","Domain")

# Solve problem with escaped double quotes (replace with single quotes)
X[] <- lapply(X, function(col) {
  if (is.character(col))
    gsub('\\"', "'", col, perl = TRUE) else col
})

# Get info on rank of each scientific name
cat("Reading 'ncbi_taxonomy_dump/nodes.dmp'\n")
Y <- fread("ncbi_taxonomy_dump/nodes.dmp")
Y <- Y[,c(1,5)]
colnames(Y) <- c("TaxonID","Rank")

# Complement the table with these
name_rank <- Y$Rank[match(X$TaxonID,Y$TaxonID)]
idx <- which(name_rank=="domain")
X$Domain[idx]<-X$Scientific.Name[idx]
idx <- which(name_rank=="kingdom")
X$Kingdom[idx]<-X$Scientific.Name[idx]
idx <- which(name_rank=="phylum")
X$Phylum[idx]<-X$Scientific.Name[idx]
idx <- which(name_rank=="class")
X$Class[idx]<-X$Scientific.Name[idx]
idx <- which(name_rank=="order")
X$Order[idx]<-X$Scientific.Name[idx]
idx <- which(name_rank=="family")
X$Family[idx]<-X$Scientific.Name[idx]
idx <- which(name_rank=="genus")
X$Genus[idx]<-X$Scientific.Name[idx]
idx <- which(name_rank=="species")
X$Species[idx]<-X$Scientific.Name[idx]

# Save full name table with ranked classification
write.table(X,"ncbi_name_ranked_classification.tsv", sep="\t", row.names=FALSE)

# Extract arthropod part
D <- data.frame(X[X$Phylum=="Arthropoda",])
rm(X)
rm(Y)

# Read in the full lineage information
cat("Reading 'ncbi_taxonomy_dump/fullnamelineage.dmp'\n")
X <- fread("ncbi_taxonomy_dump/fullnamelineage.dmp")
X <- X[,c(1,3,5)]
colnames(X) <- c("TaxonID", "Scientific.name", "Lineage")
F <- data.frame(X[grepl("Arthropoda",X$Lineage),])
rm(X)

# Assemble new data frame with info about ranks between Genus and Family
genera <- unique(D$Genus[D$Genus!=""])
idx <- match(genera,D$Scientific.Name)
E <- data.frame(list(TaxonID=D$TaxonID[idx], Genus=genera, Class=D$Class[idx], Order=D$Order[idx], Family=D$Family[idx]))
E$Subfamily <- ""
E$Supertribe <- ""
E$Tribe <- ""
E$Subtribe <- ""

cat("Processing genera\n")
for (i in 1:nrow(E)) {

    genus <- E$Genus[i]

    if (i%%1000==0) {
        cat(round(100*i/nrow(E),2),"%\n")
    }
    
    idx <- match(genus, F$Scientific.name)

    if (is.na(idx)) {
        cat("Warning: Genus '",genus,"' not found in fullnamelineage table\n", sep="")
        next;
    }

    family <- E$Family[i]
    
    if (family=="" || is.na(family)) {
        cat("Warning:Genus '",genus,"' does not have a family name\n", sep="")
        next;
    }

    a <- paste0(F$Lineage[idx], " ")
    b <- strsplit(a, "; ")[[1]]

    start <- which(b==family)
    if (length(start) != 1) {
        cat("Warning: Problem finding family '",family,"' in the full lineage of genus '", genus, "'\n", sep="");
        cat("   The full lineage is: ", b, ";\n", sep="; ")
        next;
    }

    for (j in start:length(b)) {
        name <- b[j]
        len <- nchar(name)
        if (substr(name, len-3, len)=="inae")
            E$Subfamily[i] <- name
        else if (substr(name, len-2, len)=="idi")
            E$Supertribe[i] <- name
        else if (substr(name, len-2, len)=="ini")
            E$Tribe[i] <- name
        else if (substr(name, len-2, len)=="ina")
            E$Subtribe[i] <- name
    }
}

write.table(E,"ncbi_arthropod_genus_classification_detailed.tsv", sep="\t", row.names=FALSE)

