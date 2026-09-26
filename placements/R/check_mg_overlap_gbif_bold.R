# Check if any of the IBA species encountered in Madagascar
# are recorded in GBIF as occurring in Madagascar

# Check both against all species and against barcoded species
# in GBIF recorded from Madagascar

library(data.table)

T <- readRDS("../../iba_data/cluster_taxonomy_mg.rds")
G <- read.delim("../source/gbif_mg_species_list.tsv")

# Check # IBA OTUs from MG identified to species with SINTAX
x <- T$Species.sintax[!grepl(".",T$Species.sintax,fixed=TRUE) & !grepl("_X",T$Species.sintax) & !grepl("_sp%",paste0(T$Species.sintax,"%"))]
x <- gsub("_"," ",x)

# Check # species recorded from MG in GBIF
y <- G$species[grepl(" ",G$species)]

cat("Total # Madagascan OTUs retrieved by IBA:", nrow(T), "\n")
cat("Of these, identified to species level [Sintax]:", length(x), "\n")
cat("Of these, identified to unique species [Sintax]:", length(unique(x)), "\n")
cat("Total # Madagascan species recorded in GBIF:", length(y), "\n")
cat("Of these, with unique species names:", length(unique(y)), "\n")

# Prune BOLD data table down to what we need (and can manage in R)

# Read in only the desired columns
# TODO: adjust path to bold dump
ranks <- c("phylum","species")
headers <- colnames(fread("../../bold/BOLD_Public.31-Jul-2026.tsv", nrows=0))
col_idx <- which(headers %in% c("bin_uri",ranks))
B <- fread("../../bold/BOLD_Public.31-Jul-2026.tsv",select=col_idx,quote="")

# Constrain to species in Arthropoda
B <- B[B$phylum=="Arthropoda",]
D <- B[B$species %in% unique(y),]
D <- D[!duplicated(D$bin_uri),]

# Resolve taxonomy by 80% consensus rule (used by GBIF and IBA)
resolve_taxonomy <- function(B, D, ranks) {
    B <- data.frame(B)
    res <- data.frame(B[FALSE,])
    for (i in 1:nrow(D)) {
        if (i%%100==0)
            cat("Processed",floor(100*i/nrow(D)),"%\n")
        bin <- D$bin_uri[i]
        T <- B[B$bin_uri == bin,]
        if (nrow(T) == 0)
            next
        vals <- bin
        for (j in 1:length(ranks)) {
            idx <- which(colnames(T)==ranks[j])
            x <- table(T[,idx])
            x <- x[order(x,decreasing=TRUE)]
            if (as.numeric(x[1])/nrow(T) >= 0.8)
                vals[j+1] <- names(x)[1]
            else
                vals[j+1] <- "unresolved"
        }
        names(vals) <- colnames(res)
        res <- rbind(res,t(vals))   # R assumes column vectors
    }
    return (res)
}

res <- resolve_taxonomy(B, D, ranks)

cat("Of these, barcoded:", sum(y %in% res$species), "\n")
cat("Of the unique names, barcoded:", sum(unique(y) %in% res$species), "\n")


