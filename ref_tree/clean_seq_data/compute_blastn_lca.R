# Compute lca assignments from blastn searches of Chesters sequences

# Read in functions needed
source("lca_fxns.R")

# Read in taxonomic data for expanded tree
cat("Reading in taxonomic data for expanded tree (raw data)\n")
C <- read.delim("chesters_expanded_taxonomy_raw.tsv")

# Read in blastn results
cat("Reading in blastn.tsv file\n")
B <- read.delim("blastn.tsv", header=FALSE)
colnames(B) <- c("query","accession","per.ident","acc.length","differences","gaps","query.start","query.stop","sbjct.start","sbjct.stop","e.value","score","taxon.id")

# Restrict to the relevant hits
B <- B[B$per.ident >= 95.0,]

# Read in ncbi taxonomy info
cat("Reading in ncbi taxonomy info\n")
T <- read.delim("../ncbi/ncbi_name_ranked_classification.tsv")
# NB! This file does not include the scientific name in the ranked lineage info,
# even when the name belongs to one of the recognized ranks. This means that the
# code below recoqnizes the least common ancestor down to the level just above
# the lca, if the lca is a taxon belonging to one of the recognized ranks. This
# is no problem here, though, as we are looking for mismatches at the order or
# class level. If the scientific name looks like a species name, however, we
# add it as a species name.

# Create an empty data frame for the data
D <- data.frame(list(query=character(),
                     Domain=character(),
                     Kingdom=character(),
                     Phylum=character(),
                     Class=character(),
                     Order=character(),
                     Family=character(),
                     Genus=character(),
                     Species=character()))

misses <- data.frame(list(TipLabel=C$TipLabel[!(C$TipLabel %in% unique(B$query))]))
# cat("There are", nrow(misses), "sequences without matches >= 95.0 per.ident\n")
write.table(misses,"no_close_match.tsv",sep="\t",row.names=FALSE)

cat("Now processing names\n")
for (i in 1:nrow(C)) {
    query <- C$TipLabel[i]
    if (i %% 100 == 0)
        cat("Processed",round(100*i/nrow(C),2),"%\n")
#    cat("Now processing query:",query,"\n")
    if (!(query %in% B$query)) {
        D <- rbind(D,list(query=query,
                          Domain="None",
                          Kingdom="None",
                          Phylum="None",
                          Class="None",
                          Order="None",
                          Family="None",
                          Genus="None",
                          Species="None"))
    } else {
        taxon_ids <- get_taxon_ids(B$taxon.id[B$query==query])
#        cat("The taxon ids are", taxon_ids, "\n")
        x <- lca(T,taxon_ids)
#        cat("The lca list is\n")
#        print(x)
        D <- rbind(D, c(list(query=query),x))
    }
}

write.table(D,"blastn_lca_0_95_idty.tsv", row.names=FALSE, sep="\t")

