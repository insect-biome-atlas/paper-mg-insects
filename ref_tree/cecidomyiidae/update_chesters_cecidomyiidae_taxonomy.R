source("../code/taxonomy_fxns.R")

T <- read.delim("chesters_cecidomyiidae_taxonomy.tsv")
X <- read.delim("../ncbi/ncbi_arthropod_genus_classification_detailed.tsv")

outfile <- "chesters_cecidomyiidae_taxonomy_detailed.tsv"

T$Subfamily <- X$Subfamily[match(T$Genus,X$Genus)]
T$Supertribe <- X$Supertribe[match(T$Genus,X$Genus)]

# Sort columns in right order according to taxonomic ranks
T <- T[,c("TipLabel","Kingdom","Phylum","Class","Order","Family","Subfamily","Supertribe","Genus","Species","Start","Stop")]

# Write table
write.table(T,outfile,sep="\t",row.names=FALSE)

