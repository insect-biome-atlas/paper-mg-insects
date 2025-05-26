source("../code/taxonomy_fxns.R")

T <- read.delim("chesters_cecidomyiidae_taxonomy.tsv")

outfile <- "chesters_cecidomyiidae_taxonomy_detailed.tsv"

res <- get_subfam_supertribe(T$Species)

T$Subfamily <- res$Subfamily
T$Supertribe <- res$Supertribe

# Sort columns in right order according to taxonomic ranks
T <- T[,c("TipLabel","Kingdom","Phylum","Class","Order","Family","Subfamily","Supertribe","Genus","Species","Start","Stop")]

# Write table
write.table(T,outfile,sep="\t",row.names=FALSE)

