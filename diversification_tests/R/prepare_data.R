# Prepare data for random diversification composition tests


# Read in data
# ------------

# Read in placements
P <- readRDS("../../placements/data/placement_stats.rds")

# Cluster taxonomy (only Hexapoda, and with quality filtering)
clusters <- readRDS("../../iba_data/cluster_taxonomy_mg.rds")

# Add info on life-history traits
lht <- read.delim("../../traits/clade_trait_data_mg.tsv")

# Merge cluster taxonomy and life-history traits
clusters <- merge(clusters, lht)


# Tidy data
# ---------

D <- aggregate(P$placement_age,by=list(edge_num=P$edge_num),FUN=mean)
colnames(D) <- c("edge_num","placement_age")
X <- table(P$edge_num)
D$OTUs <- as.numeric(X[match(D$edge_num,names(X))])
D$ASV <- P$cluster_rep[match(D$edge_num,P$edge_num)]
D <- D[D$ASV %in% clusters$ASV,]    # Remove non-matching clusters (those lacking life-history data)
idx <- match(D$ASV,clusters$ASV)
D$cluster <- clusters$cluster[idx]
D$Clade <- clusters$Clade[idx]
D$Family <- clusters$Family[idx]
D$Order <- clusters$Order[idx]
D$Niche <- clusters$Niche[idx]
D$Habitat <- clusters$Habitat[idx]
D$Coleoptera <- 1*(D$Order=="Coleoptera")

# Save data for post-processing
write.table(D,"../data/diversification_data.tsv",row.names=FALSE,sep="\t")

# Save data for treeppl inference
cat('{"ages":[', file='../data/div_data.json')
cat(D$placement_age, sep=',\n', file='../data/div_data.json', append=TRUE)
cat('],\n"otus":[', file='../data/div_data.json', append=TRUE)
cat(D$OTUs, sep=',\n', file='../data/div_data.json', append=TRUE)
cat(']}\n', file='../data/div_data.json', append=TRUE)

# Save data for treeppl inference with bias correction
cat('{"ages":[', file='../data/div_data_bias.json')
cat(D$placement_age[D$Order!="Coleoptera"], sep=',\n', file='../data/div_data_bias.json', append=TRUE)
cat('],\n"otus":[', file='../data/div_data_bias.json', append=TRUE)
cat(D$OTUs[D$Order!="Coleoptera"], sep=',\n', file='../data/div_data_bias.json', append=TRUE)
cat('],\n"bias_ages":[', file='../data/div_data_bias.json', append=TRUE)
cat(D$placement_age[D$Order=="Coleoptera"], sep=',\n', file='../data/div_data_bias.json', append=TRUE)
cat('],\n"bias_otus":[', file='../data/div_data_bias.json', append=TRUE)
cat(D$OTUs[D$Order=="Coleoptera"], sep=',\n', file='../data/div_data_bias.json', append=TRUE)
cat(']}\n', file='../data/div_data_bias.json', append=TRUE)

