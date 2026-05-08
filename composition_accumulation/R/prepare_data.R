# Prepare data tables for faunal composition accumulation plots

# Read in data
# ------------

# Cluster taxonomy (only Hexapoda, and with quality filtering)
mg_clusters <- readRDS("../../iba_data/cluster_taxonomy_mg.rds")
se_clusters <- readRDS("../../iba_data/cluster_taxonomy_se.rds")

# Add info on life-history traits
lht_mg <- read.delim("../../traits/clade_trait_data_mg.tsv")
lht_se <- read.delim("../../traits/clade_trait_data_se.tsv")

# Merge cluster taxonomy and life-history traits
mg_clusters <- merge(mg_clusters, lht_mg)
se_clusters <- merge(se_clusters, lht_se)

# Read in cluster read numbers for malaise traps and litter samples
mg_malaise_counts_long <- readRDS("../../iba_data/cluster_counts_malaise_long_mg.rds")
se_malaise_counts_long <- readRDS("../../../iba_data/cluster_counts_malaise_long_se.rds")
mg_litter_counts_long <- readRDS("../../iba_data/cluster_counts_litter_long_mg.rds")
se_litter_counts_long <- readRDS("../../iba_data/cluster_counts_litter_long_se.rds")

# Merge cluster read numbers
mg_counts_long <- rbind(mg_malaise_counts_long, mg_litter_counts_long)
se_counts_long <- rbind(se_malaise_counts_long, se_litter_counts_long)


# Tidy data
# ---------

# Merge counts with cluster taxonomy
taxa_counts_mg <- merge(mg_counts_long, mg_clusters, by="cluster")
taxa_counts_se <- merge(se_counts_long, se_clusters, by="cluster")

# Read in collected sample metadata
sample_meta_mg <- read.delim("../../iba_data/malaise_litter_sample_meta_mg.tsv")
sample_meta_se <- read.delim("../../iba_data/malaise_litter_sample_meta_se.tsv")


# Combine metadata and cluster data into otu table at sample level
# ----------------------------------------------------------------

otu_sample_meta_mg <- merge(taxa_counts_mg, sample_meta_mg, by="sampleID_NGI")
otu_sample_meta_se <- merge(taxa_counts_se, sample_meta_se, by="sampleID_NGI")

# Restrict Swedish data to forest sites
otu_sample_meta_se <- otu_sample_meta_se[otu_sample_meta_se$trap_habitat=="forest",]

# Save data files
write.table(otu_sample_meta_mg,"../data/otu_sample_meta_mg.tsv",row.names=FALSE,sep="\t")
write.table(otu_sample_meta_se,"../data/otu_sample_meta_se.tsv",row.names=FALSE,sep="\t")


# Aggregate counts at trap level
# ------------------------------

otu_site_meta_mg <- aggregate(read_count~trapID+cluster+Order+Niche+Habitat+trap_habitat+latitude+longitude, data=otu_sample_meta_mg, FUN=sum)
otu_site_meta_se <- aggregate(read_count~trapID+cluster+Order+Niche+Habitat+trap_habitat+latitude+longitude, data=otu_sample_meta_se, FUN=sum)

# Save data files
write.table(otu_site_meta_mg,"../data/otu_site_meta_mg.tsv",row.names=FALSE,sep="\t")
write.table(otu_site_meta_se,"../data/otu_site_meta_se.tsv",row.names=FALSE,sep="\t")

