# Prepare data tables for faunal composition accumulation plots
# using all Swedish data (forests and all other habitats included)

# Read in data
# ------------

# Cluster taxonomy (only Hexapoda, and with quality filtering)
se_clusters <- read.delim("../../iba_data/cluster_taxonomy_se.tsv")

# Add info on life-history traits
lht_se <- read.delim("../../traits/clade_trait_data_se.tsv")

# Merge cluster taxonomy and life-history traits
se_clusters <- merge(se_clusters, lht_se)

# Read in cluster read numbers for malaise traps and litter samples
se_malaise_counts_long <- readRDS("../../iba_data/cluster_counts_malaise_long_se.rds")
se_litter_counts_long  <- readRDS("../../iba_data/cluster_counts_litter_long_se.rds")

# Merge cluster read numbers
se_counts_long <- rbind(se_malaise_counts_long, se_litter_counts_long)


# Tidy data
# ---------

# Merge counts with cluster taxonomy
taxa_counts_se <- merge(se_counts_long, se_clusters, by="cluster")

# Read in collected sample metadata
sample_meta_se <- read.delim("../../iba_data/malaise_litter_sample_meta_se.tsv")


# Combine metadata and cluster data into otu table at sample level
# ----------------------------------------------------------------
otu_sample_meta_se <- merge(taxa_counts_se, sample_meta_se, by="sampleID_NGI")

# Save data file
saveRDS(otu_sample_meta_se,"../data/otu_sample_meta_se.rds")


# Aggregate counts at trap level
# ------------------------------

otu_site_meta_se <- aggregate(read_count~trapID+cluster+Order+Niche+Habitat+trap_habitat+latitude+longitude, data=otu_sample_meta_se, FUN=sum)

# Save data file
saveRDS(otu_site_meta_se,"../data/otu_site_meta_se.rds")

