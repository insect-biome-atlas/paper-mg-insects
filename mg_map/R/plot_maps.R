# Plot mg map with sample sites, diversity etc

library(vegan)
library(patchwork)
library(ggplot2)



# Read in data
# ------------

# Cluster taxonomy (only Hexapoda, and with quality filtering)
clusters <- read.delim("../../iba_data/cluster_taxonomy_mg.tsv")

# Add info on life-history traits
lht <- read.delim("../../traits/clade_trait_data_mg.tsv")

# Merge cluster taxonomy and life-history traits
clusters <- merge(clusters, lht)

# Read in cluster read numbers for malaise traps and litter samples
malaise_counts_long <- readRDS("../../iba_data/cluster_counts_malaise_long_mg.rds")
litter_counts_long <- readRDS("../../iba_data/cluster_counts_litter_long_mg.rds")

# Merge cluster read numbers
counts_long <- rbind(malaise_counts_long, litter_counts_long)


# Tidy data
# ---------

# Merge counts with cluster taxonomy
taxa_counts <- merge(counts_long, clusters, by="cluster")

# Read in collected sample metadata
sample_meta <- read.delim("../../iba_data/malaise_litter_sample_meta_mg.tsv")

# Extract trap meta
trap_meta <- unique(sample_meta[,c("trapID","latitude","longitude","trap_habitat")])

# Combine metadata and cluster data into otu table at sample level
otu_sample_meta <- merge(taxa_counts, sample_meta, by="sampleID_NGI")

# Aggregate counts at trap level
otu_site_meta <- aggregate(read_count~trapID+cluster+Order+Niche+Habitat+trap_habitat+latitude+longitude, data=otu_sample_meta, FUN=sum)


# Compute plot data
# -----------------

# Compute site species richness
site_catch <- data.frame(table(otu_site_meta$trapID))
colnames(site_catch) <- c("trapID","OTUs")
idx <- match(site_catch$trapID,trap_meta$trapID)
site_catch$latitude <- trap_meta$latitude[idx]
site_catch$longitude <- trap_meta$longitude[idx]
site_catch$trap_habitat <- trap_meta$trap_habitat[idx]

# Compute site diversity (Shannon)
site_otu_abundance <- readRDS("../../iba_data/site_otu_abundance_malaise_mg.rds")
shannon_diversity <- diversity(site_otu_abundance, index="shannon")
site_diversity <- data.frame(shannon_diversity)
site_diversity$trapID <- row.names(site_diversity)
idx <- match(site_diversity$trapID,trap_meta$trapID)
site_diversity$latitude <- trap_meta$latitude[idx]
site_diversity$longitude <- trap_meta$longitude[idx]
site_diversity$trap_habitat <- trap_meta$trap_habitat[idx]

# Compute site uniqueness (based on Bray-Curtis dissimilarities)
dist_matrix <- vegdist(site_otu_abundance, method="bray")
dist_mat <- as.matrix(dist_matrix)
x <- colSums(dist_mat) / (nrow(dist_mat) - 1) # Compute mean dissimilarity, exclude self-distance (=0)
site_uniqueness <- data.frame(list(trapID=names(x),dissimilarity=x))
idx <- match(site_uniqueness$trapID,trap_meta$trapID)
site_uniqueness$latitude <- trap_meta$latitude[idx]
site_uniqueness$longitude <- trap_meta$longitude[idx]
site_uniqueness$trap_habitat <- trap_meta$trap_habitat[idx]

# Compute site 'colonizations' (number of endemic clades)
P <- readRDS("../../placements/data/placement_stats.rds")
D <- unique(otu_site_meta[,c("cluster","trapID")])
D$ASV <- clusters$ASV[match(D$cluster,clusters$cluster)]
D <- D[D$ASV %in% P$cluster_rep,]
D$placement <- P$edge_num[match(D$ASV,P$cluster_rep)]
x <- table(P$edge_num)
x <- x[x > 4]   # Threshold for an MG radiation 
D <- D[D$placement %in% names(x),]
E <- data.frame(table(D$placement,by=D$trapID))
colnames(E) <- c("Placement","trapID","num_placements")
E <- E[E$num_placements!=0,]
site_colonizations <- data.frame(table(E$trapID))
colnames(site_colonizations) <- c("trapID","num_clades")
idx <- match(site_colonizations$trapID,trap_meta$trapID)
site_colonizations$latitude <- trap_meta$latitude[idx]
site_colonizations$longitude <- trap_meta$longitude[idx]
site_colonizations$trap_habitat <- trap_meta$trap_habitat[idx]

# Compute site 'colonization ages' (age of endemic clade)
X <- aggregate(P$placement_age,by=list(edge_num=P$edge_num),FUN=mean)
colnames(X) <- c("edge_num","placement_age")
F <- unique(D[,c("placement","trapID")])
F$placement_age <- X$placement_age[match(F$placement,X$edge_num)]
site_radiation_age <- aggregate(placement_age~trapID, data=F, FUN=median)
colnames(site_radiation_age) <- c("trapID","median_age")
idx <- match(site_radiation_age$trapID,trap_meta$trapID)
site_radiation_age$latitude <- trap_meta$latitude[idx]
site_radiation_age$longitude <- trap_meta$longitude[idx]
site_radiation_age$trap_habitat <- trap_meta$trap_habitat[idx]

# Compute site 'endemicity' (median endemic clade size)
site_endemicity <- aggregate(num_placements~trapID, data=E, FUN=mean)
colnames(site_endemicity) <- c("trapID","mean_radiation")
idx <- match(site_endemicity$trapID,trap_meta$trapID)
site_endemicity$latitude <- trap_meta$latitude[idx]
site_endemicity$longitude <- trap_meta$longitude[idx]
site_endemicity$trap_habitat <- trap_meta$trap_habitat[idx]


# Generate plots
# --------------

p1 <- ggplot(data=site_catch, aes(x=longitude, y=latitude, shape=trap_habitat, colour=OTUs)) +
        geom_point(size=3, position="jitter") +
#        scale_colour_viridis_c(option="plasma") +
        labs(title = "Richness",
             x = NULL,
             y = NULL,
             shape = "Habitat",
             colour = "# OTUs")

p2 <- ggplot(data=site_diversity, aes(x=longitude, y=latitude, shape=trap_habitat, colour=shannon_diversity)) +
        geom_point(size=3, position="jitter") +
#        scale_colour_viridis_c(option="plasma") +
        labs(title = "Diversity",
             x = NULL,
             y = NULL,
             shape = "Habitat",
             colour = "Shannon")

p3 <- ggplot(data=site_uniqueness, aes(x=longitude, y=latitude, shape=trap_habitat, colour=dissimilarity)) +
        geom_point(size=3, position="jitter") +
#        scale_colour_viridis_c(option="plasma") +
        labs(title = "Uniqueness",
             x = NULL,
             y = NULL,
             shape = "Habitat",
             colour = "Bray-Curtis")

p4 <- ggplot(data=site_colonizations, aes(x=longitude, y=latitude, shape=trap_habitat, colour=num_clades)) +
        geom_point(size=3, position="jitter") +
#        scale_colour_viridis_c(option="plasma") +
        labs(title = "Number of MG radiations",
             x = NULL,
             y = NULL,
             shape = "Habitat",
             colour = "# clades")

p5 <- ggplot(data=site_radiation_age, aes(x=longitude, y=latitude, shape=trap_habitat, colour=median_age)) +
        geom_point(size=3, position="jitter") +
#        scale_colour_viridis_c(option="plasma") +
        labs(title = "Median age of MG radiations",
             x = NULL,
             y = NULL,
             shape = "Habitat",
             colour = "Median age (Ma)")

p6 <- ggplot(data=site_endemicity, aes(x=longitude, y=latitude, shape=trap_habitat, colour=mean_radiation)) +
        geom_point(size=3, position="jitter") +
#        scale_colour_viridis_c(option="plasma") +
        labs(title = "Size of MG radiations",
             x = NULL,
             y = NULL,
             shape = "Habitat",
             colour = "Mean # OTUs")



# Put plots together and save
ggsave("../figs/Fig_mg_maps.jpg",
       width=18,
       height=14,
       plot = p1 + p2 + p3 + p4 + p5 + p6 +
            plot_layout(axis_titles="collect", ncol=3) +
            plot_annotation(tag_levels="A") #  & theme(legend.position="bottom")
       )

