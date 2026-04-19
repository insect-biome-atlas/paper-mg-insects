# Generate accumulation plots

library(vegan)
library(ggplot2)
library(patchwork)

# Set basic plot params
bs    <- 20
lsize <- 2

# Read data
# ---------

cat("Reading data\n")

# Read species * site matrix
sp_matrix_mg <- readRDS("../../iba_data/site_otu_occurrence_combined_mg.rds")
sp_matrix_mg <- 1*sp_matrix_mg  # Convert to numeric

# Read in sample metadata
site_meta_mg <- read.delim("../../iba_data/malaise_litter_sample_meta_mg.tsv")

# Read in cluster taxonomy
cluster_taxonomy_mg <- read.delim("../../iba_data/cluster_taxonomy_mg.tsv")

# Add life history data
lht <- read.delim("../../traits/clade_trait_data_mg.tsv")
cluster_taxonomy_mg$Niche <- lht$Niche[match(cluster_taxonomy_mg$Clade,lht$Clade)]
cluster_taxonomy_mg$Habitat <- lht$Habitat[match(cluster_taxonomy_mg$Clade,lht$Clade)]

# Generate niche subsets
niches <- c("Saprophage","Saprophage-parasitoid","Phytophage","Phytophage-parasitoid","Predator","Predator-parasitoid")
sp_matrix <- list()
for (i in 1:length(niches)) {
    niche <- niches[i]
    clusters <- cluster_taxonomy_mg$cluster[cluster_taxonomy_mg$Niche==niche]
    sp_matrix[[i]] <- sp_matrix_mg[,colnames(sp_matrix_mg) %in% clusters]
}


# Generate total species accumulation curves
# ------------------------------------------

cat("Generating plot data\n")
spDF <- list()
poolDF <- list()
for (i in 1:length(niches)) {
    X <- specaccum(sp_matrix[[i]], method="random", permutations=100)
    spDF[[i]] <- data.frame(species=X$richness, sd=X$sd, sites=X$sites)
    poolDF[[i]] <- data.frame(poolaccum(sp_matrix[[i]], permutations=100)$means)
}

cat("Rendering plot data\n")

# Plot fxn
accum_plot <- function(spDF, poolDF, niche) {
    ggplot(data=spDF, aes(x=sites, y=species)) +
        geom_line(linewidth=lsize) +
        geom_ribbon(aes(ymin=species-sd, ymax=species+sd), alpha=0.2, fill="blue") +
        geom_line(data=poolDF, linewidth=lsize, aes(x=N, y=Chao), linetype="dashed") +
        theme_linedraw(base_size=bs) +
        labs(x="Number of sites", y="Number of OTUs", title=niche)
}

# Make plots
p <- list()
for (i in 1:length(niches)) {
    p[[i]] <- accum_plot(spDF[[i]], poolDF[[i]], niches[i])
}

# Save figure
ggsave( file = "../figs/Fig_accumulation_mg_niches.jpg",
        width = 14,
        height = 21,
        plot = p[[1]] + p[[2]] + p[[3]] + p[[4]] + p[[5]] + p[[6]] +
            plot_layout(axis_titles="collect", ncol=2) +
            plot_annotation(tag_levels="A")
        )

