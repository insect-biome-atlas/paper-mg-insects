# Plot site specpool accumulation for taxonomic groups

library(vegan)
library(ggplot2)
library(patchwork)


# Read in and tidy data
# ---------------------

cat("Reading and formating data\n")

# Read species * site matrix
sp_matrix_se <- readRDS("../../iba_data/site_otu_occurrence_combined_complete_se.rds")

# Read in sample metadata
site_meta_se <- read.delim("../../iba_data/malaise_litter_sample_meta_se.tsv")

# Read in cluster taxonomy
cluster_taxonomy_se <- readRDS("../../iba_data/cluster_taxonomy_se.rds")

# Generate taxonomy subsets
niches <- c("Saprophage","Saprophage-parasitoid","Phytophage","Phytophage-parasitoid","Predator","Predator-parasitoid")
big_five <- c("Diptera","Hymenoptera","Coleoptera","Lepidoptera","Hemiptera")
sp_matrix <- list()
for (i in 1:length(big_five)) {
    order <- big_five[i]
    clusters <- cluster_taxonomy_se$cluster[cluster_taxonomy_se$Order==order]
    sp_matrix[[i]] <- sp_matrix_se[,colnames(sp_matrix_se) %in% clusters]
}
clusters <- cluster_taxonomy_se$cluster[!(cluster_taxonomy_se$Order %in% big_five)]
sp_matrix[[6]] <- sp_matrix_se[,colnames(sp_matrix_se) %in% clusters]

# Get total species number estimates
taxon_comp <- read.delim("../data/taxon_comp.tsv")
spp_estimate <- function(taxon) { taxon_comp$species[taxon_comp$order==taxon] }
dipt_estimate <- spp_estimate("Diptera")
hyme_estimate <- spp_estimate("Hymenoptera")
cole_estimate <- spp_estimate("Coleoptera")
lepi_estimate <- spp_estimate("Lepidoptera")
hemi_estimate <- spp_estimate("Hemiptera")
other_estimate <- spp_estimate("Other")


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
plot_accum <- function(spDF, poolDF, spp_estimate, plot_title) {
    ggplot(data=spDF, aes(x=sites, y=species)) +
        ggtitle(plot_title) +
        geom_line(linewidth=2) +
        geom_ribbon(aes(ymin=species-sd, ymax=species+sd), alpha=0.2, fill="blue") +
        geom_line(data=poolDF, linewidth=2, aes(x=N, y=Chao), linetype="dashed") +
        geom_hline(yintercept=spp_estimate, colour="red", linetype="dotted", linewidth=1) +
        theme_linedraw(base_size=20) +
        labs(x="Number of sites", y="Number of OTUs")
}

p1 <- plot_accum(spDF[[1]], poolDF[[1]], dipt_estimate, "Diptera")
p2 <- plot_accum(spDF[[2]], poolDF[[2]], hyme_estimate, "Hymenoptera")
p3 <- plot_accum(spDF[[3]], poolDF[[3]], cole_estimate, "Coleoptera")
p4 <- plot_accum(spDF[[4]], poolDF[[4]], lepi_estimate, "Lepidoptera")
p5 <- plot_accum(spDF[[5]], poolDF[[5]], hemi_estimate, "Hemiptera")
p6 <- plot_accum(spDF[[6]], poolDF[[6]], other_estimate, "Other")

# Save plots
ggsave( "../figs/Fig_site_pool_acc_taxonomic_complete_se.jpg",
        width=21,
        height=14,
        plot=p1 + p2 + p3 + p4 + p5 + p6 +
            plot_layout(axis_titles="collect", ncol=3) +
            plot_annotation(tag_levels="A")
      )

