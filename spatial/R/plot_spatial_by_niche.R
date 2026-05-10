# Load required libraries
library(ggplot2)
library(tidyverse)
library(patchwork)
source("functions.R")

# Read in and format data
# -----------------------

# Read species * site matrices
# These are already filtered so that SE data have a minimum requirement of 10 samples per site
# and SE data only comprise forest traps
sp_matrix_mg <- readRDS("../../iba_data/site_otu_occurrence_combined_mg.rds")
sp_matrix_se <- readRDS("../../iba_data/site_otu_occurrence_combined_se.rds")
sp_matrix_mg <- 1*sp_matrix_mg  # Convert to numeric
sp_matrix_se <- 1*sp_matrix_se

# Read in sample metadata
site_meta_mg <- read.delim("../../iba_data/malaise_litter_sample_meta_mg.tsv")
site_meta_se <- read.delim("../../iba_data/malaise_litter_sample_meta_se.tsv")

# Read in life history data
lht_mg <- read.delim("../../traits/clade_trait_data_mg.tsv")
lht_se <- read.delim("../../traits/clade_trait_data_se.tsv")

# Read in taxonomy data
tax_niche_mg <- readRDS("../../iba_data/cluster_taxonomy_mg.rds")
tax_niche_se <- readRDS("../../iba_data/cluster_taxonomy_se.rds")

# Add Niche information to the taxonomy
tax_niche_mg$Niche <- lht_mg$Niche[match(tax_niche_mg$Clade,lht_mg$Clade)]
tax_niche_se$Niche <- lht_se$Niche[match(tax_niche_se$Clade,lht_se$Clade)]

# Filter out singletons
sp_matrix_total_se <- sp_matrix_se[,colSums(sp_matrix_se) > 1]
sp_matrix_total_mg <- sp_matrix_mg[,colSums(sp_matrix_mg) > 1]

# Get distances between traps
sf_meta_mg <- site_meta_mg |> st_as_sf(coords=c("longitude", "latitude"), crs = 4326) 
sf_meta_se <- site_meta_se |> st_as_sf(coords=c("longitude", "latitude"), crs = 4326) 
dist_mg <- get_dists(sf_meta_mg)
dist_se <- get_dists(sf_meta_se)


# Render plots
# ------------

spatial_niche_plot <- function(spmat, dist, niche, tax_niche, show_legend) {

    spmat <- spmat[,colnames(spmat) %in% tax_niche$cluster[tax_niche$Niche==niche]]
    betapart <- partition_beta_diversity(spmat, dist)
    mono <- monotonic_gam(betapart, nK=5)

    ggplot(betapart, aes(distance , jaccard)) +
        geom_point(alpha=0.2, size=2, aes(colour=distance), show.legend=show_legend) +
        theme_linedraw(base_size=20) +
        geom_line(data=mono, aes(distance, pred_fit), lwd=2) +
        scale_colour_viridis_c(option="rocket") +
        scale_y_continuous(limits=c(0.4, 1)) +
        scale_x_continuous(limits=c(0,1510)) +
        labs(title=niche, x="Distance (km)", y="Dissimilarity (J)", colour="Distance (km)")
}

p1 <- spatial_niche_plot(sp_matrix_total_mg, dist_mg, "Saprophage", tax_niche_mg, FALSE)
p2 <- spatial_niche_plot(sp_matrix_total_mg, dist_mg, "Phytophage", tax_niche_mg, FALSE)
p3 <- spatial_niche_plot(sp_matrix_total_mg, dist_mg, "Predator", tax_niche_mg, TRUE)
p4 <- spatial_niche_plot(sp_matrix_total_mg, dist_mg, "Saprophage-parasitoid", tax_niche_mg, FALSE)
p5 <- spatial_niche_plot(sp_matrix_total_mg, dist_mg, "Phytophage-parasitoid", tax_niche_mg, FALSE)
p6 <- spatial_niche_plot(sp_matrix_total_mg, dist_mg, "Predator-parasitoid", tax_niche_mg, FALSE)

q1 <- spatial_niche_plot(sp_matrix_total_se, dist_se, "Saprophage", tax_niche_se, FALSE)
q2 <- spatial_niche_plot(sp_matrix_total_se, dist_se, "Phytophage", tax_niche_se, FALSE)
q3 <- spatial_niche_plot(sp_matrix_total_se, dist_se, "Predator", tax_niche_se, TRUE)
q4 <- spatial_niche_plot(sp_matrix_total_se, dist_se, "Saprophage-parasitoid", tax_niche_se, FALSE)
q5 <- spatial_niche_plot(sp_matrix_total_se, dist_se, "Phytophage-parasitoid", tax_niche_se, FALSE)
q6 <- spatial_niche_plot(sp_matrix_total_se, dist_se, "Predator-parasitoid", tax_niche_se, FALSE)

# Save plots
ggsave( file = "../figs/Fig_spatial_niche_mg.jpg",
        width = 23.0,
        height = 14.0,
        plot = p1 + p2 + p3 + p4 + p5 + p6 +
            plot_annotation(tag_levels="A") +
            plot_layout(ncol=3, axis_titles="collect")
      )

ggsave( file = "../figs/Fig_spatial_niche_se.jpg",
        width = 23.0,
        height = 14.0,
        plot = q1 + q2 + q3 + q4 + q5 + q6 +
            plot_annotation(tag_levels="A") +
            plot_layout(ncol=3, axis_titles="collect_x")
      )

ggsave( file = "../figs/Fig_spatial_phytopage_community_mg.jpg",
        width = 16.0,
        height = 7.0,
        plot = p2 + p5 +
            plot_annotation(tag_levels="A") +
            plot_layout(ncol=2, axis_titles="collect")
      )

