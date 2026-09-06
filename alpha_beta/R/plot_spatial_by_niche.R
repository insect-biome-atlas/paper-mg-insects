# Script for examining spatial turnover by niche
# There are too few data points to make meaningful
# turnover analysis for predator parasitoids and
# and saprophage parasitoids in MG. Also note
# that temporal turnover by niche is challenging
# to analyze as the spread becomes difficult to
# handle when the data are analyzed per sample,
# so we refrain from these analyses.

# Load required libraries
library(ggplot2)
library(tidyverse)
library(patchwork)
source("functions.R")
source("../../fig_settings/fig_colours.R")

remove_singletons <- TRUE 
lsize <- 1

# Read in and format data
# -----------------------

# Read site * species matrices
# These are already filtered so that SE data have a minimum requirement of 10 samples per site
# and SE data only comprise forest traps
site_otu_matrix_mg <- readRDS("../../iba_data/site_otu_occurrence_combined_mg.rds")
site_otu_matrix_se <- readRDS("../../iba_data/site_otu_occurrence_combined_se.rds")
site_otu_matrix_mg <- 1*site_otu_matrix_mg  # Convert to numeric
site_otu_matrix_se <- 1*site_otu_matrix_se

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
if (remove_singletons) {
    site_otu_matrix_mg <- site_otu_matrix_mg[,colSums(site_otu_matrix_mg) > 1]
    site_otu_matrix_se <- site_otu_matrix_se[,colSums(site_otu_matrix_se) > 1]
}

# Get distances between traps
sf_meta_mg <- site_meta_mg |> st_as_sf(coords=c("longitude", "latitude"), crs = 4326) 
sf_meta_se <- site_meta_se |> st_as_sf(coords=c("longitude", "latitude"), crs = 4326) 
dist_mg <- get_dists(sf_meta_mg)
dist_se <- get_dists(sf_meta_se)


# Render plots
# ------------

spatial_niche_plot <- function(spmat, dist, niche, tax_niche, colour, alpha=0.2) {

    spmat <- spmat[,colnames(spmat) %in% tax_niche$cluster[tax_niche$Niche==niche]]
    if (remove_singletons)
        spmat <- spmat[ ,colSums(spmat) > 1]
    betapart <- partition_beta_diversity(spmat, dist)
    mono <- monotonic_gam(betapart, nK=5)

    ggplot(betapart, aes(distance , jaccard)) +
        geom_point(alpha=alpha, size=2, colour=colour, show.legend=FALSE) +
        theme_linedraw(base_size=20) +
        geom_line(data=mono, aes(distance, pred_fit), lwd=lsize) +
        scale_y_continuous(limits=c(0.4, 1)) +
        scale_x_continuous(limits=c(0,1510)) +
        labs(title=niche, x="Distance (km)", y="Dissimilarity (J)")
}

p1 <- spatial_niche_plot(site_otu_matrix_mg, dist_mg, "Saprophage", tax_niche_mg, mg_col)
p2 <- spatial_niche_plot(site_otu_matrix_mg, dist_mg, "Phytophage", tax_niche_mg, mg_col)
p3 <- spatial_niche_plot(site_otu_matrix_mg, dist_mg, "Predator", tax_niche_mg, mg_col)
p4 <- spatial_niche_plot(site_otu_matrix_mg, dist_mg, "Phytophage-parasitoid", tax_niche_mg, mg_col)

q1 <- spatial_niche_plot(site_otu_matrix_se, dist_se, "Saprophage", tax_niche_se, se_col, 0.1)
q2 <- spatial_niche_plot(site_otu_matrix_se, dist_se, "Phytophage", tax_niche_se, se_col, 0.1)
q3 <- spatial_niche_plot(site_otu_matrix_se, dist_se, "Predator", tax_niche_se, se_col, 0.1)
q4 <- spatial_niche_plot(site_otu_matrix_se, dist_se, "Saprophage-parasitoid", tax_niche_se, se_col, 0.1)
q5 <- spatial_niche_plot(site_otu_matrix_se, dist_se, "Phytophage-parasitoid", tax_niche_se, se_col, 0.1)
q6 <- spatial_niche_plot(site_otu_matrix_se, dist_se, "Predator-parasitoid", tax_niche_se, se_col, 0.1)

# Save plots
ggsave( file = "../figs/Fig_spatial_niche_mg.jpg",
        width = 16.0,
        height = 14.0,
        plot = p1 + p2 + p3 + p4 +
            plot_annotation(tag_levels="A") +
            plot_layout(ncol=2, axis_titles="collect")
      )

ggsave( file = "../figs/Fig_spatial_niche_se.jpg",
        width = 23.0,
        height = 14.0,
        plot = q1 + q2 + q3 + q4 + q5 + q6 +
            plot_annotation(tag_levels="A") +
            plot_layout(ncol=3, axis_titles="collect_x")
      )

