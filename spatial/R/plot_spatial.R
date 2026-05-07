# Load required libraries
library(ggplot2)
library(tidyverse)
library(patchwork)


# Set random seed
set.seed(10)
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

# Filter out singletons
sp_matrix_total_se <- sp_matrix_se[,colSums(sp_matrix_se) > 1]
sp_matrix_total_mg <- sp_matrix_mg[,colSums(sp_matrix_mg) > 1]

# Generate MG trap forest habitat subsets
site_meta_df <- site_meta_mg[site_meta_mg$trap_habitat=="Dry_Forest",]
site_meta_rf <- site_meta_mg[site_meta_mg$trap_habitat %in% c("Montane_Rainforest","Rainforest"),]

sp_matrix_df <- sp_matrix_mg[rownames(sp_matrix_mg) %in% unique(site_meta_df$trapID),]
sp_matrix_rf <- sp_matrix_mg[rownames(sp_matrix_mg) %in% unique(site_meta_rf$trapID),]

# Filter out singletons
sp_matrix_total_df <- sp_matrix_df[,colSums(sp_matrix_df) > 1]
sp_matrix_total_rf <- sp_matrix_rf[,colSums(sp_matrix_rf) > 1]


# Turnover analysis
# -----------------

# Get distances between traps
sf_meta_mg <- site_meta_mg |> st_as_sf(coords=c("longitude", "latitude"), crs = 4326) 
sf_meta_se <- site_meta_se |> st_as_sf(coords=c("longitude", "latitude"), crs = 4326) 
dist_swe <- get_dists(sf_meta_se)
dist_mad <- get_dists(sf_meta_mg)

sf_meta_df <- site_meta_df |> st_as_sf(coords=c("longitude", "latitude"), crs = 4326)
sf_meta_rf <- site_meta_rf |> st_as_sf(coords=c("longitude", "latitude"), crs = 4326)
dist_df <- get_dists(sf_meta_df)
dist_rf <- get_dists(sf_meta_rf)

# Compute beta diversity metric
betapart_mad <- partition_beta_diversity(sp_matrix_total_mg, trap_dist=dist_mad)
betapart_swe <- partition_beta_diversity(sp_matrix_total_se, trap_dist=dist_swe)
betapart_df  <- partition_beta_diversity(sp_matrix_total_df, trap_dist=dist_df)
betapart_rf  <- partition_beta_diversity(sp_matrix_total_rf, trap_dist=dist_rf)


# Get monotonic gam fits
# ----------------------

# Fits a single cubic regression spline to the distance / turnover data
mono_mg <- monotonic_gam(betapart_mad, nK=5)
mono_se <- monotonic_gam(betapart_swe, nK=5)
mono_df <- monotonic_gam(betapart_df, nK=5)
mono_rf <- monotonic_gam(betapart_rf, nK=5)



# Render plots
# ------------

p1 <- ggplot(betapart_mad, aes(distance , jaccard)) +
    geom_point(alpha=0.2, size=2, aes(colour=distance), show.legend=FALSE) +
    theme_linedraw(base_size=20) +
    geom_line(data=mono_mg, aes(distance, pred_fit), lwd=2) +
#    geom_smooth(aes(distance, jaccard), method="loess") +
    scale_colour_viridis_c(option="rocket") +
    scale_y_continuous(limits=c(0.4, 1)) +
    labs(x="Distance (km)", y="Dissimilarity (J)", colour="Distance (km)") 

p2 <- ggplot(betapart_swe, aes(distance, jaccard)) +
    geom_point(alpha=0.1, size=2, aes(colour=distance), show.legend=FALSE) +
    theme_linedraw(base_size=20) +
    geom_line(data=mono_se, aes(distance, pred_fit), lwd=2) +
#    geom_smooth(aes(distance, jaccard), method="loess") +
    scale_colour_viridis_c(option="rocket") +
    scale_y_continuous(limits=c(0.4, 1)) +
    labs(x="Distance (km)", y="Dissimilarity (J)", colour="Distance (km)") 

p3 <- ggplot(betapart_rf, aes(distance, jaccard)) +
    geom_point(alpha=0.1, size=2, aes(colour=distance), show.legend=FALSE) +
    theme_linedraw(base_size=20) +
    geom_line(data=mono_rf, aes(distance, pred_fit), lwd=2) +
#    geom_smooth(aes(distance, jaccard), method="loess") +
    scale_colour_viridis_c(option="rocket") +
    scale_y_continuous(limits=c(0.4, 1)) +
    labs(x="Distance (km)", y="Dissimilarity (J)", colour="Distance (km)") 

p4 <- ggplot(betapart_df, aes(distance , jaccard)) +
    geom_point(alpha=0.2, size=2, aes(colour=distance), show.legend=FALSE) +
    theme_linedraw(base_size=20) +
    geom_line(data=mono_df, aes(distance, pred_fit), lwd=2) +
#    geom_smooth(aes(distance, jaccard), method="loess") +
    scale_colour_viridis_c(option="rocket") +
    scale_y_continuous(limits=c(0.4, 1)) +
    labs(x="Distance (km)", y="Dissimilarity (J)", colour="Distance (km)") 

# Save plots
ggsave( file = "../figs/Fig_spatial.jpg",
        width = 14.0,
        height = 14.0,
        plot = p1 + p2 + p3 + p4 +
            plot_annotation(tag_levels="A") +
            plot_layout(ncol=2, axis_titles="collect", guides="collect")
      )

