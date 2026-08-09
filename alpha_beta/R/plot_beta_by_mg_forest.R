# Load required libraries
library(ggplot2)
library(tidyverse)
library(patchwork)


# Set random seed
set.seed(10)

# Read in functions and fig settings
source("functions.R")
source("../../fig_settings/fig_colours.R")


# Read in and format data
# -----------------------

# Read species * site matrices
sp_matrix_mg <- readRDS("../../iba_data/site_otu_occurrence_combined_mg.rds")
sp_matrix_mg <- 1*sp_matrix_mg  # Convert to numeric

# Read in sample metadata
site_meta_mg <- read.delim("../../iba_data/malaise_litter_sample_meta_mg.tsv")

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
sf_meta_df <- site_meta_df |> st_as_sf(coords=c("longitude", "latitude"), crs = 4326)
sf_meta_rf <- site_meta_rf |> st_as_sf(coords=c("longitude", "latitude"), crs = 4326)
dist_df <- get_dists(sf_meta_df)
dist_rf <- get_dists(sf_meta_rf)

# Compute beta diversity metric
betapart_df  <- partition_beta_diversity(sp_matrix_total_df, trap_dist=dist_df)
betapart_rf  <- partition_beta_diversity(sp_matrix_total_rf, trap_dist=dist_rf)


# Get monotonic gam fits
# ----------------------

# Fits a single cubic regression spline to the distance / turnover data
mono_df <- monotonic_gam(betapart_df, nK=5)
mono_rf <- monotonic_gam(betapart_rf, nK=5)


# Render plots
# ------------

p1 <- ggplot(betapart_rf, aes(distance, jaccard)) +
    geom_point(alpha=0.6, size=2, colour=rf_col, show.legend=FALSE) +
    theme_linedraw(base_size=20) +
    theme(plot.title=element_text(hjust=0.5)) +
    geom_line(data=mono_rf, aes(distance, pred_fit), lwd=2) +
#    geom_smooth(aes(distance, jaccard), method="loess") +
    scale_y_continuous(limits=c(0.4, 1)) +
    scale_x_continuous(limits=c(0,1510)) +
    labs(title="Rainforest", x="Distance (km)", y="Dissimilarity (J)") 

p2 <- ggplot(betapart_df, aes(distance , jaccard)) +
    geom_point(alpha=0.6, size=2, colour=df_col, show.legend=FALSE) +
    theme_linedraw(base_size=20) +
    theme(plot.title=element_text(hjust=0.5)) +
    geom_line(data=mono_df, aes(distance, pred_fit), lwd=2) +
#    geom_smooth(aes(distance, jaccard), method="loess") +
    scale_y_continuous(limits=c(0.4, 1)) +
    scale_x_continuous(limits=c(0,1510)) +
    labs(title="Dry forest", x="Distance (km)", y=NULL) 


# Generate temporal beta plots
# ----------------------------

cat("Generating temporal beta plots\n")


# Read in and format data
# -----------------------

# Read species * sample matrix
sp_sample_matrix_mg <- readRDS("../../iba_data/sample_otu_abundance_malaise_mg.rds")
sp_sample_matrix_mg <- 1*(sp_sample_matrix_mg > 0)  # Convert to 0/1 occurrence

# Read in sample metadata (for mid_date field compatibility with get_temp_dists)
sample_meta_mg <- read.delim("../../iba_data/malaise_sample_meta_mg.tsv")

# Split into forest habitat subsets (above site_meta matrices have info for all samples also)
sp_sample_matrix_df <- sp_sample_matrix_mg[rownames(sp_sample_matrix_mg) %in% site_meta_df$sampleID_NGI,]
sp_sample_matrix_rf <- sp_sample_matrix_mg[rownames(sp_sample_matrix_mg) %in% site_meta_rf$sampleID_NGI,]

# Filter out singletons
sp_sample_matrix_df <- sp_sample_matrix_df[,colSums(sp_sample_matrix_df) > 1]
sp_sample_matrix_rf <- sp_sample_matrix_rf[,colSums(sp_sample_matrix_rf) > 1]


# Temporal turnover analysis
# --------------------------

# Define function to compute the within-site data
compute_temp_data <- function(D, M) {
    res <- data.frame()
    for (trap in unique(M$trapID)) {
        trap_samples <- M$sampleID_NGI[M$trapID==trap]
        E <- D[row.names(D) %in% trap_samples,]
        E <- E[rowSums(E)>0,colSums(E)>1]
        if (ncol(E)<=2)
            next
        X <- M[match(row.names(E),M$sampleID_NGI),]
        dist <- get_temp_dists(X)
        betapart <- partition_beta_temp_diversity(E, sample_dist=dist)
        res <- rbind(res, betapart)
    }
    return(res)
}

# Compute temporal within-site beta diversity metric
cat("Computing temporal data for MG dry forest sites\n")
betapart_temp_df <- compute_temp_data(sp_sample_matrix_df, sample_meta_mg)
cat("Computing temporal data for MG rainforest sites\n")
betapart_temp_rf <- compute_temp_data(sp_sample_matrix_rf, sample_meta_mg)


# Get monotonic gam fits
# ----------------------

# Fits a single cubic regression spline to the distance / turnover data
mono_temp_df <- monotonic_gam(betapart_temp_df, nK=5)
mono_temp_rf <- monotonic_gam(betapart_temp_rf, nK=5)


# Render plots
# ------------

p3 <- ggplot(betapart_temp_rf, aes(distance, jaccard)) +
    geom_point(alpha=0.02, size=2, colour=rf_col, show.legend=FALSE) +
    theme_linedraw(base_size=20) +
    geom_line(data=mono_temp_rf, aes(distance, pred_fit), lwd=2) +
    scale_y_continuous(limits=c(0.4, 1)) +
    scale_x_continuous(limits=c(0,185)) +
    labs(title=NULL,x="Days", y="Dissimilarity (J)")

p4 <- ggplot(betapart_temp_df, aes(distance, jaccard)) +
    geom_point(alpha=0.02, size=2, colour=df_col, show.legend=FALSE) +
    theme_linedraw(base_size=20) +
    geom_line(data=mono_temp_df, aes(distance, pred_fit), lwd=2) +
    scale_y_continuous(limits=c(0.4, 1)) +
    scale_x_continuous(limits=c(0,185)) +
    labs(title=NULL,x="Days", y=NULL)


# Save plots
ggsave( file = "../figs/Fig_beta_by_mg_forest.jpg",
        width = 14.0,
        height = 14.0,
        plot = p1 + p2 + p3 + p4 +
            plot_annotation(tag_levels="A") +
            plot_layout(ncol=2, axis_titles="collect_x")
      )

