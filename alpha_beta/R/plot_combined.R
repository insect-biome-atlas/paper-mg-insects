# Generate accumulation and beta plots

library(vegan)
library(ggplot2)
library(patchwork)
library(tidyverse)

source("functions.R")
source("../../fig_settings/fig_colours.R")

# Set basic plot params and analysis params
bs    <- 20
lsize1 <- 1.0
lsize2 <- 2.0
set.seed(10)
remove_singletons <- TRUE 

# Read data
# ---------

# Read OTU site occurrence tables
cat("Reading data\n")
site_otu_occurrence_combined_mg <- readRDS("../../iba_data/site_otu_occurrence_combined_mg.rds")
site_otu_occurrence_combined_se <- readRDS("../../iba_data/site_otu_occurrence_combined_se.rds")


# Generate total species accumulation curves
# ------------------------------------------

cat("Generating accumulation plot data\n")
cat("Generating specaccum permutations\n")
spMG <- specaccum(site_otu_occurrence_combined_mg, method="random", permutations=100)
spSE <- specaccum(site_otu_occurrence_combined_se, method="random", permutations=100)
cat("Generating poolaccum permutations\n")
poolMG <- poolaccum(site_otu_occurrence_combined_mg, permutations=100)
poolSE <- poolaccum(site_otu_occurrence_combined_se, permutations=100)

cat("Formatting accumulation plot data\n")
mgDF <- data.frame(species=spMG$richness, sd=spMG$sd, sites=spMG$sites)
seDF <- data.frame(species=spSE$richness, sd=spSE$sd, sites=spSE$sites)
mgPoolDF <- data.frame(poolMG$means)
sePoolDF <- data.frame(poolSE$means)
mgPoolDF <- mgPoolDF[mgPoolDF$N >=5,]

# Madagascar
p1 <- ggplot(data=mgDF, aes(x=sites, y=species)) +
        geom_line(linewidth=1.0) +
        geom_ribbon(aes(ymin=species-sd, ymax=species+sd), alpha=0.2, fill=mg_col) +
        geom_line(data=mgPoolDF, linewidth=lsize1, aes(x=N, y=Chao), linetype="dashed") +
        theme_linedraw(base_size=bs) +
        theme(plot.title=element_text(hjust=0.5)) +
        scale_y_continuous(limits=c(0,126000),
                           breaks=c(25000,50000,75000,100000,125000),
                           labels=c("25k","50k","75k","100k","125k")) +
        labs(title="Madagascar", x="Number of sites", y="Number of OTUs")

# Sweden
p2 <- ggplot(data=seDF, aes(x=sites, y=species)) +
        geom_line(linewidth=1.0) +
        geom_ribbon(aes(ymin=species-sd, ymax=species+sd) , alpha=0.2, fill=se_col) +
        geom_line(data=sePoolDF, linewidth=lsize2, aes(x=N, y=Chao), linetype="dashed") +
        theme_linedraw(base_size=bs) +
        theme(plot.title=element_text(hjust=0.5)) +
        scale_y_continuous(limits=c(0,126000),
                           breaks=c(25000,50000,75000,100000,125000),
                           labels=c("25k","50k","75k","100k","125k")) +
        labs(title="Sweden", x="Number of sites", y=NULL)


# Generate spatial beta plots
# ---------------------------

cat("Generating spatial beta plots\n")


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

# Filter out singletons or not
if (remove_singletons) {
    site_otu_matrix_mg <- site_otu_matrix_mg[,colSums(site_otu_matrix_mg) > 1]
    site_otu_matrix_se <- site_otu_matrix_se[,colSums(site_otu_matrix_se) > 1]
}


# Spatial turnover analysis
# -------------------------

# Get distances between traps
sf_meta_mg <- site_meta_mg |> st_as_sf(coords=c("longitude", "latitude"), crs = 4326) 
sf_meta_se <- site_meta_se |> st_as_sf(coords=c("longitude", "latitude"), crs = 4326) 
dist_swe <- get_dists(sf_meta_se)
dist_mad <- get_dists(sf_meta_mg)

# Compute spatial beta diversity metric
betapart_mad <- partition_beta_diversity(site_otu_matrix_mg, trap_dist=dist_mad)
betapart_swe <- partition_beta_diversity(site_otu_matrix_se, trap_dist=dist_swe)


# Get monotonic gam fits
# ----------------------

# Fits a single cubic regression spline to the distance / turnover data
mono_mg <- monotonic_gam(betapart_mad, nK=5)
mono_se <- monotonic_gam(betapart_swe, nK=5)


# Render plots
# ------------

p3 <- ggplot(betapart_mad, aes(distance, jaccard)) +
    geom_point(alpha=0.2, size=2, colour=mg_col, show.legend=FALSE) +
    theme_linedraw(base_size=20) +
    geom_line(data=mono_mg, aes(distance, pred_fit), lwd=lsize2) +
    scale_y_continuous(limits=c(0.4, 1)) +
    scale_x_continuous(limits=c(0,1510)) +
    labs(title=NULL,x="Distance (km)", y="Spatial turnover (Jaccard dissim.)") 

p4 <- ggplot(betapart_swe, aes(distance, jaccard)) +
    geom_point(alpha=0.1, size=2, colour=se_col, show.legend=FALSE) +
    theme_linedraw(base_size=20) +
    geom_line(data=mono_se, aes(distance, pred_fit), lwd=lsize2) +
    scale_y_continuous(limits=c(0.4, 1)) +
    scale_x_continuous(limits=c(0,1510)) +
    labs(title=NULL,x="Distance (km)", y=NULL) 


# Generate temporal beta plots
# ----------------------------

cat("Generating temporal beta plots\n")


# Read in and format data
# -----------------------

# Read sample * species matrices
sample_otu_matrix_mg <- readRDS("../../iba_data/sample_otu_abundance_malaise_mg.rds")
sample_otu_matrix_se <- readRDS("../../iba_data/sample_otu_abundance_malaise_se.rds")
sample_otu_matrix_mg <- 1*(sample_otu_matrix_mg > 0)  # Convert to 0/1 occurrence
sample_otu_matrix_se <- 1*(sample_otu_matrix_se > 0)

# Read in sample metadata
sample_meta_mg <- read.delim("../../iba_data/malaise_sample_meta_mg.tsv")
sample_meta_se <- read.delim("../../iba_data/malaise_sample_meta_se.tsv")

# Correc trapID error discovered late in the process (time diff==0)
sample_meta_mg$trapID[sample_meta_mg$sampleID_FIELD=="S3UVGT"]<-"TQKQRH"


# Temporal turnover analysis
# --------------------------

# Define function to compute the within-site data
compute_temp_data <- function(D, M) {
    res <- data.frame()
    for (trap in unique(M$trapID)) {
        trap_samples <- M$sampleID_NGI[M$trapID==trap]
        E <- D[row.names(D) %in% trap_samples,]
        E <- E[rowSums(E)>0, colSums(E)>0]  # Remove empty rows and columns
        if (remove_singletons) 
            E <- E[, colSums(E) > 1]
        if (nrow(E)<2 || ncol(E)<2) # Sanity check
            next
        X <- M[match(row.names(E),M$sampleID_NGI),]
        dist <- get_temp_dists(X)
        betapart <- partition_beta_temp_diversity(E, sample_dist=dist)
        res <- rbind(res, betapart)
    }
    return(res)
}

# Compute temporal within-site beta diversity metric
cat("Computing temporal data for MG\n")
betapart_temp_mad <- compute_temp_data(sample_otu_matrix_mg, sample_meta_mg)
cat("Computing temporal data for SE\n")
betapart_temp_swe <- compute_temp_data(sample_otu_matrix_se, sample_meta_se)

# Get monotonic gam fits
# ----------------------

# Fits a single cubic regression spline to the distance / turnover data
mono_temp_mg <- monotonic_gam(betapart_temp_mad, nK=5)
mono_temp_se <- monotonic_gam(betapart_temp_swe, nK=5)


# Render plots
# ------------

p5 <- ggplot(betapart_temp_mad, aes(distance, jaccard)) +
    geom_point(alpha=0.01, size=2, colour=mg_col, show.legend=FALSE) +
    theme_linedraw(base_size=20) +
    geom_line(data=mono_temp_mg, aes(distance, pred_fit), lwd=lsize2) +
    scale_y_continuous(limits=c(0.4, 1)) +
    scale_x_continuous(limits=c(0,185)) +
    labs(title=NULL,x="Days", y="Temporal turnover (Jaccard dissim.)")

p6 <- ggplot(betapart_temp_swe, aes(distance, jaccard)) +
    geom_point(alpha=0.005, size=2, colour=se_col, show.legend=FALSE) +
    theme_linedraw(base_size=20) +
    geom_line(data=mono_temp_se, aes(distance, pred_fit), lwd=lsize2) +
    scale_y_continuous(limits=c(0.4, 1)) +
    scale_x_continuous(limits=c(0,185)) +
    labs(title=NULL,x="Days", y=NULL)


# Save plots
ggsave( file = "../figs/Fig_alpha_beta.jpg",
        width = 16.0,
        height = 21.0,
        plot = p1 + p2 + p3 + p4 + p5 + p6 +
            plot_annotation(tag_levels="A") +
            plot_layout(ncol=2, axis_titles="collect_x")
      )

