# Load required libraries
library(ggplot2)
library(tidyverse)
library(patchwork)
source("functions.R")
source("../../fig_settings/fig_colours.R")

lsize <- 1
remove_singletons <- TRUE


# Read in and format data
# -----------------------

# Read sample * species matrices
sample_otu_matrix_mg <- readRDS("../../iba_data/sample_otu_abundance_malaise_mg.rds")
sample_otu_matrix_se <- readRDS("../../iba_data/sample_otu_abundance_malaise_se.rds")
sample_otu_matrix_mg <- 1*(sample_otu_matrix_mg > 0)  # Convert to 0/1 occurrence
sample_otu_matrix_se <- 1*(sample_otu_matrix_se > 0)

# Read in sample metadata (only MT samples have time diffs)
sample_meta_mg <- read.delim("../../iba_data/malaise_sample_meta_mg.tsv")
sample_meta_se <- read.delim("../../iba_data/malaise_sample_meta_se.tsv")

# Correc trapID error discovered late in the process (time diff==0)
sample_meta_mg$trapID[sample_meta_mg$sampleID_FIELD=="S3UVGT"]<-"TQKQRH"

# Read in life history data
lht_mg <- read.delim("../../traits/clade_trait_data_mg.tsv")
lht_se <- read.delim("../../traits/clade_trait_data_se.tsv")

# Read in taxonomy data
tax_niche_mg <- readRDS("../../iba_data/cluster_taxonomy_mg.rds")
tax_niche_se <- readRDS("../../iba_data/cluster_taxonomy_se.rds")

# Add Niche information to the taxonomy
tax_niche_mg$Niche <- lht_mg$Niche[match(tax_niche_mg$Clade,lht_mg$Clade)]
tax_niche_se$Niche <- lht_se$Niche[match(tax_niche_se$Clade,lht_se$Clade)]


# Render plots
# ------------

# Define function to compute the within-site data
compute_temp_data <- function(D, M) {
    res <- data.frame()
    for (trap in unique(M$trapID)) {
        trap_samples <- M$sampleID_NGI[M$trapID==trap]
        E <- D[row.names(D) %in% trap_samples,]
        E <- E[rowSums(E)>0,colSums(E)>0]   # Remove empty rows and columns
        if (remove_singletons)
            E <- E[ ,colSums(E)>1]   # Only keep species occuring in at least 2 samples
        if (nrow(E)<2 || ncol(E)<2) # Sanity check 
            next
        X <- M[match(row.names(E),M$sampleID_NGI),]
        dist <- get_temp_dists(X)
        betapart <- partition_beta_temp_diversity(E, sample_dist=dist)
        res <- rbind(res, betapart)
    }
    return(res[!is.na(res$jaccard),])
}


# Define function to render the plots
temporal_niche_plot <- function(sample_otu_mat, sample_meta, niche, tax_niche, colour) {

    sample_otu_mat <- sample_otu_mat[,colnames(sample_otu_mat) %in% tax_niche$cluster[tax_niche$Niche==niche]]
    sample_otu_mat <- sample_otu_mat[rowSums(sample_otu_mat) > 0, colSums(sample_otu_mat) > 0]     # Remove empty rows & columns
    cat("Computing temporal data for niche:",niche,"\n")
    betapart_temp <- compute_temp_data(sample_otu_mat, sample_meta)
    
    mono <- monotonic_gam(betapart_temp, nK=5)

    ggplot(betapart_temp, aes(distance, jaccard)) +
        geom_point(alpha=0.2, size=2, colour=colour, show.legend=FALSE) +
        theme_linedraw(base_size=20) +
        geom_line(data=mono, aes(distance, pred_fit), lwd=lsize) +
        scale_y_continuous(limits=c(0.4, 1)) +
        scale_x_continuous(limits=c(0,185)) +
        labs(title=niche,x="Days", y="Dissimilarity (J)")
}

p1 <- temporal_niche_plot(sample_otu_matrix_mg, sample_meta_mg, "Saprophage", tax_niche_mg, mg_col)
p2 <- temporal_niche_plot(sample_otu_matrix_mg, sample_meta_mg, "Phytophage", tax_niche_mg, mg_col)
p3 <- temporal_niche_plot(sample_otu_matrix_mg, sample_meta_mg, "Predator", tax_niche_mg, mg_col)
p4 <- temporal_niche_plot(sample_otu_matrix_mg, sample_meta_mg, "Phytophage-parasitoid", tax_niche_mg, mg_col)

q1 <- temporal_niche_plot(sample_otu_matrix_se, sample_meta_se, "Saprophage", tax_niche_se, se_col)
q2 <- temporal_niche_plot(sample_otu_matrix_se, sample_meta_se, "Phytophage", tax_niche_se, se_col)
q3 <- temporal_niche_plot(sample_otu_matrix_se, sample_meta_se, "Predator", tax_niche_se, se_col)
q4 <- temporal_niche_plot(sample_otu_matrix_se, sample_meta_se, "Saprophage-parasitoid", tax_niche_se, se_col)
q5 <- temporal_niche_plot(sample_otu_matrix_se, sample_meta_se, "Phytophage-parasitoid", tax_niche_se, se_col)
q6 <- temporal_niche_plot(sample_otu_matrix_se, sample_meta_se, "Predator-parasitoid", tax_niche_se, se_col)

# Save plots
ggsave( file = "../figs/Fig_temporal_niche_mg.jpg",
        width = 16.0,
        height = 14.0,
        plot = p1 + p2 + p3 + p4 +
            plot_annotation(tag_levels="A") +
            plot_layout(ncol=2, axis_titles="collect")
      )

ggsave( file = "../figs/Fig_temporal_niche_se.jpg",
        width = 23.0,
        height = 14.0,
        plot = q1 + q2 + q3 + q4 + q5 + q6 +
            plot_annotation(tag_levels="A") +
            plot_layout(ncol=3, axis_titles="collect_x")
      )

