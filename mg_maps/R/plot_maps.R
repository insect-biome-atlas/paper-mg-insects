# Plot mg map with sample sites, diversity etc

library(vegan)
library(patchwork)
library(ggplot2)
library(ggnewscale)
library(tidyverse)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(MetBrewer)
library(terra)
library(tidyterra)
library(geodata)
library(scales)
library(gridExtra)
library(ggspatial)


# Read in data
# ------------

# Cluster taxonomy (only Hexapoda, and with quality filtering)
clusters <- readRDS("../../iba_data/cluster_taxonomy_mg.rds")

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

# Read in site metadata
sites_meta <- read.delim("../../iba_data/sites_metadata_MG.tsv")
sites_meta$trap_habitat[sites_meta$trap_habitat=="Dry_Forest"] <- "Dry forest"
sites_meta$trap_habitat[sites_meta$trap_habitat=="Montane_Rainforest"] <- "Montane forest"
sites_meta$trap_habitat[sites_meta$trap_habitat=="Rainforest"] <- "Wet forest"

# Combine metadata and cluster data into otu table at sample level
otu_sample_meta <- merge(taxa_counts, sample_meta, by="sampleID_NGI")

# Aggregate counts at trap level (we will later compute medians at site/location level)
otu_site_meta <- aggregate(read_count~trapID+cluster+Order+Niche+Habitat+trap_habitat+latitude+longitude, data=otu_sample_meta, FUN=sum)


# Compute plot data
# -----------------

# Function to compute site-level data from trap level data
compute_site_data <- function(D, value_col_name) {
    D$siteID <- sites_meta$siteID[match(D$trapID,sites_meta$trapID)]
    D <- aggregate(D[,value_col_name], by=list(siteID=D$siteID), FUN=median)
    colnames(D)[2] <- value_col_name
    idx <- match(D$siteID,sites_meta$siteID)
    D$latitude <- sites_meta$latitude[idx]
    D$longitude <- sites_meta$longitude[idx]
    D$trap_habitat <- sites_meta$trap_habitat[idx]    
    D
}

# Compute site species richness
site_catch <- data.frame(table(otu_site_meta$trapID))
colnames(site_catch) <- c("trapID","OTUs")
site_catch <- compute_site_data(site_catch,"OTUs")

# Compute site diversity (Shannon)
site_otu_abundance <- readRDS("../../iba_data/site_otu_abundance_malaise_mg.rds")
shannon_diversity <- diversity(site_otu_abundance, index="shannon")
site_diversity <- data.frame(shannon_diversity)
site_diversity$trapID <- row.names(site_diversity)
site_diversity <- compute_site_data(site_diversity,"shannon_diversity")

# Compute site uniqueness (based on Bray-Curtis dissimilarities)
dist_matrix <- vegdist(site_otu_abundance, method="bray")
dist_mat <- as.matrix(dist_matrix)
x <- colSums(dist_mat) / (nrow(dist_mat) - 1) # Compute mean dissimilarity, exclude self-distance (=0)
site_uniqueness <- data.frame(list(trapID=names(x),dissimilarity=x))
site_uniqueness <- compute_site_data(site_uniqueness,"dissimilarity")

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
site_colonizations <- compute_site_data(site_colonizations,"num_clades")

# Compute site 'colonization ages' (age of endemic clade)
X <- aggregate(P$placement_age,by=list(edge_num=P$edge_num),FUN=mean)
colnames(X) <- c("edge_num","placement_age")
F <- unique(D[,c("placement","trapID")])
F$placement_age <- X$placement_age[match(F$placement,X$edge_num)]
site_radiation_age <- aggregate(placement_age~trapID, data=F, FUN=median)
colnames(site_radiation_age) <- c("trapID","median_age")
site_radiation_age <- compute_site_data(site_radiation_age,"median_age")

# Compute site 'endemicity' (mean endemic clade size)
site_endemicity <- aggregate(num_placements~trapID, data=E, FUN=mean)
colnames(site_endemicity) <- c("trapID","mean_radiation")
site_endemicity <- compute_site_data(site_endemicity,"mean_radiation")


# Re-generate hillshade map (saving sf objects introduces some annoying errors)
# -----------------------------------------------------------------------------

# Get elevation, slope, and aspect data
mad_elev        <- elevation_30s("MDG" , path = ".") 
mad_slope       <- terrain(mad_elev, "slope", unit = "radians")
mad_aspect      <- terrain(mad_elev, "aspect", unit = "radians")
pal_greys <- hcl.colors(1000, "Grays")
  
# Decide on which slopes to shade
mad_hshade        <- shade(mad_slope, mad_aspect, 30, 270)
names(mad_hshade) <- "shades" # add names to index variable

# Get an index of which values to shade and to what degree
index <- mad_hshade %>%
    mutate(index_col = scales::rescale(shades, to = c(1, length(pal_greys)))) %>% 
    mutate(index_col = round(index_col)) %>%
    pull(index_col)
  
# Get a palette of grey colours to act as the "shade"
shade_cols <- pal_greys[index]

# Make the plot
mad_hshade_plot <- ggplot() +
    geom_spatraster(data = mad_hshade, fill = shade_cols, maxcell = Inf,alpha = .8)   # Avoid resampling with maxcell
  
# Scale min max values of elevation raster
elev_limits <- minmax(mad_elev) %>% as.vector() # Get min max
elev_limits <- c(floor(elev_limits[1] / 500), ceiling(elev_limits[2] / 500)) * 500 # Rounded to lower and upper 500
elev_limits <- pmax(elev_limits, 0) # Set min to 0


# Generate main figure
# --------------------

# Define common figure settings
symbol_size <- 5
symbol_set <- c(21,24,22)
pal_option <- "mako"
legend_pos <- c(0.25,0.93)
dissim_legend_pos <- c(0.30,0.93)   # Axis text of Bray-Curtis bar necessitates right shift
compound_legend_pos <- c(0.25,0.85)
title_size <- 18
legend_text_size <- 11
legend_title_size <- 14
legend_key_width <- unit(0.8,"cm")

p1 <- mad_hshade_plot + # PLot hillshaded map
    geom_spatraster(data = mad_elev, maxcell = Inf , show.legend = FALSE) +
    ggspatial::annotation_scale(location = 'br', width_hint = .4, text_cex = 1.2) +
    scale_fill_hypso_tint_c(limits = elev_limits , palette = "dem_poster",alpha = 0.1,direction = 1) +
    new_scale_fill() +
    geom_point(data=site_catch, 
               aes(x=longitude, y=latitude, fill=OTUs , shape = trap_habitat),
               size = symbol_size) +
    labs(title = "(A) Richness",
         x = NULL,
         y = NULL,
         fill = "# OTUs",
         shape = "Habitat") +
    theme_linedraw() +
    theme(legend.background = element_rect(fill = "white"),
          legend.position = compound_legend_pos,
          plot.title = element_text(size=title_size),
          legend.text = element_text(size=legend_text_size),
          legend.title = element_text(size=legend_title_size),
          legend.key.width = legend_key_width,
          legend.spacing.y = unit(0.0,"cm"),
          axis.text = element_text(size = 12)) +
    guides(fill = guide_colorbar(direction = "horizontal",
                                 title.position = "top", order=1),
           shape = guide_legend(title.position = "top", order=2, override.aes = list(size=3))) +
    scale_shape_manual(values = symbol_set) +
    scale_fill_viridis_c(option = pal_option)

p2 <- mad_hshade_plot + # Plot hillshaded map
    geom_spatraster(data = mad_elev, maxcell = Inf , show.legend = FALSE) +
    scale_fill_hypso_tint_c(limits = elev_limits , palette = "dem_poster",alpha = 0.1,direction = 1) +
    new_scale_fill() +
    geom_point(data=site_diversity, 
               aes(x=longitude, y=latitude, shape=trap_habitat , fill = shannon_diversity),
               size = symbol_size)  +
    labs(title = "(B) Diversity",
         x = NULL,
         y = NULL,
         fill = "Shannon diversity",
         shape = "Habitat") +
    theme_linedraw() +
    theme(legend.background = element_rect(fill = "white"),
          legend.position = legend_pos,
          plot.title = element_text(size=title_size),
          legend.text = element_text(size=legend_text_size),
          legend.title = element_text(size=legend_title_size),
          legend.key.width = legend_key_width,
          axis.text = element_blank()) +
    guides(fill = guide_colorbar(direction = "horizontal",
                                 title.position = "top"),
           shape = "none") +
    scale_shape_manual(values = symbol_set) +
    scale_fill_viridis_c(option = pal_option)

p3 <- mad_hshade_plot + # PLot hillshaded map
    geom_spatraster(data = mad_elev, maxcell = Inf , show.legend = FALSE) +
    scale_fill_hypso_tint_c(limits = elev_limits , palette = "dem_poster",alpha = 0.1,direction = 1) +
    new_scale_fill() +
    geom_point(data=site_uniqueness, 
               aes(x=longitude, y=latitude, shape=trap_habitat , fill = dissimilarity),
               size = symbol_size)  +
    geom_point(size=3, position="jitter") +
    labs(title = "(C) Uniqueness",
         x = NULL,
         y = NULL,
         fill = "Bray-Curtis dissim.",
         shape = "Habitat") +
    theme_linedraw() +
    theme(legend.background = element_rect(fill = "white"),
          legend.position = dissim_legend_pos,
          plot.title = element_text(size=title_size),
          legend.text = element_text(size=legend_text_size),
          legend.title = element_text(size=legend_title_size),
          legend.key.width = legend_key_width,
          axis.text = element_blank()) +
    guides(fill = guide_colorbar(direction = "horizontal",
                                 title.position = "top"),
           shape = "none") +
    scale_shape_manual(values = symbol_set) +
    scale_fill_viridis_c(option = pal_option)
  
p4 <- mad_hshade_plot + # PLot hillshaded map
    geom_spatraster(data = mad_elev, maxcell = Inf , show.legend = FALSE) +
    scale_fill_hypso_tint_c(limits = elev_limits , palette = "dem_poster",alpha = 0.1,direction = 1) +
    new_scale_fill() +
    geom_point(data=site_colonizations, 
               aes(x=longitude, y=latitude, shape=trap_habitat , fill = num_clades),
               size = symbol_size)  +
    geom_point(size=3, position="jitter") +
    labs(title = "(D) Number of radiations",
         x = NULL,
         y = NULL,
         shape = "Habitat",
         fill = "# clades") +
    theme_linedraw() +
    theme(legend.background = element_rect(fill = "white"),
          legend.position = legend_pos, 
          plot.title = element_text(size=title_size),
          legend.text = element_text(size=legend_text_size),
          legend.title = element_text(size=legend_title_size),
          legend.key.width = legend_key_width,
          axis.text = element_blank()) +
    guides(fill = guide_colorbar(direction = "horizontal",
                                 title.position = "top"),
           shape = "none") +
    scale_shape_manual(values = symbol_set) +
    scale_fill_viridis_c(option = pal_option)

p5 <- mad_hshade_plot + # PLot hillshaded map
    geom_spatraster(data = mad_elev, maxcell = Inf , show.legend = FALSE) +
    scale_fill_hypso_tint_c(limits = elev_limits , palette = "dem_poster",alpha = 0.1,direction = 1) +
    new_scale_fill() +
    geom_point(data=site_endemicity, 
               aes(x=longitude, y=latitude, shape=trap_habitat , fill = mean_radiation),
               size = symbol_size)  +
    labs(title = "(E) Size of radiations",
         x = NULL,
         y = NULL,
         shape = "Habitat",
         fill = "Mean # OTUs") +
    theme_linedraw() +
    theme(legend.background = element_rect(fill = "white"),
          legend.position = legend_pos,
          plot.title = element_text(size=title_size),
          legend.text = element_text(size=legend_text_size),
          legend.title = element_text(size=legend_title_size),
          legend.key.width = legend_key_width,
          axis.text = element_blank()) +
    guides(fill = guide_colorbar(direction = "horizontal",
                                 title.position = "top"),
           shape = "none") +
    scale_shape_manual(values = symbol_set) +
    scale_fill_viridis_c(option = pal_option)

p6 <- mad_hshade_plot + # PLot hillshaded map
    geom_spatraster(data = mad_elev, maxcell = Inf , show.legend = FALSE) +
    scale_fill_hypso_tint_c(limits = elev_limits , palette = "dem_poster",alpha = 0.1,direction = 1) +
    new_scale_fill() +
    geom_point(data=site_radiation_age, 
               aes(x=longitude, y=latitude, shape=trap_habitat , fill = median_age),
               size = symbol_size)   +
    labs(title = "(F) Age of radiations",
         x = NULL,
         y = NULL,
         shape = "Habitat",
         fill = "Median age (Ma)") +
    theme_linedraw() +
    theme(legend.background = element_rect(fill = "white"),
          legend.position = legend_pos,
          plot.title = element_text(size=title_size),
          legend.text = element_text(size=legend_text_size),
          legend.title = element_text(size=legend_title_size),
          legend.key.width = legend_key_width,
          axis.text = element_blank(),
          legend.spacing.y = unit(0, "cm")) +
    guides(fill = guide_colorbar(direction = "horizontal",
                                 title.position = "top"),
           shape = "none") +
    scale_shape_manual(values = symbol_set) +
    scale_fill_viridis_c(option = pal_option)

# Plot using patchwork
ggsave("../figs/Fig_mg_maps.jpg",
       width = 21.0,
       height = 15.0,
       plot = ((p1 + p2 + p3) / (p4 + p5 + p6)))


# Generate supplementary figure
# -----------------------------
  
# Add placement age to base data frame
D$placement_age <- F$placement_age[match(D$placement,F$placement)]
  
# Generate plot dataset from age subset input
plot_data <- function(D) {
    E <- data.frame(table(D$placement,by=D$trapID))
    colnames(E) <- c("Placement","trapID","num_placements")
    E <- E[E$num_placements!=0,]
    res <- data.frame(table(E$trapID))
    colnames(res) <- c("trapID","num_clades")
    compute_site_data(res,"num_clades")
}
  
# Plot function
plot_numclades <- function(D, plot_title, lgd_pos=legend_pos, axis_txt=element_blank()) {
    mad_hshade_plot + # PLot hillshaded map
    geom_spatraster(data = mad_elev, maxcell = Inf , show.legend = FALSE) +
#    ggspatial::annotation_scale(location = 'tl',width_hint = .4,text_cex = 1) +
    scale_fill_hypso_tint_c(limits = elev_limits , palette = "dem_poster",alpha = 0.1,direction = 1) +
    new_scale_fill() +
    geom_point(data=D, aes(x=longitude, y=latitude, fill=num_clades, shape=trap_habitat), size = symbol_size) +
    labs(title = plot_title,
         x = NULL,
         y = NULL,
         shape = "Habitat",
         fill = "# clades") +
    theme_linedraw() +
    theme(legend.background = element_rect(fill = "white"),
          legend.position = lgd_pos,
          plot.title = element_text(size=title_size),
          legend.text = element_text(size=legend_text_size),
          legend.title = element_text(size=legend_title_size),
          legend.key.width = legend_key_width,
          axis.text = axis_txt,
          legend.spacing.y = unit(0, "cm")) +
    guides(fill = guide_colorbar(direction = "horizontal",
                                 title.position = "top", order =1),
           shape = 'none') +
    scale_shape_manual(values = symbol_set) +
    scale_fill_viridis_c(option = pal_option)
}

E1 <- plot_data(D[D$placement_age < 23.0,])
E2 <- plot_data(D[D$placement_age > 23.0,])
E3 <- plot_data(D[D$placement_age < 34.0,])
E4 <- plot_data(D[D$placement_age > 34.0,])

supp1 <- plot_numclades(E1, "(A) Young clades (< 23.0 Ma)", compound_legend_pos, element_text(size=12)) +
            ggspatial::annotation_scale(location = 'br',width_hint = .4,text_cex = 1) +
            guides(shape = guide_legend(title.position = "top", order=2, override.aes = list(size=3)))
supp2 <- plot_numclades(E2, "(B) Old clades (> 23.0 Ma)")
supp3 <- plot_numclades(E3, "(C) Young clades (< 34.0 Ma)")
supp4 <- plot_numclades(E4, "(D) Old clades (> 34.0 Ma)")

# Put plots together and save
ggsave("../figs/Fig_mg_maps_old_vs_young.jpg",
       width=21.0,
       height=15.0,
       plot = (supp1 + supp2) / (supp3 + supp4))

