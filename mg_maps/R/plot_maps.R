# Plot mg map with sample sites, diversity etc

library(vegan)
library(patchwork)
library(ggplot2)
library(patchwork)
library(ggnewscale)


# Read in data
# ------------
# mdg hillshade plot 
mad_hshade_plot <- readRDS("mg_maps/R/mdg_hs_map") %>% 
  
  # Cluster taxonomy (only Hexapoda, and with quality filtering)
  clusters <- readRDS("iba_data/cluster_taxonomy_mg.rds")
  
  # Add info on life-history traits
  lht <- read.delim("traits/clade_trait_data_mg.tsv")
  
  # Merge cluster taxonomy and life-history traits
  clusters <- merge(clusters, lht)
  
  # Read in cluster read numbers for malaise traps and litter samples
  malaise_counts_long <- readRDS("iba_data/cluster_counts_malaise_long_mg.rds")
  litter_counts_long <- readRDS("iba_data/cluster_counts_litter_long_mg.rds")
  
  # Merge cluster read numbers
  counts_long <- rbind(malaise_counts_long, litter_counts_long)
  
  
  
  # Tidy data
  # ---------
  
  # Merge counts with cluster taxonomy
  taxa_counts <- merge(counts_long, clusters, by="cluster")
  
  # Read in collected sample metadata
  sample_meta <- read.delim("iba_data/malaise_litter_sample_meta_mg.tsv")
  
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
  site_otu_abundance <- readRDS("iba_data/site_otu_abundance_malaise_mg.rds")
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
  P <- readRDS("placements/data/placement_stats.rds")
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
  
  
  
  # Re-generate maps (saving sf objects introduces some annoying errors) --------
  
  # same but for madagascar -------------------------------------------------
  
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
  
  
  
  # -------------------------------------------------------------------------
  
  
  
  # Generate main figure
  # --------------------
  
  p1 <- mad_hshade_plot + # PLot hillshaded map
    geom_spatraster(data = mad_elev, maxcell = Inf , show.legend = FALSE) +
    ggspatial::annotation_scale(location = 'tl',width_hint = .4,text_cex = 1)+
    scale_fill_hypso_tint_c(limits = elev_limits , palette = "dem_poster",alpha = 0.1,direction = 1)+
    new_scale_fill() +
    geom_point(data=site_catch, 
               aes(x=longitude, y=latitude, fill=OTUs , shape = trap_habitat),
               size = 3) +
    labs(title = "Richness",
         x = NULL,
         y = NULL,
         fill = "# OTUs",
         shape = "Habitat")+
    theme_linedraw()+
    theme(legend.background = element_rect(fill = "transparent"),
          legend.position = c(0.25, 0.90), 
          axis.text = element_text(size = 12)) +
    guides(fill = guide_colorbar(direction = "horizontal",
                                 title.position = "top"),
           shape = "none")+
    scale_shape_manual(values = c(24, 21, 25))+
    scale_fill_viridis_c(option = "mako")
  
  p2 <- mad_hshade_plot + # PLot hillshaded map
    geom_spatraster(data = mad_elev, maxcell = Inf , show.legend = FALSE) +
    ggspatial::annotation_scale(location = 'tl',width_hint = .4,text_cex = 1)+
    scale_fill_hypso_tint_c(limits = elev_limits , palette = "dem_poster",alpha = 0.1,direction = 1)+
    new_scale_fill()+
    geom_point(data=site_diversity, 
               aes(x=longitude, y=latitude, shape=trap_habitat , fill = shannon_diversity),
               size = 3)  +
    labs(title = "Diversity",
         x = NULL,
         y = NULL,
         fill = "Shannon",
         shape = "Habitat")+
    theme_linedraw()+
    theme(legend.background = element_rect(fill = "transparent"),
          legend.position = c(0.25, 0.90),
          axis.text = element_text(size = 12)) +
    guides(fill = guide_colorbar(direction = "horizontal",
                                 title.position = "top"),
           shape = "none")+
    scale_shape_manual(values = c(24, 21, 25))+
    scale_fill_viridis_c(option = "mako")
  
  p3 <- mad_hshade_plot + # PLot hillshaded map
    geom_spatraster(data = mad_elev, maxcell = Inf , show.legend = FALSE) +
    ggspatial::annotation_scale(location = 'tl',width_hint = .4,text_cex = 1)+
    scale_fill_hypso_tint_c(limits = elev_limits , palette = "dem_poster",alpha = 0.1,direction = 1)+
    new_scale_fill()+
    geom_point(data=site_uniqueness, 
               aes(x=longitude, y=latitude, shape=trap_habitat , fill = dissimilarity),
               size = 3)  +
    geom_point(size=3, position="jitter") +
    labs(title = "Diversity",
         x = NULL,
         y = NULL,
         fill = "Uniqueness",
         shape = "Habitat")+
    theme_linedraw()+
    theme(legend.background = element_rect(fill = "transparent"),
          legend.position = c(0.25, 0.90), 
          axis.text = element_text(size = 12)) +
    guides(fill = guide_colorbar(direction = "horizontal",
                                 title.position = "top"),
           shape = "none")+
    scale_shape_manual(values = c(24, 21, 25))+
    scale_fill_viridis_c(option = "mako")
  
  p4 <- mad_hshade_plot + # PLot hillshaded map
    geom_spatraster(data = mad_elev, maxcell = Inf , show.legend = FALSE) +
    ggspatial::annotation_scale(location = 'tl',width_hint = .4,text_cex = 1)+
    scale_fill_hypso_tint_c(limits = elev_limits , palette = "dem_poster",alpha = 0.1,direction = 1)+
    new_scale_fill()+
    geom_point(data=site_colonizations, 
               aes(x=longitude, y=latitude, shape=trap_habitat , fill = num_clades),
               size = 3)  +
    geom_point(size=3, position="jitter") +
    labs(title = "Number of Radiations",
         x = NULL,
         y = NULL,
         shape = "Habitat",
         fill = "# Clades")+
    theme_linedraw()+
    theme(legend.background = element_rect(fill = "transparent"),
          legend.position = c(0.25, 0.90), 
          axis.text = element_text(size = 12)) +
    guides(fill = guide_colorbar(direction = "horizontal",
                                 title.position = "top"),
           shape = "none")+
    scale_shape_manual(values = c(24, 21, 25))+
    scale_fill_viridis_c(option = "mako")
  
  
  
  
  p5 <- mad_hshade_plot + # PLot hillshaded map
    geom_spatraster(data = mad_elev, maxcell = Inf , show.legend = FALSE) +
    ggspatial::annotation_scale(location = 'tl',width_hint = .4,text_cex = 1)+
    scale_fill_hypso_tint_c(limits = elev_limits , palette = "dem_poster",alpha = 0.1,direction = 1)+
    new_scale_fill()+
    geom_point(data=site_endemicity, 
               aes(x=longitude, y=latitude, shape=trap_habitat , fill = mean_radiation),
               size = 3)  +
    labs(title = "Size of Radiations",
         x = NULL,
         y = NULL,
         shape = "Habitat",
         fill = "Mean # OTUs")+
    theme_linedraw()+
    theme(legend.background = element_rect(fill = "transparent"),
          legend.position = c(0.25, 0.90),
          axis.text = element_text(size = 12)) +
    guides(fill = guide_colorbar(direction = "horizontal",
                                 title.position = "top"),
           shape = "none")+
    scale_shape_manual(values = c(24, 21, 25))+
    scale_fill_viridis_c(option = "mako")
  
  
  
  p6 <-  mad_hshade_plot + # PLot hillshaded map
    geom_spatraster(data = mad_elev, maxcell = Inf , show.legend = FALSE) +
    ggspatial::annotation_scale(location = 'tl',width_hint = .4,text_cex = 1)+
    scale_fill_hypso_tint_c(limits = elev_limits , palette = "dem_poster",alpha = 0.1,direction = 1)+
    new_scale_fill()+
    geom_point(data=site_radiation_age, 
               aes(x=longitude, y=latitude, shape=trap_habitat , fill = median_age),
               size = 3)   +
    labs(title = "Median age of Radiations",
         x = NULL,
         y = NULL,
         shape = "Habitat",
         fill = "Median age (Ma)")+
    theme_linedraw()+
    theme(legend.background = element_rect(fill = "transparent"),
          legend.position = c(0.25, 0.82),
          axis.text = element_text(size = 12),
          legend.spacing.y = unit(0, "cm")) +
    guides(fill = guide_colorbar(direction = "horizontal",
                                 title.position = "top"),
           shape = guide_legend(title.position = "top"))+
    scale_shape_manual(values = c(24, 21, 25))+
    scale_fill_viridis_c(option = "mako")
  
  # Put plots together and save
  # ggsave("../figs/Fig_mg_maps.jpg",
  #        width=18,
  #        height=14,
  #        plot = p1 + p2 + p3 + p4 + p5 + p6 +
  #             plot_layout(axis_titles="collect", ncol=3) +
  #             plot_annotation(tag_levels="A") #  & theme(legend.position="bottom")
  #        )
  
  # Plot using patchwork
  
  tiff("mg_maps/figs/Fig_mg_maps.tiff", width = 1400, height = 1000, compression = "lzw")
  ((p1 + p2 + p3) / (p4 + p5 + p6)) 
  
  dev.off()
  
  browseURL("mg_maps/figs/Fig_mg_maps.tiff")
  
  # Generate supplementary figure
  # -----------------------------
  
  # Add placement age to base data frame
  D$placement_age <- F$placement_age[match(D$placement,F$placement)]
  
  # Generate plot dataset from age subset input
  plot_data <- function(D, trap_meta) {
    E <- data.frame(table(D$placement,by=D$trapID))
    colnames(E) <- c("Placement","trapID","num_placements")
    E <- E[E$num_placements!=0,]
    res <- data.frame(table(E$trapID))
    colnames(res) <- c("trapID","num_clades")
    idx <- match(res$trapID,trap_meta$trapID)
    res$latitude <- trap_meta$latitude[idx]
    res$longitude <- trap_meta$longitude[idx]
    res$trap_habitat <- trap_meta$trap_habitat[idx]
    return (res)
  }
  
  # Plot function
  plot_numclades <- function(D, plot_title) {
    mad_hshade_plot + # PLot hillshaded map
      geom_spatraster(data = mad_elev, maxcell = Inf , show.legend = FALSE) +
      ggspatial::annotation_scale(location = 'tl',width_hint = .4,text_cex = 1)+
      scale_fill_hypso_tint_c(limits = elev_limits , palette = "dem_poster",alpha = 0.1,direction = 1)+
      new_scale_fill()+
      geom_point(data=D, aes(x=longitude, y=latitude, fill=num_clades, shape=trap_habitat), size = 3) +
      labs(title = plot_title,
           x = NULL,
           y = NULL,
           shape = "Habitat",
           fill = "# Clades")+
      theme_linedraw()+
      theme(legend.background = element_rect(fill = "transparent"),
            legend.position = c(0.25, 0.85),
            axis.text = element_text(size = 12),
            legend.spacing.y = unit(0, "cm"),
            legend.key.height = unit(0.4, "cm")) +
      guides(fill = guide_colorbar(direction = "horizontal",
                                   title.position = "top", order =1),
             shape = 'none')+
      scale_shape_manual(values = c(24, 21, 25))+
      scale_fill_viridis_c(option = "mako")
  }
  
  E1 <- plot_data(D[D$placement_age < 23.0,], trap_meta)
  E2 <- plot_data(D[D$placement_age > 23.0,], trap_meta)
  E3 <- plot_data(D[D$placement_age < 34.0,], trap_meta)
  E4 <- plot_data(D[D$placement_age > 34.0,], trap_meta)
  
  supp1 <- plot_numclades(E1, "Number of young clades (< 23.0 Ma)")
  supp2 <- plot_numclades(E2, "Number of old clades (> 23.0 Ma)")
  supp3 <- plot_numclades(E3, "Number of young clades (< 34.0 Ma)")
  supp4 <- plot_numclades(E4, "Number of old clades (> 34.0 Ma)") + guides(shape = guide_legend(title.position = "top"))
  
  # Put plots together and save
  tiff("mg_maps/figs/Fig_mg_maps_old_vs_young.tiff", width = 1400, height = 1000, compression = "lzw")
  (supp1 + supp2) / (supp3 + supp4)
  dev.off()
  
  browseURL("mg_maps/figs/Fig_mg_maps_old_vs_young.tiff")
  
  