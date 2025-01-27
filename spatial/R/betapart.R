# Load required libraries
library(tidyverse)
library(data.table)
library(sf)
library(units)
library(betapart)
library(mgcv)
library(patchwork)

# Set random seed
set.seed(10)
source("R/functions.R")

# ---------------------------------------------------------------------------------------------
# tidy data -----------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------

# load species tables  ----------------------- 
OTU_DT_se <- readRDS("data/species_table_se.rds")
OTU_DT_mg <- readRDS("data/species_table_mg.rds")

# Meta data
# TODO: Change to local path to figshare repo for metadata
metadata_path <- "~/dev/figshare-repos/iba/raw_data/"
site_meta_mg <- fread(paste0(metadata_path,"sites_metadata_MG.tsv"))
site_meta_se <- fread(paste0(metadata_path,"sites_metadata_SE.tsv")) 

# Get species * site matrices  ----------------------- 
sp_matrix_se <- OTU_DT_se |>  
          pivot_wider(names_from = cluster, values_from = pres , values_fill = 0) |> 
          column_to_rownames("trapID")

sp_matrix_mg <- OTU_DT_mg |> 
                pivot_wider(names_from = cluster, values_from = pres , values_fill = 0) |> 
                column_to_rownames("trapID")


# Filter singletons ----------------------- 
sp_matrix_total_se <- sp_matrix_se[,colSums(sp_matrix_se) > 0]
sp_matrix_total_mg <- sp_matrix_mg[,colSums(sp_matrix_mg) > 0]

# ---------------------------------------------------------------------------------------------
# Turnover analysis ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------

# get distances between traps ---------------------------------------------

sf_meta_se <- site_meta_se |> st_as_sf(coords = c("longitude_WGS84" , "latitude_WGS84") , crs = 4326) 
sf_meta_mg <- site_meta_mg |> st_as_sf(coords = c("longitude_WGS84" , "latitude_WGS84") , crs = 4326) 

dist_swe <- get_dists(sf_meta_se) 
dist_mad <- get_dists(sf_meta_mg)

# partition ---------------------------------------------------------------

betapart_swe <- partition_beta_diversity(sp_matrix_total_se, trap_dist = dist_swe)
betapart_mad <- partition_beta_diversity(sp_matrix_total_mg , trap_dist = dist_mad)

# monotonic gam fits --------------------------------------------------------------------------


# Fits a single cubic regression spline to the distance / turnover data
mono_se <- monotonic_gam(betapart_swe , nK = 5)
mono_mg <- monotonic_gam(betapart_mad , nK = 5)

# Render plots
p1 <- ggplot(betapart_swe , aes(distance , jaccard))+
  geom_point(alpha = .05 , size = 2 , aes(colour = distance) , show.legend = FALSE)+
  theme_linedraw(base_size = 20)+
  geom_line(data=mono_se , aes(distance,pred_fit) , lwd=2)+
  scale_colour_viridis_c(option="rocket")+
  scale_y_continuous(limits = c(0.4 , 1))+
  labs(x = "Distance (km)" , y = "Dissimilarity (J)" , colour = "Distance (km)") 

p2 <- ggplot(betapart_mad , aes(distance , jaccard))+
  geom_point(alpha = .2 , size = 2 , aes(colour = distance) , show.legend = FALSE)+
  theme_linedraw(base_size = 20)+
  geom_line(data=mono_mg , aes(distance,pred_fit) , lwd=2)+
  scale_colour_viridis_c(option="rocket")+
  scale_y_continuous(limits = c(0.4 , 1))+
  labs(x = "Distance (km)" , y = "Dissimilarity (J)" , colour = "Distance (km)") 

# Print Plots
print(p1 + p2) 

