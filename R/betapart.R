# Load required libraries
library(tidyverse)
library(data.table)
library(sf)
library(units)
library(betapart)
library(mgcv)

# Set random seed
set.seed(10)

# ---------------------------------------------------------------------------------------------
# data ----------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------

# Spike in IDs 
spike_ins_mg  <- fread("data/Biological_spikes_MG_taxonomy.tsv")
spike_ins_se  <- fread("data/Biological_spikes_SE_taxonomy.tsv")

# ----------------- Get and clean cluster data
swe_clusters <- fread("data/non_numts_cleaned_SE.tsv") |>
  _[representative == 1 ,] |> 
  _[Phylum == "Arthropoda" ,] |> 
  _[!str_detect(Family , "_X*"),] |> 
  _[!str_detect(Species , paste(spike_ins_se$Species,collapse="|")),]

mg_clusters  <- fread("data/non_numts_cleaned_MG.tsv") |>
  _[representative == 1 ,] |> 
  _[Phylum == "Arthropoda" ,] |> 
  _[!str_detect(Family , "_X*"),] |> 
  _[!str_detect(Species , paste(spike_ins_mg$Species,collapse="|")),]

# --------------- Get sample IDs 
seq_meta_se <- fread("data/co1_sequencing_metadata_se.tsv")
seq_meta_mg <- fread("data/co1_sequencing_metadata_mg.tsv")
trap_meta_se <- fread("data/samples_metadata_malaise_SE_2019.tsv")
trap_meta_mg <- fread("data/samples_metadata_malaise_MG_2019.tsv")
site_meta_se <- fread("data/sites_metadata_SE_2019.tsv") |> st_as_sf(coords = c("longitude_WGS84" , "latitude_WGS84") , crs = 4326) 
site_meta_mg <- fread("data/sites_metadata_MG_2019.tsv") |> st_as_sf(coords = c("longitude_WGS84" , "latitude_WGS84") , crs = 4326) 


# Assemble the metadata
full_meta_DT_se <-  merge(seq_meta_se , trap_meta_se , all=TRUE) 
full_meta_DT_se <- merge(full_meta_DT_se , site_meta_se,by="trapID") |> _[lab_sample_type %in% "sample",]
full_meta_DT_mg <-  merge(seq_meta_mg , trap_meta_mg , all=TRUE) 
full_meta_DT_mg <- merge(full_meta_DT_mg , site_meta_mg,by="trapID") |> _[lab_sample_type %in% "sample",]

# ------------ Get read numbers 
swe_counts   <- fread("data/non_numts_cluster_counts_cleaned_SE.tsv")
mg_counts    <- fread("data/non_numts_cluster_counts_cleaned_MG.tsv")

# ---------------------------------------------------------------------------------------------
# tidy data -----------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------

# First get counts of each species in each sample
swe_counts_long <- swe_counts |> 
                   melt(id.vars = 1 , variable.name = "sampleID_NGI" , value.name = "read_count") |> 
                  _[read_count > 0,]

mg_counts_long <- mg_counts |> 
                  melt(id.vars = 1 , variable.name = "sampleID_NGI" , value.name = "read_count") |> 
                  _[read_count > 0,]

# Merge with cluster IDS ----------------------- 
se_counts_DT <- swe_counts_long[swe_clusters, on = "cluster"] |> 
                _[!str_detect(Order , "unclassified"),] 

mg_counts_DT <- mg_counts_long[mg_clusters, on = "cluster"] |> 
                _[!str_detect(Order , "unclassified"),] 

# Combine metadata and OTU data ----------------------- 
OTU_DT_se <- merge(se_counts_DT , full_meta_DT_se, by = "sampleID_NGI") |> 
             _[lab_sample_type == "sample", ] |> 
             group_by(trapID , cluster) |> 
             summarise(pres = 1*(sum(read_count)>1))

OTU_DT_mg <- merge(mg_counts_DT , full_meta_DT_mg, by = "sampleID_NGI") |> 
              _[lab_sample_type == "sample", ] |> 
              group_by(trapID , cluster) |> 
              summarise(pres = 1*(sum(read_count)>1))

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

dist_swe <- get_dists(site_meta_se) 
dist_mad <- get_dists(site_meta_mg)

# partition ---------------------------------------------------------------

betapart_swe <- partition_beta_diversity(sp_matrix_total_se, trap_dist = dist_swe)
betapart_mad <- partition_beta_diversity(sp_matrix_total_mg , trap_dist = dist_mad)

# monotonic gam fits --------------------------------------------------------------------------


# Fits a single cubic regression spline to the distance / turnover data
mono_se <- monotonic_gam(betapart_swe , nK = 5)
mono_mg <- monotonic_gam(betapart_mad , nK = 5)

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

p1 + p2 

