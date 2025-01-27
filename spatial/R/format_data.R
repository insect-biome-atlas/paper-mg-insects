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

# Following data is required from the figshare

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
site_meta_se <- fread("data/sites_metadata_SE_2019.tsv") 
site_meta_mg <- fread("data/sites_metadata_MG_2019.tsv") 

# Assemble the metadata
full_meta_DT_se <-  merge(seq_meta_se , trap_meta_se , all=TRUE) 
full_meta_DT_se <-  merge(full_meta_DT_se , site_meta_se,by="trapID") |> _[lab_sample_type %in% "sample",]
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


# save ----------------------------------------------------------------------------------------

saveRDS(OTU_DT_se , "data/species_table_se.rds")
saveRDS(OTU_DT_mg , "data/species_table_mg.rds")
