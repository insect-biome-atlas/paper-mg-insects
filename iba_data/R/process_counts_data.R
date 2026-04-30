# Process counts data
library(data.table)

# Read in raw counts data
malaise_mg <- readRDS("../raw_counts_malaise_mg.rds")
malaise_se <- readRDS("../raw_counts_malaise_se.rds")
litter_mg <- readRDS("../raw_counts_litter_mg.rds")
litter_se <- readRDS("../raw_counts_litter_se.rds")
malaise_spec_mg <- readRDS("../raw_counts_malaise_spec_mg.rds")
malaise_spec_se <- readRDS("../raw_counts_malaise_spec_se.rds")
malaise_litter_mg <- readRDS("../raw_counts_uncal_malaise_litter_mg.rds")
malaise_litter_se <- readRDS("../raw_counts_uncal_malaise_litter_se.rds")
homogenate_se <- readRDS("../raw_counts_uncal_homogenate_se.rds")
lysate_homogenate_se <- readRDS("../raw_counts_uncal_lysate_homogenate_se.rds")

# Remove taxonomy annotations. They are the sintax+
# annotations, not the phylogenetic annotations, and
# they are redundant anyway. Only keep "cluster" and
# sample counts.
min_columns <- function(D) {
    idx <- which(grepl("_",colnames(D)))
    idx <- c(which(colnames(D)=="cluster"),idx)
    D <- D[,idx]
    idx <- which(colnames(D)=="BOLD_bin")
    return (D[,-idx])
}
malaise_mg <- min_columns(malaise_mg)
malaise_se <- min_columns(malaise_se)
litter_mg <- min_columns(litter_mg)
litter_se <- min_columns(litter_se)
malaise_spec_mg <- min_columns(malaise_spec_mg)
malaise_spec_se <- min_columns(malaise_spec_se)
malaise_litter_mg <- min_columns(malaise_litter_mg) 
malaise_litter_se <- min_columns(malaise_litter_se)
homogenate_se <- min_columns(homogenate_se)
lysate_homogenate_se <- min_columns(lysate_homogenate_se)

# Read in taxonomy data
T1 <- readRDS("../cluster_taxonomy_mg.rds")
T2 <- readRDS("../cluster_taxonomy_se.rds")

# Only keep the hexapod data (the ones in T1,T2)
malaise_mg <- malaise_mg[malaise_mg$cluster %in% T1$cluster,]
malaise_se <- malaise_se[malaise_se$cluster %in% T2$cluster,]
litter_mg <- litter_mg[litter_mg$cluster %in% T1$cluster,]
litter_se <- litter_se[litter_se$cluster %in% T2$cluster,]
malaise_spec_mg <- malaise_spec_mg[malaise_spec_mg$cluster %in% T1$cluster,]
malaise_spec_se <- malaise_spec_se[malaise_spec_se$cluster %in% T2$cluster,]
malaise_litter_mg <- malaise_litter_mg[malaise_litter_mg$cluster %in% T1$cluster,]
malaise_litter_se <- malaise_litter_se[malaise_litter_se$cluster %in% T2$cluster,]
homogenate_se <- homogenate_se[homogenate_se$cluster %in% T2$cluster,]
lysate_homogenate_se <- lysate_homogenate_se[lysate_homogenate_se$cluster %in% T2$cluster,]

# Convert counts to occurrence
malaise_litter_mg[,-1] <- (malaise_litter_mg[,-1] > 0)
malaise_litter_se[,-1] <- (malaise_litter_se[,-1] > 0)

# Save fat tables
saveRDS(malaise_mg,"../cluster_counts_malaise_mg.rds")
saveRDS(malaise_se,"../cluster_counts_malaise_se.rds")
saveRDS(litter_mg,"../cluster_counts_litter_mg.rds")
saveRDS(litter_se,"../cluster_counts_litter_se.rds")
saveRDS(homogenate_se,"../cluster_counts_homogenate_se.rds")
saveRDS(lysate_homogenate_se,"../cluster_counts_lysate_homogenate_se.rds")
saveRDS(malaise_spec_mg,"../cluster_counts_spec_malaise_mg.rds")
saveRDS(malaise_spec_se,"../cluster_counts_spec_malaise_se.rds")
saveRDS(malaise_litter_mg,"../cluster_occurrence_malaise_litter_mg.rds")
saveRDS(malaise_litter_se,"../cluster_occurrence_malaise_litter_se.rds")

# Make long tables (data.table is handy for this)
malaise_long_mg <- data.table(malaise_mg) |>
                   melt(id.vars = 1 , variable.name = "sampleID_NGI" , value.name = "read_count") |>
                   _[read_count > 0,]
malaise_long_se <- data.table(malaise_se) |>
                   melt(id.vars = 1 , variable.name = "sampleID_NGI" , value.name = "read_count") |>
                   _[read_count > 0,]
litter_long_mg <- data.table(litter_mg) |>
                  melt(id.vars = 1 , variable.name = "sampleID_NGI" , value.name = "read_count") |>
                  _[read_count > 0,]
litter_long_se <- data.table(litter_se) |>
                  melt(id.vars = 1 , variable.name = "sampleID_NGI" , value.name = "read_count") |>
                  _[read_count > 0,]
homogenate_long_se <- data.table(homogenate_se) |>
                      melt(id.vars = 1 , variable.name = "sampleID_NGI" , value.name = "read_count") |>
                      _[read_count > 0,]

# Save long tables
saveRDS(malaise_long_mg,"../cluster_counts_malaise_long_mg.rds")
saveRDS(malaise_long_se,"../cluster_counts_malaise_long_se.rds")
saveRDS(litter_long_mg,"../cluster_counts_litter_long_mg.rds")
saveRDS(litter_long_se,"../cluster_counts_litter_long_se.rds")
saveRDS(homogenate_long_se,"../cluster_counts_homogenate_long_se.rds")

