# Process counts data
library(data.table)

# Read in raw counts data
malaise_mg <- read.delim("raw_counts_malaise_mg.tsv")
malaise_se <- read.delim("raw_counts_malaise_se.tsv")
litter_mg <- read.delim("raw_counts_litter_mg.tsv")
litter_se <- read.delim("raw_counts_litter_se.tsv")
malaise_spec_mg <- read.delim("raw_counts_malaise_spec_mg.tsv")
malaise_spec_se <- read.delim("raw_counts_malaise_spec_se.tsv")
malaise_litter_mg <- read.delim("raw_counts_uncal_malaise_litter_mg.tsv")
malaise_litter_se <- read.delim("raw_counts_uncal_malaise_litter_se.tsv")

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

# Read in taxonomy data
T1 <- read.delim("cluster_taxonomy_mg.tsv")
T2 <- read.delim("cluster_taxonomy_se.tsv")

# Only keep the hexapod data (the ones in T1,T2)
malaise_mg <- malaise_mg[malaise_mg$cluster %in% T1$cluster,]
malaise_se <- malaise_se[malaise_se$cluster %in% T2$cluster,]
litter_mg <- litter_mg[litter_mg$cluster %in% T1$cluster,]
litter_se <- litter_se[litter_se$cluster %in% T2$cluster,]
malaise_spec_mg <- malaise_spec_mg[malaise_spec_mg$cluster %in% T1$cluster,]
malaise_spec_se <- malaise_spec_se[malaise_spec_se$cluster %in% T2$cluster,]
malaise_litter_mg <- malaise_litter_mg[malaise_litter_mg$cluster %in% T1$cluster,]
malaise_litter_se <- malaise_litter_se[malaise_litter_se$cluster %in% T2$cluster,]

# Convert counts to occurrence
malaise_litter_mg[,-1] <- (malaise_litter_mg[,-1] > 0)
malaise_litter_se[,-1] <- (malaise_litter_se[,-1] > 0)

# Write fat tables
write.tsv <- function(D,file) { write.table(D,file,row.names=FALSE,sep="\t") }
write.tsv(malaise_mg,"cluster_counts_malaise_mg.tsv")
write.tsv(malaise_se,"cluster_counts_malaise_se.tsv")
write.tsv(litter_mg,"cluster_counts_litter_mg.tsv")
write.tsv(litter_se,"cluster_counts_litter_se.tsv")
write.tsv(malaise_spec_mg,"cluster_counts_spec_malaise_mg.tsv")
write.tsv(malaise_spec_se,"cluster_counts_spec_malaise_se.tsv")
write.tsv(malaise_litter_mg,"cluster_occurrence_malaise_litter_mg.tsv")
write.tsv(malaise_litter_se,"cluster_occurrence_malaise_litter_se.tsv")

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

# Write long tables
write.tsv(malaise_long_mg,"cluster_counts_malaise_long_mg.tsv")
write.tsv(malaise_long_se,"cluster_counts_malaise_long_se.tsv")
write.tsv(litter_long_mg,"cluster_counts_litter_long_mg.tsv")
write.tsv(litter_long_se,"cluster_counts_litter_long_se.tsv")

