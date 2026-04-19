# Compute permuted site accumulation data for habitats
# in total Swedish data

# Preamble
# --------

# Function for filtering traps based on number of samples
filter_traps <- function(sample_meta, min_samples) {
    X <- unique(sample_meta[,c("trapID","sampleID_NGI")])
    trap_samples <- table(X$trapID)
    names(trap_samples[trap_samples >= min_samples])
}

# Generate reproducible random numbers
set.seed(10)


# Read in prepared data
# ---------------------

otu_site_meta_se <- readRDS("../data/otu_site_meta_se.rds")

# Uncomment the lines below to do more stringent filtering on sampling effort (now set at 10 traps)
# Setting the effort at 20 traps would remove some of the northern sites altogether
# traps_to_keep <- filter_se_traps(20)    # Only keep se traps with 20 or more samples
# otu_site_meta_se <- otu_site_meta_se[otu_site_meta_se$trapID %in% traps_to_keep,]


# Compute accumulation data for habitats
# --------------------------------------------

generate_habitat_accumulation_data <- function(D, habitat) {
    res <- data.frame()
    sites <- unique(D$trapID)
    num_sites <- length(sites)
    for (n in 1:num_sites) {
        for (k in 1:100) {
            idx <- sample(1:num_sites,size=n)
            traps <- sites[idx]
            X <- D[D$trapID %in% traps,]
            habitat_otu_prop <- length(unique(X$cluster[X$Habitat==habitat])) / length(unique(X$cluster))
            res <- rbind(res,list(samples=n, habitat_otu_prop=habitat_otu_prop))
        }
    }
    return (res)
}

cat("Computing data for Plants\n")
plants <- generate_habitat_accumulation_data(otu_site_meta_se, "Plants")

cat("Computing data for Soil\n")
soil <- generate_habitat_accumulation_data(otu_site_meta_se, "Soil")

cat("Computing data for Water\n")
water <- generate_habitat_accumulation_data(otu_site_meta_se, "Water")

cat("Computing data for Wood\n")
wood <- generate_habitat_accumulation_data(otu_site_meta_se, "Wood")

cat("Computing data for Temporary habitats\n")
temporary <- generate_habitat_accumulation_data(otu_site_meta_se, "Temporary habitats")

cat("Computing data for Fungi\n")
fungi <- generate_habitat_accumulation_data(otu_site_meta_se, "Fungi")


# Write resulting data files
# --------------------------

write.tsv <- function(D, file) { write.table(D, file, row.names=FALSE, sep="\t") }
write.tsv(plants,"../data/site_acc_plants_prop.tsv")
write.tsv(soil,"../data/site_acc_soil_prop.tsv")
write.tsv(water,"../data/site_acc_water_prop.tsv")
write.tsv(wood,"../data/site_acc_wood_prop.tsv")
write.tsv(temporary,"../data/site_acc_temporary_prop.tsv")
write.tsv(fungi,"../data/site_acc_fungi_prop.tsv")

