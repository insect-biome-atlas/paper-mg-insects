# Compute permuted site accumulation data for habitats

source("filter_se_traps.R")

# Generate reproducible random numbers
set.seed(10)


# Read in prepared data
# ---------------------

otu_site_meta_mg <- read.delim("../data/otu_site_meta_mg.tsv")
otu_site_meta_se <- read.delim("../data/otu_site_meta_se.tsv")
traps_to_keep <- filter_se_traps(20)    # Only keep se traps with 20 or more samples
otu_site_meta_se <- otu_site_meta_se[otu_site_meta_se$trapID %in% traps_to_keep,]


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
plants_mg <- generate_habitat_accumulation_data(otu_site_meta_mg, "Plants")
plants_mg$country <- "Madagascar"
plants_se <- generate_habitat_accumulation_data(otu_site_meta_se, "Plants")
plants_se$country <- "Sweden"
plants <- rbind(plants_mg, plants_se)

cat("Computing data for Soil\n")
soil_mg <- generate_habitat_accumulation_data(otu_site_meta_mg, "Soil")
soil_mg$country <- "Madagascar"
soil_se <- generate_habitat_accumulation_data(otu_site_meta_se, "Soil")
soil_se$country <- "Sweden"
soil <- rbind(soil_mg, soil_se)

cat("Computing data for Water\n")
water_mg <- generate_habitat_accumulation_data(otu_site_meta_mg, "Water")
water_mg$country <- "Madagascar"
water_se <- generate_habitat_accumulation_data(otu_site_meta_se, "Water")
water_se$country <- "Sweden"
water <- rbind(water_mg, water_se)

cat("Computing data for Wood\n")
wood_mg <- generate_habitat_accumulation_data(otu_site_meta_mg, "Wood")
wood_mg$country <- "Madagascar"
wood_se <- generate_habitat_accumulation_data(otu_site_meta_se, "Wood")
wood_se$country <- "Sweden"
wood <- rbind(wood_mg, wood_se)

cat("Computing data for Temporary habitats\n")
temporary_mg <- generate_habitat_accumulation_data(otu_site_meta_mg, "Temporary habitats")
temporary_mg$country <- "Madagascar"
temporary_se <- generate_habitat_accumulation_data(otu_site_meta_se, "Temporary habitats")
temporary_se$country <- "Sweden"
temporary <- rbind(temporary_mg, temporary_se)

cat("Computing data for Fungi\n")
fungi_mg <- generate_habitat_accumulation_data(otu_site_meta_mg, "Fungi")
fungi_mg$country <- "Madagascar"
fungi_se <- generate_habitat_accumulation_data(otu_site_meta_se, "Fungi")
fungi_se$country <- "Sweden"
fungi <- rbind(fungi_mg, fungi_se)

write.tsv <- function(D, file) { write.table(D, file, row.names=FALSE, sep="\t") }
write.tsv(plants,"../data/site_acc_plants_prop.tsv")
write.tsv(soil,"../data/site_acc_soil_prop.tsv")
write.tsv(water,"../data/site_acc_water_prop.tsv")
write.tsv(wood,"../data/site_acc_wood_prop.tsv")
write.tsv(temporary,"../data/site_acc_temporary_prop.tsv")
write.tsv(fungi,"../data/site_acc_fungi_prop.tsv")

