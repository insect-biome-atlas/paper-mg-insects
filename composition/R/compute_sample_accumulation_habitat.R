# Compute permuted sample accumulation data for habitats

# Generate reproducible random numbers
set.seed(10)

# Read in and prepare data
# ------------------------

# Read in data
otu_sample_meta_mg <- read.delim("../data/otu_sample_meta_mg.tsv")
otu_sample_meta_se <- read.delim("../data/otu_sample_meta_se.tsv")


# Restrict to traps with sufficient number of temporal samples
X1 <- unique(otu_sample_meta_mg[,c("trapID","sampleID_NGI")])
X2 <- unique(otu_sample_meta_se[,c("trapID","sampleID_NGI")])
x1 <- table(X1$trapID)
x2 <- table(X2$trapID)
x1 <- x1[x1>=40]
x2 <- x2[x2>=20]
otu_sample_meta_mg <- otu_sample_meta_mg[otu_sample_meta_mg$trapID %in% names(x1),]
otu_sample_meta_se <- otu_sample_meta_se[otu_sample_meta_se$trapID %in% names(x2),]


# Compute accumulation data for major habitats
# --------------------------------------------

generate_habitat_accumulation_data <- function(D, habitat, num_samples) {
    pb <- txtProgressBar(min=0, max=num_samples, style=3, width=50)
    res <- data.frame()
    sites <- unique(D$trapID)
    for (n in 1:num_samples) {
        setTxtProgressBar(pb, n)
        for (k in 1:100) {
            samples <- character()
            for (site in sites) {
                samples <- c(samples, sample(unique(D$sampleID_NGI[D$trapID==site]),size=n))
            }
            X <- D[D$sampleID_NGI %in% samples,]
            habitat_otu_prop <- length(unique(X$cluster[X$Habitat==habitat])) / length(unique(X$cluster))
            res <- rbind(res,list(samples=n, habitat_otu_prop=habitat_otu_prop))
        }
    }
    cat("\n")
    return (res)
}

cat("Computing data for Plants\n")
habitat <- "Plants"
plant_mg <- generate_habitat_accumulation_data(otu_sample_meta_mg, habitat, 40)
plants_mg$country <- "Madagascar"
plants_se <- generate_habitat_accumulation_data(otu_sample_meta_se, habitat, 20)
plants_se$country <- "Sweden"
plants <- rbind(plants_mg, plants_se)

cat("Computing data for Soil\n")
habitat <- "Soil"
soil_mg <- generate_habitat_accumulation_data(otu_sample_meta_mg, habitat, 40)
soil_mg$country <- "Madagascar"
soil_se <- generate_habitat_accumulation_data(otu_sample_meta_se, habitat, 20)
soil_se$country <- "Sweden"
soil <- rbind(soil_mg, soil_se)

cat("Computing data for Water\n")
habitat <- "Water"
water_mg <- generate_habitat_accumulation_data(otu_sample_meta_mg, habitat, 40)
water_mg$country <- "Madagascar"
water_se <- generate_habitat_accumulation_data(otu_sample_meta_se, habitat, 20)
water_se$country <- "Sweden"
water <- rbind(water_mg, water_se)

cat("Computing data for Wood\n")
habitat <- "Wood"
wood_mg <- generate_habitat_accumulation_data(otu_sample_meta_mg, habitat, 40)
wood_mg$country <- "Madagascar"
wood_se <- generate_habitat_accumulation_data(otu_sample_meta_se, habitat, 20)
wood_se$country <- "Sweden"
wood <- rbind(wood_mg, wood_se)

cat("Computing data for Temporary habitats\n")
habitat <- "Temporary habitats"
temporary_mg <- generate_habitat_accumulation_data(otu_sample_meta_mg, habitat, 40)
temporary_mg$country <- "Madagascar"
temporary_se <- generate_habitat_accumulation_data(otu_sample_meta_se, habitat, 20)
temporary_se$country <- "Sweden"
temporary <- rbind(temporary_mg, temporary_se)

cat("Computing data for Fungi\n")
habitat <- "Fungi"
fungi_mg <- generate_habitat_accumulation_data(otu_sample_meta_mg, habitat, 40)
fungi_mg$country <- "Madagascar"
fungi_se <- generate_habitat_accumulation_data(otu_sample_meta_se, habitat, 20)
fungi_se$country <- "Sweden"
fungi <- rbind(fungi_mg, fungi_se)

write.tsv <- function(D, file) { write.table(D, file, row.names=FALSE, sep="\t") }
write.tsv(plants,"../data/sample_acc_plants_prop.tsv")
write.tsv(soil,"../data/sample_acc_soil_prop.tsv")
write.tsv(water,"../data/sample_acc_water_prop.tsv")
write.tsv(wood,"../data/sample_acc_wood_prop.tsv")
write.tsv(temporary,"../data/sample_acc_temporary_prop.tsv")
write.tsv(fungi,"../data/sample_acc_fungi_prop.tsv")

