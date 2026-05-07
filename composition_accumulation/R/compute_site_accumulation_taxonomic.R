# Compute permuted site accumulation data for taxonomic groups

source("filter_se_traps.R")

# Generate reproducible random numbers
set.seed(10)


# Read in prepared data
# ---------------------

otu_site_meta_mg <- read.delim("../data/otu_site_meta_mg.tsv")
otu_site_meta_se <- read.delim("../data/otu_site_meta_se.tsv")
traps_to_keep <- filter_se_traps(20)    # Only keep se traps with 20 or more samples
otu_site_meta_se <- otu_site_meta_se[otu_site_meta_se$trapID %in% traps_to_keep,]


# Compute accumulation data for taxonomic groups
# ----------------------------------------------

generate_taxonomic_accumulation_data <- function(D, taxa) {
    res <- data.frame()
    sites <- unique(D$trapID)
    num_sites <- length(sites)
    for (n in 1:num_sites) {
        for (k in 1:100) {
            idx <- sample(1:num_sites,size=n)
            traps <- sites[idx]
            X <- D[D$trapID %in% traps,]
            taxon_otu_prop <- length(unique(X$cluster[X$Order %in% taxa])) / length(unique(X$cluster))
            res <- rbind(res,list(samples=n, taxon_otu_prop=taxon_otu_prop))
        }
    }
    return (res)
}

cat("Computing data for Diptera\n")
dipt_mg <- generate_taxonomic_accumulation_data(otu_site_meta_mg, "Diptera")
dipt_mg$country <- "Madagascar"
dipt_se <- generate_taxonomic_accumulation_data(otu_site_meta_se, "Diptera")
dipt_se$country <- "Sweden"
dipt <- rbind(dipt_mg, dipt_se)

cat("Computing data for Hymenoptera\n")
hyme_mg <- generate_taxonomic_accumulation_data(otu_site_meta_mg, "Hymenoptera")
hyme_mg$country <- "Madagascar"
hyme_se <- generate_taxonomic_accumulation_data(otu_site_meta_se, "Hymenoptera")
hyme_se$country <- "Sweden"
hyme <- rbind(hyme_mg, hyme_se)

cat("Computing data for Coleoptera\n")
cole_mg <- generate_taxonomic_accumulation_data(otu_site_meta_mg, "Coleoptera")
cole_mg$country <- "Madagascar"
cole_se <- generate_taxonomic_accumulation_data(otu_site_meta_se, "Coleoptera")
cole_se$country <- "Sweden"
cole <- rbind(cole_mg, cole_se)

cat("Computing data for Lepidoptera\n")
lepi_mg <- generate_taxonomic_accumulation_data(otu_site_meta_mg, "Lepidoptera")
lepi_mg$country <- "Madagascar"
lepi_se <- generate_taxonomic_accumulation_data(otu_site_meta_se, "Lepidoptera")
lepi_se$country <- "Sweden"
lepi <- rbind(lepi_mg, lepi_se)

cat("Computing data for Hemiptera\n")
hemi_mg <- generate_taxonomic_accumulation_data(otu_site_meta_mg, "Hemiptera")
hemi_mg$country <- "Madagascar"
hemi_se <- generate_taxonomic_accumulation_data(otu_site_meta_se, "Hemiptera")
hemi_se$country <- "Sweden"
hemi <- rbind(hemi_mg, hemi_se)

cat("Computing data for Other\n")
big_five <- c("Diptera","Hymenoptera","Coleoptera","Lepidoptera","Hemiptera")
orders <- unique(otu_site_meta_mg$Order)
other_mg <- generate_taxonomic_accumulation_data(otu_site_meta_mg, orders[!(orders %in% big_five)])
other_mg$country <- "Madagascar"
orders <- unique(otu_site_meta_se$Order)
other_se <- generate_taxonomic_accumulation_data(otu_site_meta_se, orders[!(orders %in% big_five)])
other_se$country <- "Sweden"
other <- rbind(other_mg, other_se)


# Write resulting data files
# --------------------------

write.tsv <- function(D, file) { write.table(D, file, row.names=FALSE, sep="\t") }
write.tsv(dipt,"../data/site_acc_diptera_prop.tsv")
write.tsv(hyme,"../data/site_acc_hymenoptera_prop.tsv")
write.tsv(cole,"../data/site_acc_coleoptera_prop.tsv")
write.tsv(lepi,"../data/site_acc_lepidoptera_prop.tsv")
write.tsv(hemi,"../data/site_acc_hemiptera_prop.tsv")
write.tsv(other,"../data/site_acc_other_orders_prop.tsv")

