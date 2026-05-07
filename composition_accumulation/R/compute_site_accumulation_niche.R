# Compute permuted site accumulation data for niches

source("filter_se_traps.R")

# Generate reproducible random numbers
set.seed(10)


# Read in prepared data
# ---------------------

otu_site_meta_mg <- read.delim("../data/otu_site_meta_mg.tsv")
otu_site_meta_se <- read.delim("../data/otu_site_meta_se.tsv")
traps_to_keep <- filter_se_traps(20)    # Only keep se traps with 20 or more samples
otu_site_meta_se <- otu_site_meta_se[otu_site_meta_se$trapID %in% traps_to_keep,]


# Compute accumulation data for major niches
# --------------------------------------------

generate_niche_accumulation_data <- function(D, niche, niche_group) {
    res <- data.frame()
    sites <- unique(D$trapID)
    num_sites <- length(sites)
    for (n in 1:num_sites) {
        for (k in 1:100) {
            idx <- sample(1:num_sites,size=n)
            traps <- sites[idx]
            X <- D[D$trapID %in% traps,]
            niche_otu_prop <- length(unique(X$cluster[X$Niche==niche])) / length(unique(X$cluster[X$Niche %in% niche_group]))
            res <- rbind(res,list(samples=n, niche_otu_prop=niche_otu_prop))
        }
    }
    return (res)
}

cat("Computing data for saprophages\n")
niche <- "Saprophage"
niche_group <- c("Saprophage","Phytophage","Predator")
sapro_mg <- generate_niche_accumulation_data(otu_site_meta_mg, niche, niche_group)
sapro_mg$country <- "Madagascar"
sapro_se <- generate_niche_accumulation_data(otu_site_meta_se, niche, niche_group)
sapro_se$country <- "Sweden"
sapro <- rbind(sapro_mg, sapro_se)

cat("Computing data for phytophages\n")
niche <- "Phytophage"
niche_group <- c("Saprophage","Phytophage","Predator")
phyto_mg <- generate_niche_accumulation_data(otu_site_meta_mg, niche, niche_group)
phyto_mg$country <- "Madagascar"
phyto_se <- generate_niche_accumulation_data(otu_site_meta_se, niche, niche_group)
phyto_se$country <- "Sweden"
phyto <- rbind(phyto_mg, phyto_se)

cat("Computing data for predators\n")
niche <- "Predator"
niche_group <- c("Saprophage","Phytophage","Predator")
pred_mg <- generate_niche_accumulation_data(otu_site_meta_mg, niche, niche_group)
pred_mg$country <- "Madagascar"
pred_se <- generate_niche_accumulation_data(otu_site_meta_se, niche, niche_group)
pred_se$country <- "Sweden"
pred <- rbind(pred_mg, pred_se)

write.tsv <- function(D, file) { write.table(D, file, row.names=FALSE, sep="\t") }
write.tsv(sapro,"../data/site_acc_saprophage_prop.tsv")
write.tsv(phyto,"../data/site_acc_phytophage_prop.tsv")
write.tsv(pred,"../data/site_acc_predator_prop.tsv")


# Compute accumulation data for parasitoid:host ratios
# ----------------------------------------------------

generate_ph_accumulation_data <- function(D, parasitoid, host) {
    res <- data.frame()
    sites <- unique(D$trapID)
    num_sites <- length(sites)
    for (n in 1:num_sites) {
        for (k in 1:100) {
            idx <- sample(1:num_sites,size=n)
            traps <- sites[idx]
            X <- D[D$trapID %in% traps,]
            ph_otu_ratio <- length(unique(X$cluster[X$Niche==parasitoid])) / length(unique(X$cluster[X$Niche==host]))
            res <- rbind(res,list(samples=n, ph_otu_ratio=ph_otu_ratio))
        }
    }
    return (res)
}

cat("Computing data for saprophage parasitoids\n")
parasitoid <- "Saprophage-parasitoid"
host <- "Saprophage"
sapro_ph_mg <- generate_ph_accumulation_data(otu_site_meta_mg, parasitoid, host)
sapro_ph_mg$country <- "Madagascar"
sapro_ph_se <- generate_ph_accumulation_data(otu_site_meta_se, parasitoid, host)
sapro_ph_se$country <- "Sweden"
sapro_ph <- rbind(sapro_ph_mg, sapro_ph_se)

cat("Computing data for phytophage parasitoids\n")
parasitoid <- "Phytophage-parasitoid"
host <- "Phytophage"
phyto_ph_mg <- generate_ph_accumulation_data(otu_site_meta_mg, parasitoid, host)
phyto_ph_mg$country <- "Madagascar"
phyto_ph_se <- generate_ph_accumulation_data(otu_site_meta_se, parasitoid, host)
phyto_ph_se$country <- "Sweden"
phyto_ph <- rbind(phyto_ph_mg, phyto_ph_se)

cat("Computing data for predator parasitoids\n")
parasitoid <- "Predator-parasitoid"
host <- "Predator"
pred_ph_mg <- generate_ph_accumulation_data(otu_site_meta_mg, parasitoid, host)
pred_ph_mg$country <- "Madagascar"
pred_ph_se <- generate_ph_accumulation_data(otu_site_meta_se, parasitoid, host)
pred_ph_se$country <- "Sweden"
pred_ph <- rbind(pred_ph_mg, pred_ph_se)

write.tsv(sapro_ph,"../data/site_acc_saprophage_ph_ratio.tsv")
write.tsv(phyto_ph,"../data/site_acc_phytophage_ph_ratio.tsv")
write.tsv(pred_ph,"../data/site_acc_predator_ph_ratio.tsv")

