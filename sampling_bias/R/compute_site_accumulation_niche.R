# Compute permuted site accumulation data for niches
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
sapro <- generate_niche_accumulation_data(otu_site_meta_se, niche, niche_group)

cat("Computing data for phytophages\n")
niche <- "Phytophage"
niche_group <- c("Saprophage","Phytophage","Predator")
phyto <- generate_niche_accumulation_data(otu_site_meta_se, niche, niche_group)

cat("Computing data for predators\n")
niche <- "Predator"
niche_group <- c("Saprophage","Phytophage","Predator")
pred <- generate_niche_accumulation_data(otu_site_meta_se, niche, niche_group)


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
sapro_ph <- generate_ph_accumulation_data(otu_site_meta_se, parasitoid, host)

cat("Computing data for phytophage parasitoids\n")
parasitoid <- "Phytophage-parasitoid"
host <- "Phytophage"
phyto_ph <- generate_ph_accumulation_data(otu_site_meta_se, parasitoid, host)

cat("Computing data for predator parasitoids\n")
parasitoid <- "Predator-parasitoid"
host <- "Predator"
pred_ph <- generate_ph_accumulation_data(otu_site_meta_se, parasitoid, host)


# Write resulting data files
# ------------------------------

write.tsv <- function(D, file) { write.table(D, file, row.names=FALSE, sep="\t") }
write.tsv(sapro,"../data/site_acc_saprophage_prop.tsv")
write.tsv(phyto,"../data/site_acc_phytophage_prop.tsv")
write.tsv(pred,"../data/site_acc_predator_prop.tsv")
write.tsv(sapro_ph,"../data/site_acc_saprophage_ph_ratio.tsv")
write.tsv(phyto_ph,"../data/site_acc_phytophage_ph_ratio.tsv")
write.tsv(pred_ph,"../data/site_acc_predator_ph_ratio.tsv")

