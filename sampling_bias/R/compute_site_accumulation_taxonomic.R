# Compute permuted site accumulation data for taxonomic groups
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
dipt <- generate_taxonomic_accumulation_data(otu_site_meta_se, "Diptera")

cat("Computing data for Hymenoptera\n")
hyme <- generate_taxonomic_accumulation_data(otu_site_meta_se, "Hymenoptera")

cat("Computing data for Coleoptera\n")
cole <- generate_taxonomic_accumulation_data(otu_site_meta_se, "Coleoptera")

cat("Computing data for Lepidoptera\n")
lepi <- generate_taxonomic_accumulation_data(otu_site_meta_se, "Lepidoptera")

cat("Computing data for Hemiptera\n")
hemi <- generate_taxonomic_accumulation_data(otu_site_meta_se, "Hemiptera")

cat("Computing data for Other\n")
big_five <- c("Diptera","Hymenoptera","Coleoptera","Lepidoptera","Hemiptera")
orders <- unique(otu_site_meta_se$Order)
other <- generate_taxonomic_accumulation_data(otu_site_meta_se, orders[!(orders %in% big_five)])


# Write resulting data files
# --------------------------

write.tsv <- function(D, file) { write.table(D, file, row.names=FALSE, sep="\t") }
write.tsv(dipt,"../data/site_acc_diptera_prop.tsv")
write.tsv(hyme,"../data/site_acc_hymenoptera_prop.tsv")
write.tsv(cole,"../data/site_acc_coleoptera_prop.tsv")
write.tsv(lepi,"../data/site_acc_lepidoptera_prop.tsv")
write.tsv(hemi,"../data/site_acc_hemiptera_prop.tsv")
write.tsv(other,"../data/site_acc_other_orders_prop.tsv")

