# Compute permuted sample accumulation data for taxon groups

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


# Compute accumulation data for taxonomic groups
# ----------------------------------------------

generate_taxonomic_accumulation_data <- function(D, taxa, num_samples) {
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
            taxon_otu_prop <- length(unique(X$cluster[X$Order %in% taxa])) / length(unique(X$cluster))
            res <- rbind(res,list(samples=n, taxon_otu_prop=taxon_otu_prop))
        }
    }
    cat("\n")
    return (res)
}

cat("Computing data for Diptera\n")
dipt_mg <- generate_taxonomic_accumulation_data(otu_sample_meta_mg, "Diptera", 40)
dipt_mg$country <- "Madagascar"
dipt_se <- generate_taxonomic_accumulation_data(otu_sample_meta_se, "Diptera", 20)
dipt_se$country <- "Sweden"
dipt <- rbind(dipt_mg, dipt_se)

cat("Computing data for Hymenoptera\n")
hyme_mg <- generate_taxonomic_accumulation_data(otu_sample_meta_mg, "Hymenoptera", 40)
hyme_mg$country <- "Madagascar"
hyme_se <- generate_taxonomic_accumulation_data(otu_sample_meta_se, "Hymenoptera", 20)
hyme_se$country <- "Sweden"
hyme <- rbind(hyme_mg, hyme_se)

cat("Computing data for Coleoptera\n")
cole_mg <- generate_taxonomic_accumulation_data(otu_sample_meta_mg, "Coleoptera", 40)
cole_mg$country <- "Madagascar"
cole_se <- generate_taxonomic_accumulation_data(otu_sample_meta_se, "Coleoptera", 20)
cole_se$country <- "Sweden"
cole <- rbind(cole_mg, cole_se)

cat("Computing data for Lepidoptera\n")
lepi_mg <- generate_taxonomic_accumulation_data(otu_sample_meta_mg, "Lepidoptera", 40)
lepi_mg$country <- "Madagascar"
lepi_se <- generate_taxonomic_accumulation_data(otu_sample_meta_se, "Lepidoptera", 20)
lepi_se$country <- "Sweden"
lepi <- rbind(lepi_mg, lepi_se)

cat("Computing data for Hemiptera\n")
hemi_mg <- generate_taxonomic_accumulation_data(otu_sample_meta_mg, "Hemiptera", 40)
hemi_mg$country <- "Madagascar"
hemi_se <- generate_taxonomic_accumulation_data(otu_sample_meta_se, "Hemiptera", 20)
hemi_se$country <- "Sweden"
hemi <- rbind(hemi_mg, hemi_se)

cat("Computing data for Other\n")
big_five <- c("Diptera","Hymenoptera","Coleoptera","Lepidoptera","Hemiptera")
orders <- unique(otu_sample_meta_mg$Order)
other_mg <- generate_taxonomic_accumulation_data(otu_sample_meta_mg, orders[!(orders %in% big_five)], 40)
other_mg$country <- "Madagascar"
orders <- unique(otu_sample_meta_se$Order)
other_se <- generate_taxonomic_accumulation_data(otu_sample_meta_se, orders[!(orders %in% big_five)], 20)
other_se$country <- "Sweden"
other <- rbind(other_mg, other_se)


# Write resulting data files
# --------------------------

write.tsv <- function(D, file) { write.table(D, file, row.names=FALSE, sep="\t") }
write.tsv(dipt,"../data/sample_acc_diptera_prop.tsv")
write.tsv(hyme,"../data/sample_acc_hymenoptera_prop.tsv")
write.tsv(cole,"../data/sample_acc_coleoptera_prop.tsv")
write.tsv(lepi,"../data/sample_acc_lepidoptera_prop.tsv")
write.tsv(hemi,"../data/sample_acc_hemiptera_prop.tsv")
write.tsv(other,"../data/sample_acc_other_orders_prop.tsv")

