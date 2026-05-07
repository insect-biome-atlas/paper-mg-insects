# Compute permuted sample accumulation data for niches

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


# Compute accumulation data for major niches
# --------------------------------------------

generate_niche_accumulation_data <- function(D, niche, niche_group, num_samples) {
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
            niche_otu_prop <- length(unique(X$cluster[X$Niche==niche])) / length(unique(X$cluster[X$Niche %in% niche_group]))
            res <- rbind(res,list(samples=n, niche_otu_prop=niche_otu_prop))
        }
    }
    cat("\n")
    return (res)
}

cat("Computing data for saprophages\n")
niche <- "Saprophage"
niche_group <- c("Saprophage","Phytophage","Predator")
sapro_mg <- generate_niche_accumulation_data(otu_sample_meta_mg, niche, niche_group, 40)
sapro_mg$country <- "Madagascar"
sapro_se <- generate_niche_accumulation_data(otu_sample_meta_se, niche, niche_group, 20)
sapro_se$country <- "Sweden"
sapro <- rbind(sapro_mg, sapro_se)

cat("Computing data for phytophages\n")
niche <- "Phytophage"
niche_group <- c("Saprophage","Phytophage","Predator")
phyto_mg <- generate_niche_accumulation_data(otu_sample_meta_mg, niche, niche_group, 40)
phyto_mg$country <- "Madagascar"
phyto_se <- generate_niche_accumulation_data(otu_sample_meta_se, niche, niche_group, 20)
phyto_se$country <- "Sweden"
phyto <- rbind(phyto_mg, phyto_se)

cat("Computing data for predators\n")
niche <- "Predator"
niche_group <- c("Saprophage","Phytophage","Predator")
pred_mg <- generate_niche_accumulation_data(otu_sample_meta_mg, niche, niche_group, 40)
pred_mg$country <- "Madagascar"
pred_se <- generate_niche_accumulation_data(otu_sample_meta_se, niche, niche_group, 20)
pred_se$country <- "Sweden"
pred <- rbind(pred_mg, pred_se)

write.tsv <- function(D, file) { write.table(D, file, row.names=FALSE, sep="\t") }
write.tsv(sapro,"../data/sample_acc_saprophage_prop.tsv")
write.tsv(phyto,"../data/sample_acc_phytophage_prop.tsv")
write.tsv(pred,"../data/sample_acc_predator_prop.tsv")


# Compute accumulation data for parasitoid:host ratios
# ----------------------------------------------------

generate_ph_accumulation_data <- function(D, parasitoid, host, num_samples) {
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
            ph_otu_ratio <- length(unique(X$cluster[X$Niche==parasitoid])) / length(unique(X$cluster[X$Niche==host]))
            res <- rbind(res,list(samples=n, ph_otu_ratio=ph_otu_ratio))
        }
    }
    cat("\n")
    return (res)
}

cat("Computing data for saprophage parasitoids\n")
parasitoid <- "Saprophage-parasitoid"
host <- "Saprophage"
sapro_ph_mg <- generate_ph_accumulation_data(otu_sample_meta_mg, parasitoid, host, 40)
sapro_ph_mg$country <- "Madagascar"
sapro_ph_se <- generate_ph_accumulation_data(otu_sample_meta_se, parasitoid, host, 20)
sapro_ph_se$country <- "Sweden"
sapro_ph <- rbind(sapro_ph_mg, sapro_ph_se)

cat("Computing data for phytophage parasitoids\n")
parasitoid <- "Phytophage-parasitoid"
host <- "Phytophage"
phyto_ph_mg <- generate_ph_accumulation_data(otu_sample_meta_mg, parasitoid, host, 40)
phyto_ph_mg$country <- "Madagascar"
phyto_ph_se <- generate_ph_accumulation_data(otu_sample_meta_se, parasitoid, host, 20)
phyto_ph_se$country <- "Sweden"
phyto_ph <- rbind(phyto_ph_mg, phyto_ph_se)

cat("Computing data for predator parasitoids\n")
parasitoid <- "Predator-parasitoid"
host <- "Predator"
pred_ph_mg <- generate_ph_accumulation_data(otu_sample_meta_mg, parasitoid, host, 40)
pred_ph_mg$country <- "Madagascar"
pred_ph_se <- generate_ph_accumulation_data(otu_sample_meta_se, parasitoid, host, 20)
pred_ph_se$country <- "Sweden"
pred_ph <- rbind(pred_ph_mg, pred_ph_se)

write.tsv(sapro_ph,"../data/sample_acc_saprophage_ph_ratio.tsv")
write.tsv(phyto_ph,"../data/sample_acc_phytophage_ph_ratio.tsv")
write.tsv(pred_ph,"../data/sample_acc_predator_ph_ratio.tsv")

