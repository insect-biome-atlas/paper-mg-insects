# Generate OTU tables expected by the vegan package

# Read in abundance data
malaise_mg <- read.delim("../cluster_counts_spec_malaise_mg.tsv", row.names=1)
malaise_se <- read.delim("../cluster_counts_spec_malaise_se.tsv", row.names=1)

# Read in occcurrence data
malaise_litter_mg <- read.delim("../cluster_occurrence_malaise_litter_mg.tsv", row.names=1)
malaise_litter_se <- read.delim("../cluster_occurrence_malaise_litter_se.tsv", row.names=1)

# Read in sample metadata
malaise_sample_meta_mg <- read.delim("../malaise_sample_meta_mg.tsv")
malaise_sample_meta_se <- read.delim("../malaise_sample_meta_se.tsv")
malaise_litter_sample_meta_mg <- read.delim("../malaise_litter_sample_meta_mg.tsv")
malaise_litter_sample_meta_se <- read.delim("../malaise_litter_sample_meta_se.tsv")

# Only keep samples with complete metadata (using the fact that the above files only contain complete cases)
idx1 <- which(colnames(malaise_mg) %in% malaise_sample_meta_mg$sampleID_NGI)
malaise_mg <- malaise_mg[,idx1]
idx2 <- which(colnames(malaise_se) %in% malaise_sample_meta_se$sampleID_NGI)
malaise_se <- malaise_se[,idx2]
idx1 <- which(colnames(malaise_litter_mg) %in% malaise_litter_sample_meta_mg$sampleID_NGI)
malaise_litter_mg <- malaise_litter_mg[,idx1]
idx2 <- which(colnames(malaise_litter_se) %in% malaise_litter_sample_meta_se$sampleID_NGI)
malaise_litter_se <- malaise_litter_se[,idx2]

# Keep complete data from Sweden as a separate dataset
malaise_litter_complete_se <- malaise_litter_se # Save complete data from Sweden in a separate dataset

# Only keep Swedish forest data samples in main dataset
idx <- which(colnames(malaise_se) %in% malaise_sample_meta_se$sampleID_NGI[malaise_sample_meta_se$trap_habitat=="forest"])
malaise_litter_complete_se <- malaise_litter_se # Save complete data from Sweden in a separate dataset
malaise_se <- malaise_se[,idx]
malaise_litter_se <- malaise_litter_se[,idx]

# For abundance data, only keep Swedish samples from traps with at least 20 samples
X <- unique(malaise_sample_meta_se[,c("trapID","sampleID_NGI")])
trap_samples <- table(X$trapID)
traps_to_keep <- names(trap_samples[trap_samples >= 20])
samples_to_keep <- X$sampleID_NGI[X$trapID %in% traps_to_keep]
idx <- which(colnames(malaise_se) %in% samples_to_keep)
malaise_se <- malaise_se[,idx]

# For occurrence data, only keep Swedish samples from traps with at least 10 samples
traps_to_keep <- names(trap_samples[trap_samples >= 10])
samples_to_keep <- X$sampleID_NGI[X$trapID %in% traps_to_keep]
idx <- which(colnames(malaise_litter_se) %in% samples_to_keep)
malaise_litter_se <- malaise_litter_se[,idx]
idx <- which(colnames(malaise_litter_complete_se) %in% samples_to_keep)
malaise_litter_complete_se <- malaise_litter_complete_se[,idx]

# Remove empty rows
malaise_mg <- malaise_mg[rowSums(malaise_mg)>0,]
malaise_se <- malaise_se[rowSums(malaise_se)>0,]
malaise_litter_mg <- malaise_litter_mg[rowSums(malaise_litter_mg)>0,]
malaise_litter_se <- malaise_litter_se[rowSums(malaise_litter_se)>0,]
malaise_litter_complete_se <- malaise_litter_complete_se[rowSums(malaise_litter_complete_se)>0,]

# Transpose tables so that rows are samples
sample_otu_malaise_mg <- data.frame(t(malaise_mg))
sample_otu_malaise_se <- data.frame(t(malaise_se))
sample_otu_occurrence_combined_mg <- data.frame(t(malaise_litter_mg))
sample_otu_occurrence_combined_se <- data.frame(t(malaise_litter_se))
sample_otu_occurrence_combined_complete_se <- data.frame(t(malaise_litter_complete_se))

# Save OTU sample abundance tables for vegan
saveRDS(sample_otu_malaise_mg,"../sample_otu_abundance_malaise_mg.rds")
saveRDS(sample_otu_malaise_se,"../sample_otu_abundance_malaise_se.rds")

# Summarize counts for sites
make_col_rownames <- function(D, col) {
    idx <- which(colnames(D)==col)
    rownames(D) <- D[,idx]
    D[,-idx]
}
sites_mg <- malaise_sample_meta_mg$trapID[match(rownames(sample_otu_malaise_mg),malaise_sample_meta_mg$sampleID_NGI)]
sites_se <- malaise_sample_meta_se$trapID[match(rownames(sample_otu_malaise_se),malaise_sample_meta_se$sampleID_NGI)]
site_otu_malaise_mg <- aggregate(sample_otu_malaise_mg, by=list(trapID=sites_mg), FUN=sum)
site_otu_malaise_se <- aggregate(sample_otu_malaise_se, by=list(trapID=sites_se), FUN=sum)
site_otu_malaise_mg <- make_col_rownames(site_otu_malaise_mg,"trapID")
site_otu_malaise_se <- make_col_rownames(site_otu_malaise_se,"trapID")

# Save OTU site abundance tables for vegan
saveRDS(site_otu_malaise_mg,"../site_otu_abundance_malaise_mg.tsv")
saveRDS(site_otu_malaise_se,"../site_otu_abundance_malaise_se.tsv")

# Summarize occurrence for sites
sites_mg <- malaise_litter_sample_meta_mg$trapID[match(rownames(sample_otu_occurrence_combined_mg),malaise_litter_sample_meta_mg$sampleID_NGI)]
sites_se <- malaise_litter_sample_meta_se$trapID[match(rownames(sample_otu_occurrence_combined_se),malaise_litter_sample_meta_se$sampleID_NGI)]
site_otu_occurrence_combined_mg <- aggregate(sample_otu_occurrence_combined_mg, by=list(trapID=sites_mg), FUN=function(x){1*(sum(x)!=0)})
site_otu_occurrence_combined_se <- aggregate(sample_otu_occurrence_combined_se, by=list(trapID=sites_se), FUN=function(x){1*(sum(x)!=0)})
site_otu_occurrence_combined_mg <- make_col_rownames(site_otu_occurrence_combined_mg,"trapID")
site_otu_occurrence_combined_se <- make_col_rownames(site_otu_occurrence_combined_se,"trapID")
site_otu_occurrence_malaise_mg <- data.frame(1*(site_otu_malaise_mg > 0))
site_otu_occurrence_malaise_se <- data.frame(1*(site_otu_malaise_se > 0))

sites_se <- malaise_sample_meta_se$trapID[match(rownames(sample_otu_occurrence_combined_complete_se),malaise_litter_sample_meta_se$sampleID_NGI)]
site_otu_occurrence_combined_complete_se <- aggregate(sample_otu_occurrence_combined_complete_se, by=list(trapID=sites_se), FUN=function(x){sum(x)!=0})
site_otu_occurrence_combined_complete_se <- make_col_rownames(site_otu_occurrence_combined_complete_se,"trapID")

# Save OTU site occurrence tables for vegan
saveRDS(site_otu_occurrence_combined_mg,"../site_otu_occurrence_combined_mg.rds")
saveRDS(site_otu_occurrence_combined_se,"../site_otu_occurrence_combined_se.rds")
saveRDS(site_otu_occurrence_malaise_mg,"../site_otu_occurrence_malaise_mg.rds")
saveRDS(site_otu_occurrence_malaise_se,"../site_otu_occurrence_malaise_se.rds")
saveRDS(site_otu_occurrence_combined_complete_se,"../site_otu_occurrence_combined_complete_se.rds")
