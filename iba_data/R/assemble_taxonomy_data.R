# Assemble taxonomy data based on phylogenetic annotation
# Collect some filtering and annotation statistics in the process

library(data.table)


# TODO: Change to your local path to the FigShare processed
# data and metadata repos
path <- "~/dev/figshare-repos/iba/processed_data/v4-prel/"
metapath <- "~/dev/figshare-repos/iba/raw_data/v6/"

# Make a subset for OTUs in SE Malaise and litter samples for filter stats
C <- fread(paste0(path,"cluster_counts_SE.tsv"))
M <- read.delim(paste0(metapath,"CO1_sequencing_metadata_SE.tsv"))
M <- M[M$dataset %in% c("CO1_lysate_2019_SE","CO1_soil_litter_2019_SE"),]
M <- M[M$lab_sample_type=="sample" & M$sequencing_successful == TRUE,]
M <- M[!grepl("_S",M$sampleID_FIELD),]
idx <- which(colnames(C) %in% M$sampleID_NGI)
idx <- c(1,idx)
C <- C[,..idx]
C <- C[rowSums(C[,-1])>0,]
se_clusters <- C$cluster

F1 <- data.frame(list(step=character(),clusters=numeric()))
F2 <- data.frame(list(step=character(),clusters=numeric()))
T1 <- read.delim(paste0(path,"cluster_taxonomy_MG.tsv"))
T2 <- read.delim(paste0(path,"cluster_taxonomy_SE.tsv"))
F1 <- rbind(F1,list(step="after_chimera_removal",clusters=sum(T1$representative==1)))
F2 <- rbind(F2,list(step="after_chimera_removal",clusters=sum(T2$representative==1 & T2$cluster %in% se_clusters)))
T1 <- read.delim(paste0(path,"noise_filtered_cluster_taxonomy_MG.tsv"))
T2 <- read.delim(paste0(path,"noise_filtered_cluster_taxonomy_SE.tsv"))
F1 <- rbind(F1,list(step="after_neeat_filter",clusters=sum(T1$representative==1)))
F2 <- rbind(F2,list(step="after_neeat_filter",clusters=sum(T2$representative==1 & T2$cluster %in% se_clusters)))
T1 <- read.delim(paste0(path,"cleaned_noise_filtered_cluster_taxonomy_MG.tsv"))
T2 <- read.delim(paste0(path,"cleaned_noise_filtered_cluster_taxonomy_SE.tsv"))
T1 <- T1[T1$representative==1,] # should be redundant...
T2 <- T2[T2$representative==1,] # should be redundant...
rename_tax <- function(D, suffix) {
    ranks <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species","BOLD_bin")
    idx <- which(colnames(D) %in% ranks)
    new_names <- paste0(colnames(D)[idx],suffix)
    colnames(D)[idx] <- new_names
    return (D)
}
T1 <- rename_tax(T1,".default")
T2 <- rename_tax(T2,".default")
F1 <- rbind(F1,list(step="after_clean_filter",clusters=nrow(T1)))
F2 <- rbind(F2,list(step="after_clean_filter",clusters=nrow(T2[T2$cluster %in% se_clusters,])))

# Phylogenetic annotation
E1 <- read.delim(paste0(path,"asv_taxonomy_epang_MG.tsv"))
E2 <- read.delim(paste0(path,"asv_taxonomy_epang_SE.tsv"))

# Sintax annotation
S1 <- read.delim(paste0(path,"asv_taxonomy_sintax_MG.tsv"))
S2 <- read.delim(paste0(path,"asv_taxonomy_sintax_SE.tsv"))
S1 <- rename_tax(S1,".sintax")
S2 <- rename_tax(S2,".sintax")

# Vsearch annotation
V1 <- read.delim(paste0(path,"asv_taxonomy_vsearch_MG.tsv"))
V2 <- read.delim(paste0(path,"asv_taxonomy_vsearch_SE.tsv"))
V1 <- rename_tax(V1,".vsearch")
V2 <- rename_tax(V2,".vsearch")

# Read in and remove spikeins
spike1 <- read.delim(paste0(path,"spikeins_tax_MG.tsv"))
spike2 <- read.delim(paste0(path,"spikeins_tax_SE.tsv"))
T1 <- T1[!(T1$cluster %in% spike1$cluster),]
T2 <- T2[!(T2$cluster %in% spike2$cluster),]
F1 <- rbind(F1,list(step="after_spike_removal",clusters=nrow(T1)))
F2 <- rbind(F2,list(step="after_spike_removal",clusters=nrow(T2[T2$cluster %in% se_clusters,])))

# Merge data frames
D1 <- merge(T1,E1,by="ASV")
D2 <- merge(T2,E2,by="ASV")
D1 <- merge(D1,S1,by="ASV")
D2 <- merge(D2,S2,by="ASV")
D1 <- merge(D1,V1,by="ASV")
D2 <- merge(D2,V2,by="ASV")

# Filter out hexapods
hexapods <- c("Collembola","Diplura","Protura","Insecta")
D1 <- D1[D1$Class %in% hexapods,]
D2 <- D2[D2$Class %in% hexapods,]
F1 <- rbind(F1,list(step="hexapods",clusters=nrow(D1)))
F2 <- rbind(F2,list(step="hexapods",clusters=nrow(D2[D2$cluster %in% se_clusters,])))


# Filter out likely errors based on a combination of criteria
# -----------------------------------------------------------

# 1. Remove poor-quality phylogenetic placements
D1 <- D1[D1$aLWR > 0.5,]
D2 <- D2[D2$aLWR > 0.5,]
F1 <- rbind(F1,list(step="placement_robust",clusters=nrow(D1)))
F2 <- rbind(F2,list(step="placement_robust",clusters=nrow(D2[D2$cluster %in% se_clusters,])))

# 2. Remove if sintax has a solid match to a BOLD bin in a different order (family-level classification is too uncertain to use here)
D1 <- D1[!(grepl("BOLD",D1$BOLD_bin.sintax) & D1$Order != D1$Order.sintax),]
D2 <- D2[!(grepl("BOLD",D2$BOLD_bin.sintax) & D2$Order != D2$Order.sintax),]
F1 <- rbind(F1,list(step="sintax_BOLD_conflict_removed",clusters=nrow(D1)))
F2 <- rbind(F2,list(step="sintax_BOLD_conflict_removed",clusters=nrow(D2[D2$cluster %in% se_clusters,])))

# 3. Remove if vsearch is certain that this is the wrong phylum
D1 <- D1[D1$Phylum==D1$Phylum.vsearch,]
D2 <- D2[D2$Phylum==D2$Phylum.vsearch,]
F1 <- rbind(F1,list(step="vsearch_phylum_conflict_removed",clusters=nrow(D1)))
F2 <- rbind(F2,list(step="vsearch_phylum_conflict_removed",clusters=nrow(D2[D2$cluster %in% se_clusters,])))

# Only keep clusters assigned to clade
D1 <- D1[!grepl("unclassified",D1$Clade),]
D2 <- D2[!grepl("unclassified",D2$Clade),]
F1 <- rbind(F1,list(step="assigned_to_clade",clusters=nrow(D1)))
F2 <- rbind(F2,list(step="assigned_to_clade",clusters=nrow(D2[D2$cluster %in% se_clusters,])))

# Write final cluster taxonomy data
saveRDS(D1,"../cluster_taxonomy_mg.rds")
saveRDS(D2,"../cluster_taxonomy_se.rds")

# Write filtering stats
F1$remaining_frac <- F1$clusters / F1$clusters[1]
F2$remaining_frac <- F2$clusters / F2$clusters[1]
write.tsv(F1,"../filter_annotation_stats_mg.tsv")
write.tsv(F2,"../filter_annotation_stats_lysate_litter_se.tsv")

