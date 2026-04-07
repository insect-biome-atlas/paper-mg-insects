# Assemble taxonomy data based on phylogenetic annotation

# TODO: Change to your local path to the FigShare processed data repo
path <- "~/dev/figshare-repos/iba/processed_data/v4-prel/"
T1 <- read.delim(paste0(path,"cleaned_noise_filtered_cluster_taxonomy_MG.tsv"))
T2 <- read.delim(paste0(path,"cleaned_noise_filtered_cluster_taxonomy_SE.tsv"))
T1 <- T1[T1$representative==1,]
T2 <- T2[T2$representative==1,]
rename_tax <- function(D, suffix) {
    ranks <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species","BOLD_bin")
    idx <- which(colnames(D) %in% ranks)
    new_names <- paste0(colnames(D)[idx],suffix)
    colnames(D)[idx] <- new_names
    return (D)
}
T1 <- rename_tax(T1,".default")
T2 <- rename_tax(T2,".default")

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

# Filter out likely errors based on a combination of criteria

# 1. Remove poor-quality phylogenetic placements
D1 <- D1[D1$aLWR > 0.5,]
D2 <- D2[D2$aLWR > 0.5,]

# 2. Remove if sintax has a solid match to a BOLD bin in a different order (family-level classification is too uncertain to use here)
D1 <- D1[!(grepl("BOLD",D1$BOLD_bin.sintax) & D1$Order != D1$Order.sintax),]
D2 <- D2[!(grepl("BOLD",D2$BOLD_bin.sintax) & D2$Order != D2$Order.sintax),]

# 3. Remove if vsearch is certain that this is the wrong phylum
D1 <- D1[D1$Phylum==D1$Phylum.vsearch,]
D2 <- D2[D2$Phylum==D2$Phylum.vsearch,]

# Only keep clusters assigned to clade
D1 <- D1[!grepl("unclassified",D1$Clade),]
D2 <- D2[!grepl("unclassified",D2$Clade),]

# Write final cluster taxonomy data
write.tsv <- function(D,file) { write.table(D,file,row.names=FALSE,sep="\t") }
write.tsv(D1,"cluster_taxonomy_mg.tsv")
write.tsv(D2,"cluster_taxonomy_se.tsv")

