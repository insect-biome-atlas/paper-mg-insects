# Analyze different filtering criteria on problemtaic clusters
library(data.table)

# Read in metadata
cat("Reading in metadata\n")
M1 <- read.delim("~/dev/figshare-repos/iba/raw_data/v6/CO1_sequencing_metadata_SE.tsv")
M2 <- read.delim("~/dev/figshare-repos/iba/raw_data/v6/CO1_sequencing_metadata_MG.tsv")
T1 <- read.delim("~/dev/figshare-repos/iba/processed_data/v3/cluster_taxonomy_SE.tsv")
T2 <- read.delim("~/dev/figshare-repos/iba/processed_data/v3/cluster_taxonomy_MG.tsv")
R2 <- T2[T2$representative==1,]
E2 <- read.delim("~/dev/figshare-repos/iba/processed_data/v3/asv_taxonomy_sintax_MG.tsv")
RE2 <- E2[E2$ASV %in% R2$ASV,]

# Read in cluster counts
cat("Reading in cluster counts\n")
C1 <- read.delim("~/dev/figshare-repos/iba/processed_data/v3/cluster_counts_SE.tsv")
C2 <- read.delim("~/dev/figshare-repos/iba/processed_data/v3/cluster_counts_MG.tsv")
cat("Completed reading in cluster counts\n\n")

# Get data subsets

# Define function for adding total reads and no. samples
add_reads_samples <- function(D) {
  reads <- rowSums(D[,-1])
  samples <- rowSums(D[,-1]!=0)
  D$reads <- reads
  D$samples <- samples
  D
}

# Define function for getting data for specific samples
get_sample_data <- function(D,samples) {
  D <- D[,c(1,which(colnames(D) %in% samples))]
  D <- add_reads_samples(D)
  D[D$samples>0,]
}

# Define function for getting specific samples
get_samples <- function(D, sample_type, dataset) {
  
  D$sampleID_NGI[(D$lab_sample_type %in% sample_type) & D$sequencing_successful & (D$dataset==dataset)]
}

# Get mg litter sample data
foo <- function(st) { get_sample_data(C2, get_samples(M2, st, "CO1_litter_2019_MG")) }
mg_litter <- foo("sample")
mg_litter_pos_con <- foo("positive_control")
mg_litter_control <- foo(c("buffer_blank", "extraction_neg", "pcr_neg"))

# Get mg malaise sample data
foo <- function(st) { get_sample_data(C2, get_samples(M2, st, "CO1_lysate_2019_MG")) }
mg_lysate <- foo("sample")
mg_lysate_pos_con <- foo("positive_control")
mg_lysate_control <- foo(c("buffer_blank", "extraction_neg", "pcr_neg"))


# Find global set of potentially problematic clusters
# Define problem as MG sample clusters that share an ASV with Swedish ASVs
problem_clusters <- unique(T2$cluster[T2$ASV %in% T1$ASV])

# Restrict to hexapods
hexapods <- c("Collembola","Diplura","Protura","Insecta")
problem_clusters <- unique(T2$cluster[T2$ASV %in% T1$ASV & T2$Class %in% hexapods & !(grepl("unclassified.Insecta",T2$Order))])

evaluate_det <- function(D) {
  D <- D[D$cluster %in% problem_clusters,]
  cat("Clusters:", nrow(D), "\n")
  cat("Samples:\n")
  print(summary(D$samples))
  cat("Prop. samples:\n")
  print(summary(D$samples/(ncol(D)-3))) 
  cat("Reads:\n")
  print(summary(D$reads))
  cat("Reads/sample:\n")
  print(summary(D$reads/D$samples))
  cat("\n")
}

evaluate <- function(D) {
  ntotclusters <- nrow(D)
  D <- D[D$cluster %in% problem_clusters,]
  cat("Clusters:", nrow(D), "\n")
  cat("Prop of clusters:", nrow(D)/ntotclusters,"\n")  
  cat("Samples (mean):", mean(D$samples), "\n")
  cat("Sample prop. (mean of proportions):", mean(D$samples/(ncol(D)-3)), "\n")
  cat("Reads (mean of totals):", mean(D$reads), "\n")
  cat("Reads/sample (mean of reads/sample):", mean(D$reads/D$samples), "\n")
  cat("\n")
}

# Spikein clusters
spikein_clusters <- c("Miridae_cluster1",
                      "Miridae_cluster2",
                      "Miridae_cluster3",
                      "Anthocoridae_cluster1",
                      "Cecidomyiidae_cluster1",
                      "Braconidae_cluster1",
                      "Coccinellidae_cluster1")

# Get baseline values
cat("Baseline for litter samples\n")
cat("Total of:", ncol(mg_litter)-3,"samples\n")
evaluate(mg_litter)

cat("Baseline for lysate samples\n")
cat("Total of:", ncol(mg_lysate)-3,"samples\n")
evaluate(mg_lysate)

cat("Filter negative control and spikein occurrences for litter samples\n")
X1 <- mg_litter[!(mg_litter$cluster %in% mg_litter_control$cluster),]
X1 <- X1[!(X1$cluster %in% spikein_clusters),]
evaluate(X1)
W1 <- X1[X1$cluster %in% problem_clusters,]
W1$max_reads <- numeric(nrow(W1))
for (i in 1:nrow(W1)) W1$max_reads[i] <- as.numeric(max(W1[i,2:(ncol(W1)-2)]))
W1$rep_asv <- R2$ASV[match(W1$cluster,R2$cluster)]
W1$Order <- R2[match(W1$cluster,R2$cluster),"Order"]
W1$Species <- R2[match(W1$cluster,R2$cluster),"Species"]
write.table(W1, "mg_litter_problemclusters_wo_controls.tsv",row.names=FALSE,sep="\t")

cat("Filter negative control and spikein occurrences for lysate samples\n")
X2 <- mg_lysate[!(mg_lysate$cluster %in% mg_lysate_control$cluster),]
X2 <- X2[!(X2$cluster %in% spikein_clusters),]
evaluate(X2)
W2 <- X2[X2$cluster %in% problem_clusters,]
W2$max_reads <- numeric(nrow(W2))
for (i in 1:nrow(W2)) W2$max_reads[i] <- as.numeric(max(W2[i,2:(ncol(W1)-2)]))
W2$rep_asv <- R2$ASV[match(W2$cluster,R2$cluster)]
W2$Order <- R2[match(W2$cluster,R2$cluster),"Order"]
W2$Species <- R2[match(W2$cluster,R2$cluster),"Species"]
write.table(W2, "mg_lysate_problemclusters_wo_controls.tsv",row.names=FALSE,sep="\t")

cat("Filter negative and positive control occurrences for litter samples\n")
X1 <- X1[!(X1$cluster %in% mg_litter_pos_con$cluster),]
evaluate(X1)

cat("Filter negative and positive control occurrences for lysate samples\n")
X2 <- X2[!(X2$cluster %in% mg_lysate_pos_con$cluster),]
evaluate(X2)

cat("Filter negative control occurrences above 5% for mg samples\n")
# Compute proportion of total control samples for mg
D <- rbind(mg_litter_control[,c("cluster","reads","samples")],mg_lysate_control[,c("cluster","reads","samples")])
E <- aggregate(D[,2:3],by=list(cluster=D$cluster),FUN=sum)
ncontrols <- ncol(mg_litter_control)+ncol(mg_lysate_control)-6
E$prop_samples <- E$samples/ncontrols
discard_clusters <- E$cluster[E$prop_samples>0.05]
# Evaluate results
cat("Results for litter samples\n")
X1 <- mg_litter[!(mg_litter$cluster %in% discard_clusters),]
evaluate(X1)
cat("Results for lysate samples\n")
X2 <- mg_lysate[!(mg_lysate$cluster %in% discard_clusters),]
X2 <- X2[!(X2$cluster %in% spikein_clusters),]
evaluate(X2)

# Define function for evaluation ratios between max values
eval_ratios <- function(D, P) {
  ret <- D[D$cluster %in% P$cluster,]
  ret$pos_con_max_ratio <- 0.0
  matches <- match(ret$cluster,P$cluster)
  for (i in 1:nrow(ret)) {
    if (!is.na(matches[i])) {
      max_pos <- max(P[matches[i],2:(ncol(P)-3)])
      max_sample <- max(ret[i,2:(ncol(ret)-3)])
      ret$pos_con_max_ratio[i] <- max_pos/max_sample
    }
  }
  ret
}

cat("Filter negative control occurrences above 5% for litter samples\n")
mg_litter_control$prop_samples <- mg_litter_control$samples/(ncol(mg_litter_control)-3)
discard_clusters1 <- mg_litter_control$cluster[mg_litter_control$prop_samples>0.05]
X1 <- mg_litter[!(mg_litter$cluster %in% discard_clusters1),]
evaluate(X1)
Y1 <- eval_ratios(X1,mg_litter_pos_con)
cat("No. clusters shared with positive controls:", nrow(Y1),"\n")
cat("Summary of read ratios\n")
print(summary(Y1$pos_con_max_ratio))
cat("No. clusters with read ratios (pos_con to sample) > 2:", sum(Y1$pos_con_max_ratio>2.0),"\n")
cat("No. clusters with read ratios (pos_con to sample) > 5:", sum(Y1$pos_con_max_ratio>5.0),"\n")

cat("Filter negative control occurrences above 5% for lysate samples\n")
mg_lysate_control$prop_samples <- mg_lysate_control$samples/(ncol(mg_lysate_control)-3)
discard_clusters2 <- mg_lysate_control$cluster[mg_lysate_control$prop_samples>0.05]
X2 <- mg_lysate[!(mg_lysate$cluster %in% discard_clusters2),]
X2 <- X2[!(X2$cluster %in% spikein_clusters),]
evaluate(X2)
Y2 <- eval_ratios(X2, mg_lysate_pos_con)
cat("No. clusters shared with positive controls:", nrow(Y2),"\n")
cat("Summary of read ratios\n")
print(summary(Y2$pos_con_max_ratio))
cat("No. clusters with read ratios (pos_con to sample) > 2:", sum(Y2$pos_con_max_ratio>2.0),"\n")
cat("No. clusters with read ratios (pos_con to sample) > 5:", sum(Y2$pos_con_max_ratio>5.0),"\n")
