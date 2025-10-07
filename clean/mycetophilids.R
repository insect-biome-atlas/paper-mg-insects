# Script for examining strange overlaps in mycetophilid species
# reported by Jostein Kjaerandsen

library(ape)

# Species concerned
species <- c(
"Allocotocera pulchella",
"Allodia triangularis",
"Boletina edwardsi",
"Boletina gripha",
"Boletina griphoides",
"Boletina nigricans",
"Boletina onegensis",
"Brevicornu sericoma",
"Brevicornu verralli",
"Coelosia tenella",
"Cordyla flaviceps",
"Cordyla parvipalpis",
"Ectrepesthoneura hirta",
"Exechia confinis",
"Exechia dorsalis",
"Exechia festiva",
"Exechia fusca",
"Exechia macula",
"Exechia pseudocincta",
"Leia cylindrica",
"Leia winthemi",
"Mycetophila ichneumonea",
"Mycetophila signatoides",
"Mycomya affinis",
"Mycomya annulata",
"Mycomya fimbriata",
"Mycomya shermani",
"Mycomya trilineata",
"Opsion nigra",
"Phronia egregia",
"Rymosia domestica",
"Rymosia fasciata",
"Rymosia placida",
"Rymosia signatipes",
"Sciophila hirta",
"Sciophila salassea",
"Synapha vitripennis",
"Syntemna hungarica",
"Tarnania fenestralis",
"Tarnania tarnanii",
"Trichonta atricauda",
"Zygomyia notata")

# Counts of these species
samples <- c(
1,
1,
2,
1,
1,
3,
2,
2,
10,
3,
2,
1,
2,
1,
1,
2,
6,
1,
1,
1,
2,
3,
1,
4,
3,
3,
6,
3,
3,
1,
3,
1,
3,
3,
2,
1,
1,
1,
1,
2,
1,
2)

write.table(data.frame(list(species=species, samples=samples)),"mycetophilid_mystery_species.tsv",sep="\t",row.names=FALSE)
# Check how many of these taxa are in the controls

# Read in metadata
cat("Reading in metadata\n")
M1 <- read.delim("~/dev/figshare-repos/iba/raw_data/v6/CO1_sequencing_metadata_SE.tsv")
M2 <- read.delim("~/dev/figshare-repos/iba/raw_data/v6/CO1_sequencing_metadata_MG.tsv")

# Read in cluster counts
cat("Reading in cluster counts\n")
C1 <- read.delim("~/dev/figshare-repos/iba/processed_data/v3/cluster_counts_SE.tsv")
C2 <- read.delim("~/dev/figshare-repos/iba/processed_data/v3/cluster_counts_MG.tsv")
cat("Completed reading in cluster counts\n")

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
    add_reads_samples(D)
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

# Get se lysate data
foo <- function(st) { get_sample_data(C1, get_samples(M1, st, "CO1_lysate_2019_SE")) }
se_lysate <- foo("sample")
se_lysate_pos_con <- foo("positive_control")
se_lysate_control <- foo(c("buffer_blank", "extraction_neg", "pcr_neg"))

# Get se homogenate data
foo <- function(st) { get_sample_data(C1, get_samples(M1, st, "CO1_homogenate_2019_SE")) }
se_homogenate <- foo("sample")
se_homogenate_pos_con <- foo("positive_control")
se_homogenate_control <- foo(c("buffer_blank", "buffer_blank_art_spikes", "extraction_neg", "pcr_neg"))

# Get se soil data
foo <- function(st) { get_sample_data(C1, get_samples(M1, st, "CO1_soil_litter_2019_SE")) }
se_litter <- foo("sample")
se_litter_pos_con <- foo("positive_control")
se_litter_control <- foo(c("buffer_blank", "buffer_blank_art_spikes", "extraction_neg", "pcr_neg"))

# Data frame for results
X <- data.frame(list(
            species=character(),
            cluster_se=character(),
            cluster_mg=character(),
            asv_se=character(),
            asv_mg=character(),
            se_lysate=numeric(),
            se_homogenate=numeric(),
            se_litter=numeric(),
            mg_lysate=numeric(),
            mg_litter=numeric(),
            se_lysate_pos_con=numeric(),
            se_lysate_control=numeric(),
            se_homogenate_pos_con=numeric(),
            se_homogenate_control=numeric(),
            se_litter_pos_con=numeric(),
            se_litter_control=numeric(),
            mg_lysate_pos_con=numeric(),
            mg_lysate_control=numeric(),
            mg_litter_pos_con=numeric(),
            mg_litter_control=numeric(),
            mg_litter_samples=numeric(),
            mg_litter_pos_con_samples=numeric()
            ))

# Read in taxonomy info
T1 <- read.delim("~/dev/figshare-repos/iba/processed_data/v3/cluster_taxonomy_SE.tsv")
T2 <- read.delim("~/dev/figshare-repos/iba/processed_data/v3/cluster_taxonomy_MG.tsv")
T1 <- T1[T1$representative==1,]
T2 <- T2[T2$representative==1,]

get_reads <- function(D, cluster) {
    idx <- match(cluster,D$cluster)
    if (is.na(idx))
        0
    else
        D$reads[idx]
}

get_nsamples <- function(D, cluster) {
    idx <- match(cluster,D$cluster)
    if (is.na(idx))
        0
    else
        D$samples[idx]
}

for (i in 1:length(species)) {

    sp <- species[i]
    clusters_se <- T1$cluster[T1$Species==sp]
    clusters_mg <- T2$cluster[T2$Species==sp]

    if (length(clusters_se) == 0 || length(clusters_mg) == 0) {
        if (length(clusters_se)==0) cat("No clusters in SE for", sp, "\n")
        if (length(clusters_mg)==0) cat("No clusters in MG for", sp, "\n")
        next
    }

    for (j in 1:length(clusters_se)) {
        cluster_se <- clusters_se[j]
        for (k in 1:length(clusters_mg)) {
            cluster_mg <- clusters_mg[k]
            cat("Processing",cluster_se,"and",cluster_mg,"\n")
            X <- rbind(X, list(
                    species = sp,
                    cluster_se = cluster_se,
                    cluster_mg = cluster_mg,
                    asv_se = T1$ASV[T1$cluster==cluster_se],
                    asv_mg = T2$ASV[T2$cluster==cluster_mg],
                    se_lysate     = get_reads(se_lysate, cluster_se),
                    se_homogenate = get_reads(se_homogenate, cluster_se),
                    se_litter     = get_reads(se_litter, cluster_se),
                    mg_lysate     = get_reads(mg_lysate, cluster_mg),
                    mg_litter     = get_reads(mg_litter, cluster_mg),
                    se_lysate_pos_con     = get_reads(se_lysate_pos_con, cluster_se),
                    se_lysate_control     = get_reads(se_lysate_control, cluster_se),
                    se_homogenate_pos_con = get_reads(se_homogenate_pos_con, cluster_se),
                    se_homogenate_control = get_reads(se_homogenate_control, cluster_se),
                    se_litter_pos_con     = get_reads(se_homogenate_pos_con, cluster_se),
                    se_litter_control     = get_reads(se_homogenate_control, cluster_se),
                    mg_lysate_pos_con     = get_reads(mg_lysate_pos_con, cluster_mg),
                    mg_lysate_control     = get_reads(mg_lysate_control, cluster_mg),
                    mg_litter_pos_con     = get_reads(mg_litter_pos_con, cluster_mg),
                    mg_litter_control     = get_reads(mg_litter_control, cluster_mg),
                    mg_litter_samples     = get_nsamples(mg_litter, cluster_mg),
                    mg_litter_pos_con_samples = get_nsamples(mg_litter_pos_con, cluster_mg)
                ))
        }
    }
}

# Write result table
write.table(X, "mycetophilid_data.tsv", row.names=FALSE, sep="\t")

# Get sequences occurring in the table 
seq1 <- read.FASTA("~/dev/figshare-repos/iba/processed_data/v3/cluster_reps_SE.fasta")
seq2 <- read.FASTA("~/dev/figshare-repos/iba/processed_data/v3/cluster_reps_MG.fasta")

rename_seqs <- function(seq) {
    old_names <- names(seq)
    words <- strsplit(old_names, split=" ")
    new_names <- character(length(old_names))
    for (i in 1:length(seq))
        new_names[i] <- words[[i]][1]
    names(seq) <- new_names
    seq
}

seq1 <- rename_seqs(seq1)
seq2 <- rename_seqs(seq2)

# Merge sequences, make sure they are unique and make them in-frame
seqs <- c(seq1[unique(X$asv_se)], seq2[unique(X$asv_mg)])
seqs <- seqs[unique(names(seqs))]
for (i in 1:length(seqs))
  seqs[[i]] <- seqs[[i]][-1]
write.FASTA(seqs, "mycetophilid_seqs.fasta")

seqs_aa <- trans(seqs, code=5, codonstart=1)
write.FASTA(seqs_aa, "mycetophilid_seqs_aa.fasta")
