# Function for filtering se traps based on number of samples
filter_se_traps <- function(min_samples) {
    D <- read.delim("../data/otu_sample_meta_se.tsv")
    X <- unique(D[,c("trapID","sampleID_NGI")])
    trap_samples <- table(X$trapID)
    names(trap_samples[trap_samples >= min_samples])
}

