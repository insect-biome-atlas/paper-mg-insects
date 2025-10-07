library(data.table)

# Questions:
# Is it one or a few samples?
# Is it positive control or mock samples?
# Is it general contamination from SE samples?

D <- fread("~/dev/figshare-repos/iba/raw_data/v6/CO1_asv_counts_MG.tsv")
M <- read.delim("~/dev/figshare-repos/iba/raw_data/v6/CO1_sequencing_metadata_MG.tsv")
T <- read.delim("~/dev/figshare-repos/iba/processed_data/v3/asv_taxonomy_MG.tsv")

b <- unique(M$sequencing_batch)

batch <- list()
for (i in 1:length(b)) {
    samples <- M$sampleID_NGI[M$sequencing_batch==b[i]]
    cols <- c(1,which(colnames(D) %in% samples))
    X <- D[,..cols]
    batch[[i]] <- X[which(rowSums(X[,-1])!=0),]
}


for (i in 1:length(batch)) {

    extract <- function(D, batch_name, sample_type) {
        samples <- M$sampleID_NGI[M$sequencing_batch==batch_name & M$lab_sample_type==sample_type]
        if (length(samples) > 0) {
            idx <- c(1,which(colnames(D) %in% samples))
            X <- D[,..idx]
            X <- X[rowSums(X[,-1])!=0,]
            return (X)
        } else { return (data.frame(ASV_ID=character())) }
    }
    
    lab_samples <- extract(batch[[i]], b[i], "sample")
    pos_control <- extract(batch[[i]], b[i], "positive_control")
    buffer_blank <- extract(batch[[i]], b[i], "buffer_blank")
    extraction_neg <- extract(batch[[i]], b[i], "extraction_neg")
    pcr_neg <- extract(batch[[i]], b[i], "pcr_neg")
    mock_sample <- extract(batch[[i]], b[i], "mock_sample")

    cat("\nBatch",b[i],"\n")
    cat("Sample asvs:",nrow(lab_samples),"\n")

    results <- function(X,sample_type) {
        cat(sample_type, "asvs:", nrow(X),"\n")
        cat("Shared asvs:",sum(lab_samples$ASV_ID %in% X$ASV_ID),"\n")
        cat("Prop of sample asvs:",sum(lab_samples$ASV_ID %in% X$ASV_ID) / nrow(lab_samples),"\n")
        cat("Prop of", sample_type, "asvs:",sum(lab_samples$ASV_ID %in% X$ASV_ID) / nrow(X),"\n")
    }

    annotation <- function(X) {
        if (nrow(X) == 0) {
            cat ("Empty set: no annotations\n")
            return ()
        }
        cat("Top ten annotations in read order:\n")
        tot_reads <- rowSums(X[,-1])
        X <- X[order(tot_reads,decreasing=TRUE),]
        print(T[match(X$ASV_ID[1:10],T$ASV),])
        cat("Top ten annotations in sample order:\n")
        tot_cols <- rowSums(X[,-1]!=0)
        X <- X[order(tot_cols,decreasing=TRUE),]
        print(T[match(X$ASV_ID[1:10],T$ASV),])
        cat("Top ten annotations for overlap with lab samples in read order:\n")
        X <- X[X$ASV_ID %in% lab_samples$ASV_ID,]
        tot_reads <- rowSums(X[,-1])
        X <- X[order(tot_reads,decreasing=TRUE),]
        print(T[match(X$ASV_ID[1:10],T$ASV),])
        cat("Top ten annotations for overlap with lab samples in sample order:\n")
        tot_cols <- rowSums(X[,-1]!=0)
        X <- X[order(tot_cols,decreasing=TRUE),]
        print(T[match(X$ASV_ID[1:10],T$ASV),])
    }

    results(pos_control, "Positive control")
    annotation(pos_control)
    results(buffer_blank, "Buffer blank")
    annotation(buffer_blank)
    results(extraction_neg, "Extraction neg")
    annotation(extraction_neg)
    results(pcr_neg, "PCR neg")
    annotation(pcr_neg)
    results(mock_sample, "Mock sample")
    annotation(mock_sample)

    control_asvs <- unique(c(buffer_blank$ASV_ID,extraction_neg$ASV_ID, pcr_neg$ASV_ID))

    X <- lab_samples[lab_samples$ASV_ID %in% pos_control$ASV_ID,]
    cat("Num samples:", ncol(X)-1, "\n")
    cat("Shared asvs across lab samples and pos controls: ", nrow(X), "\n")
    X <- X[!X$ASV_ID %in% control_asvs,]
    cat("Shared asvs across lab samples and pos controls w/o blanks/negs: ", nrow(X), "\n")
    a <- colSums(X[,-1])
    cat("Distribution of shared ASVs across samples\n")
    print(summary(a))
}

#rm(D)

#D1 <- fread("~/dev/figshare-repos/iba/raw_data/v6/CO1_asv_counts_SE.tsv")
#M1 <- read.delim("~/dev/figshare-repos/iba/raw_data/v6/CO1_sequencing_metadata_SE.tsv")

#b1 <- unique(M1$sequencing_batch)

#batch1 <- list()
#for (i in 1:length(b1)) {
#    samples <- M1$sampleID_NGI[M1$sequencing_batch==b1[i]]
#    cols <- which(colnames(D1) %in% samples)
#    X <- D1[,..cols]
#    row.names(X) <- D1$ASV_ID
#    batch1[[i]] <- data.frame(X[which(rowSums(X)!=0),])
#}

#rm(D1)

