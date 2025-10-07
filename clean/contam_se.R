M <- read.delim("~/dev/figshare-repos/iba/raw_data/v6/CO1_sequencing_metadata_SE.tsv")
T <- read.delim("~/dev/figshare-repos/iba/processed_data/v3/asv_taxonomy_SE.tsv")

b <- unique(M$sequencing_batch)

for (i in 1:length(b)) {

    extract <- function(D, batch_name, sample_type) {
        samples <- M$sampleID_NGI[M$sequencing_batch==batch_name & M$lab_sample_type==sample_type]
        idx <- c(1,which(colnames(D) %in% samples))
        if (length(idx) > 1) {
            X <- D[,idx]
            if (length(idx)==2)
                X <- X[X[,2]!=0,]
            else
                X <- X[rowSums(X[,-1])!=0,]
            return (X)
        } else { return (data.frame(ASV_ID=character())) }
    }
  
    filename <- paste0(b[i],"_se.tsv")
    D <- read.delim(filename)
    lab_samples <- extract(D, b[i], "sample")
    pos_control <- extract(D, b[i], "positive_control")
    buffer_blank <- extract(D, b[i], "buffer_blank")
    extraction_neg <- extract(D, b[i], "extraction_neg")
    pcr_neg <- extract(D, b[i], "pcr_neg")
    mock_sample <- extract(D, b[i], "mock_sample")

    cat("\nResults for",b[i],"\n")
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
        if (ncol(X)==2)
            tot_reads <- X[,2]
        else
            tot_reads <- rowSums(X[,-1])
        X <- X[order(tot_reads,decreasing=TRUE),]
        print(T[match(X$ASV_ID[1:10],T$ASV),])
        cat("Top ten annotations in sample order:\n")
        if (ncol(X)==2)
            tot_cols <- X[,2]!=0
        else
            tot_cols <- rowSums(X[,-1]!=0)
        X <- X[order(tot_cols,decreasing=TRUE),]
        print(T[match(X$ASV_ID[1:10],T$ASV),])

        cat("Top ten annotations for overlap with lab samples in read order:\n")
        X <- X[X$ASV_ID %in% lab_samples$ASV_ID,]
        if (ncol(X)==2)
            tot_reads <- X[,2]
        else
            tot_reads <- rowSums(X[,-1])
        X <- X[order(tot_reads,decreasing=TRUE),]
        print(T[match(X$ASV_ID[1:10],T$ASV),])
        cat("Top ten annotations for overlap with lab samples in sample order:\n")
        if (ncol(X)==2)
            tot_cols <- X[,2]!=0
        else
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
    cat("Distribution of total read numbers of shared ASVs across samples\n")
    print(summary(a))
}

