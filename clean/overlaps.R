library(data.table)

# Read in sample metadata
M1 <- read.delim("~/dev/figshare-repos/iba/raw_data/v6/CO1_sequencing_metadata_SE.tsv")
M2 <- read.delim("~/dev/figshare-repos/iba/raw_data/v6/CO1_sequencing_metadata_MG.tsv")

# Define useful functions
extract <- function(D,sample_type, M) {
    idx <- c(1,which(colnames(D) %in% M$sampleID_NGI[M$lab_sample_type %in% sample_type]))
    if (length(idx) < 2)
        return (data.frame(ASV_ID=character()))
    X <- D[,idx]
    if (ncol(X) > 2) {
        X <- X[rowSums(X[,-1])!=0,]
        reads <- rowSums(X[,-1])
        samples <- rowSums(X[,-1]!=0)
        X$reads <- reads
        X$samples <- samples
    }
    else {
        X <- X[X[,2]!=0,]
        X$reads <- X[,2]
        X$samples <- 1
    }
    return(X)
}

eval_res <- function(X, P) {

    cat("Asvs shared with pos_con: ", sum(X$ASV_ID %in% P$ASV_ID), "\n", sep="")
    cat("Prop asvs shared with pos_con: ", sum(X$ASV_ID %in% P$ASV_ID) / nrow(X), "\n", sep="")
    cat("Mean reads of asvs shared with pos_con: ", mean(X$reads[X$ASV_ID %in% P$ASV_ID]), "\n", sep="")
    cat("Median reads of asvs shared with pos_con: ", median(X$reads[X$ASV_ID %in% P$ASV_ID]), "\n", sep="")

    cat("Asvs shared with se: ", sum(X$ASV_ID %in% se_asv), "\n", sep="")
    cat("Prop asvs shared with se: ", sum(X$ASV_ID %in% se_asv) / nrow(X), "\n", sep="")
    cat("Mean reads of asvs shared with se: ", mean(X$reads[X$ASV_ID %in% se_asv]), "\n", sep="")
    cat("Median reads of asvs shared with se: ", median(X$reads[X$ASV_ID %in% se_asv]), "\n", sep="")
}

remove_controls <- function(X, C, cutoff) {
    rmv_asvs <- C$ASV_ID[C$samples/(ncol(C)-3) > cutoff]
    X[!(X$ASV_ID %in% rmv_asvs),]
}

# Collect data
pos_con1 <- list()
buf_bla1 <- list()
ext_neg1 <- list()
pcr_neg1 <- list()
fie_sam1 <- list()
control1 <- list()

batches1 <- unique(M1$sequencing_batch)
batches1 <- batches1[order(batches1)]
for (i in 1:length(batches1)) {

    batch <- batches1[i]
    D <- read.delim(paste0(batch,"_se.tsv"))
    pos_con1[[i]] <- extract(D, "positive_control", M1)
    buf_bla1[[i]] <- extract(D, "buffer_blank", M1)
    ext_neg1[[i]] <- extract(D, "extraction_neg", M1)
    pcr_neg1[[i]] <- extract(D, "pcr_neg", M1)
    fie_sam1[[i]] <- extract(D, "sample", M1)
    control1[[i]] <- extract(D, c("buffer_blank","extraction_neg","pcr_neg"), M1)
}

pos_con2 <- list()
buf_bla2 <- list()
ext_neg2 <- list()
pcr_neg2 <- list()
fie_sam2 <- list()
control2 <- list()

batches2 <- unique(M2$sequencing_batch)
batches2 <- batches2[order(batches2)]
for (i in 1:length(batches2)) {

    batch <- batches2[i]
    D <- read.delim(paste0(batch,"_mg.tsv"))
    pos_con2[[i]] <- extract(D, "positive_control", M2)
    buf_bla2[[i]] <- extract(D, "buffer_blank", M2)
    ext_neg2[[i]] <- extract(D, "extraction_neg", M2)
    pcr_neg2[[i]] <- extract(D, "pcr_neg", M2)
    fie_sam2[[i]] <- extract(D, "sample", M2)
    control2[[i]] <- extract(D, c("buffer_blank","extraction_neg","pcr_neg"), M2)
}

# Assemble all field sample ASVs for SE and MG
se_asv <- character()
for (i in 1:length(batches1)) {
    se_asv <- unique(c(se_asv, fie_sam1[[i]]$ASV_ID))
}

mg_asv <- character()
for (i in 1:length(batches2)) {
    mg_asv <- unique(c(mg_asv, fie_sam2[[i]]$ASV_ID))
}

# Analyze batchwise overlaps for MG samples
for (i in 1:length(batches2)) {

    cat("\nResults for ", batches2[[i]], "\n", sep="")
    cat("No. samples: ", ncol(fie_sam2[[i]]) - 2, "\n", sep="")
    cat("No. asvs: ", nrow(fie_sam2[[i]]), "\n", sep="")
    cat("No. pos control samples: ", ncol(pos_con2[[i]]) - 3, "\n", sep="")
    cat("No. pos control asvs: ", nrow(pos_con2[[i]]), "\n", sep="")

    cat("\nAfter removing all control asvs\n")
    X <- remove_controls(fie_sam2[[i]], control2[[i]], 0.0)
    eval_res(X, pos_con2[[i]])

    cat("\nWithout removing any asvs\n")
    eval_res(fie_sam2[[i]],pos_con2[[i]])
}


# Analyze batch16 with the SE samples included
D1 <- fread("~/dev/figshare-repos/iba/raw_data/v6/CO1_asv_counts_SE.tsv", nrows=0)
D2 <- fread("~/dev/figshare-repos/iba/raw_data/v6/CO1_asv_counts_MG.tsv", nrows=0)


nas2zeros <- function(X) {
    for (i in 1:nrow(X))
        for (j in 1:ncol(X))
            if (is.na(X[i,j]))
                X[i,j] <- 0
    X
}

# TODO: This results in duplicate rows for repeated ASV_IDs
control16 <- merge(control2[[4]], control1[[9]], all.x=TRUE, all.y=TRUE, no.dups=FALSE)
control16 <- nas2zeros(control16)
pos_con16 <- merge(pos_con2[[4]], pos_con1[[9]], all.x=TRUE, all.y=TRUE, no.dups=FALSE)
pos_con16 <- nas2zeros(pos_con16)
control16 <- aggregate(control16[,-1],by=list(ASV_ID=control16$ASV_ID),FUN=sum)
pos_con16 <- aggregate(pos_con16[,-1],by=list(ASV_ID=pos_con16$ASV_ID),FUN=sum)

cat("\nResults for batch16 with merged control samples\n")
cat("No. samples: ", ncol(fie_sam2[[4]]) - 2, "\n", sep="")
cat("No. asvs: ", nrow(fie_sam2[[4]]), "\n", sep="")
cat("No. pos control samples: ", ncol(pos_con16) - 3, "\n", sep="")
cat("No. pos control asvs: ", nrow(pos_con16), "\n", sep="")

cat("\nAfter removing all control asvs\n")
X <- remove_controls(fie_sam2[[4]], control16, 0.0)
eval_res(X, control16)

cat("\nWithout removing any asvs\n")
eval_res(fie_sam2[[4]],control16)

