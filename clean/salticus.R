M <- read.delim("~/dev/figshare-repos/iba/raw_data/v6/CO1_sequencing_metadata_SE.tsv")
T <- read.delim("~/dev/figshare-repos/iba/processed_data/v3/asv_taxonomy_SE.tsv")

extract <- function(D, sample_type, M) {
    idx <- c(1,which(colnames(D) %in% M$sampleID_NGI[M$lab_sample_type %in% sample_type]))
    if (length(idx) < 2)
        return (data.frame(ASV_ID=character()))
    D[,idx]
}

fie_sam <- list()
control <- list()

batches <- unique(M$sequencing_batch)
batches <- batches[order(batches)]
for (i in 1:length(batches)) {

    batch <- batches[i]
    cat("Processing",batch,"\n")
    D <- read.delim(paste0(batch,"_se.tsv"))
    fie_sam[[i]] <- extract(D, "sample", M)
    control[[i]] <- extract(D, c("buffer_blank","extraction_neg","pcr_neg"), M)
}

salticus_samples <- data.frame(Sample=character(), ASV_ID=character(), Genus=character(), Species=character(), BOLD_bin=character(), reads=numeric())
salticus_controls <- data.frame(Sample=character(), Sample_type=character(), ASV_ID=character(), Genus=character(), Species=character(), BOLD_bin=character(), reads=numeric())

for (i in 1:length(batches)) {

    batch <- fie_sam[[i]]
    cat("Processing samples for",batches[i],"\n")
    X <- batch[T$Species[match(batch$ASV_ID,T$ASV)]=="Salticus scenicus",]

    for(j in 2:ncol(X)) {

        idx <- min(which(X[,j]==max(X[,j])))
        reads <- X[idx,j]
        if (reads > 0) {
            Sample <- colnames(X)[j]
            ASV_ID <- X$ASV_ID[idx]
            Genus <- T$Genus[match(ASV_ID, T$ASV)]
            Species <- T$Species[match(ASV_ID, T$ASV)]
            BOLD_bin <- T$BOLD_bin[match(ASV_ID, T$ASV)]
            if (Species!="Salticus scenicus")
                cat("ERROR: ", Species, "for ASV", ASV_ID, "and Sample", Sample, "\n")
            salticus_samples <- rbind(salticus_samples, list(Sample=Sample, ASV_ID=ASV_ID, Genus=Genus, Species=Species, BOLD_bin=BOLD_bin, reads=reads))
        }
    }
}

for (i in 1:length(batches)) {

    batch <- control[[i]]
    cat("Processing controls for",batches[i],"\n")
    X <- batch[T$Species[match(batch$ASV_ID,T$ASV)]=="Salticus scenicus",]

    for(j in 2:ncol(X)) {

        idx <- min(which(X[,j]==max(X[,j])))
        reads <- X[idx,j]
        if (reads > 0) {
            Sample <- colnames(X)[j]
            Sample_type <- M$lab_sample_type[match(Sample,M$sampleID_NGI)]
            ASV_ID <- X$ASV_ID[idx]
            Genus <- T$Genus[match(ASV_ID, T$ASV)]
            Species <- T$Species[match(ASV_ID, T$ASV)]
            BOLD_bin <- T$BOLD_bin[match(ASV_ID, T$ASV)]
            salticus_controls <- rbind(salticus_controls, list(Sample=Sample, Sample_type=Sample_type, ASV_ID=ASV_ID, Genus=Genus, Species=Species, BOLD_bin=BOLD_bin, reads=reads))
        }
    }
}

