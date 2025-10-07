library(data.table)

filename <- "~/dev/figshare-repos/iba/raw_data/v6/CO1_asv_counts_MG.tsv"

D <- fread("~/dev/figshare-repos/iba/raw_data/v6/CO1_asv_counts_MG.tsv", nrows=0)
M <- read.delim("~/dev/figshare-repos/iba/raw_data/v6/CO1_sequencing_metadata_MG.tsv")
T <- read.delim("~/dev/figshare-repos/iba/processed_data/v3/asv_taxonomy_MG.tsv")

b <- unique(M$sequencing_batch)

for (i in 1:length(b)) {
    cat("Now processing ", b[i],"\n",sep="") 
    samples <- M$sampleID_NGI[M$sequencing_batch==b[i]]
    cols <- c(1,which(colnames(D) %in% samples))
    X <- fread(filename,select=cols)
    write.table(data.frame(X[which(rowSums(X[,-1])!=0),]),file=paste0(b[i],"_mg.tsv"),sep="\t",row.names=FALSE)
}

