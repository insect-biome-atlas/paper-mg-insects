library(rentrez)
library(ape)

# TODO: Set to the name of the downloaded hit table from NCBI, generated from
# the sequence file (mycetophilid_seqs.fasta)
D <- read.delim("AX9U00H1015-Alignment-HitTable.csv", sep=",", header=FALSE)
D <- D[,1:3]
colnames(D) <- c("ASV_ID","genbank_id","identity")
D <- D[!duplicated(D$ASV_ID),]

D$taxonomy <- character(nrow(D))
D$organism <- character(nrow(D))
for (i in 1:nrow(D)) {
  
  cat("Fetching NCBI record for ASV_ID", i, ":", D$ASV_ID[i],"\n")

  X <- entrez_fetch(db="nuccore", id=D$genbank_id[i], rettype="gb", retmode="xml")
  temp <- strsplit(X,split="GBSeq_taxonomy")
  D$taxonomy[i] <- substr(temp[[1]][2],start=2,stop=nchar(temp[[1]][2])-2)
  temp <- strsplit(X,split="GBSeq_organism")
  D$organism[i] <- substr(temp[[1]][2],start=2,stop=nchar(temp[[1]][2])-2)
}

# Write table
write.table(D,"mycetophilid_seqs_ncbi_hits.tsv", sep="\t", row.names=FALSE)
