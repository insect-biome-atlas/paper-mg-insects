# Extend mycetophilid data table with inter-sequence differences
library(ape)

D <- read.delim("mycetophilid_data.tsv")

# Compute distances
seqs <- read.FASTA("mycetophilid_seqs_aligned.fasta")
distances <- dist.dna(seqs,model="N",as.matrix=TRUE)

for (i in 1:nrow(D)) {
  D$dist[i] <- distances[D$asv_se[i],D$asv_mg[i]]
}

write.table(D,"mycetophilid_data_and_dists.tsv",sep="\t",row.names=FALSE)
