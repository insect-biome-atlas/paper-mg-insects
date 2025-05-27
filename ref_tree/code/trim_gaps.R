library(ape)

# Retrieve the command line arguments
args <- commandArgs(TRUE)
infile <- args[1]
outfile <- args[2]
max_gap_prop <- args[3]

seqs <- read.FASTA(infile)

seqs <- as.matrix(seqs)

seqs <- del.colgapsonly(seqs, threshold=max_gap_prop)

write.FASTA(seqs,file=outfile)

