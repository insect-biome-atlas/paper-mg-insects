# Generate unaligned aa sequence files

# Retrieve the command line arguments
args <- commandArgs(TRUE)
infile <- args[1]
outfile <- args[2]

# Libraries needed
library(seqinr)

# Read in fasta sequences
seqs <- read.fasta(infile)

# Note: the trans function does not quite handle ambiguities correctly, generating a few
# warnings in the pal2nal step. This should have negligible effect on the analysis.
aa_seqs <- lapply(seqs,FUN=function(x){translate(x,numcode=5,ambiguous=TRUE)})
warnings()

write.fasta(aa_seqs, names=names(aa_seqs), file.out=outfile)

