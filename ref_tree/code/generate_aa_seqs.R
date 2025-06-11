# Generate unaligned aa sequence files

# Retrieve the command line arguments
args <- commandArgs(TRUE)
infile <- args[1]
outfile <- args[2]

# Libraries needed
library(ape)

# Read in fasta sequences
seqs <- read.FASTA(infile)

# Note: the trans function does not quite handle ambiguities correctly, generating a few
# warnings in the pal2nal step. This should have negligible effect on the analysis.
aa_seqs <- trans(seqs, code=5)
warnings()

write.FASTA(aa_seqs, file=outfile)

