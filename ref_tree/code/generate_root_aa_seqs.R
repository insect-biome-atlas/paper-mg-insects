# Generate unaligned aa sequence files
# This version of the script deals with input sequences
# with different translation tables

# Retrieve the command line arguments
args <- commandArgs(TRUE)
infile <- args[1]
metafile <- args[2]
outfile <- args[3]

# Libraries needed
library(ape)

# Read in fasta sequences
seqs <- read.FASTA(infile)

# Read in metadata (translation table)
meta <- read.delim(metafile)

# Note: the trans function does not quite handle ambiguities correctly, generating a few
# warnings in the pal2nal step. This should have negligible effect on the analysis.

# Translate by genetic code (translation table)
seqs_code <- function(seqs, code) { seqs[names(seqs) %in% meta$TipLabel[meta$TranslationTable==code]] }
aa_seqs_code <- function(seqs, code) { trans(seqs_code(seqs, code), code=code) }

aa_seqs5 <- aa_seqs_code(seqs, 5)
warnings()

aa_seqs1 <- aa_seqs_code(seqs, 1)
warnings()

aa_seqs2 <- aa_seqs_code(seqs, 2)
warnings()

aa_seqs11 <- trans(seqs_code(seqs, 11), code=1)  # APE only implements codes 1-6. Code 11 is identical to 1 except for start codons.
warnings()

write.FASTA(aa_seqs5, file=outfile)
write.FASTA(aa_seqs1, file=outfile, append=TRUE)
write.FASTA(aa_seqs2, file=outfile, append=TRUE)
write.FASTA(aa_seqs11, file=outfile, append=TRUE)

