library(rentrez)

fetch_fasta <- function(D, seq_file) {

    for (i in 1:nrow(D)) {

        acsn <- D$GenBank[i]

        res <- entrez_search(db="nucleotide", term=paste0(acsn,"[ACCN]"))  

        seq <- entrez_fetch(db="nucleotide", id=res$ids[1], rettype="fasta")
        cat(file=seq_file, append=TRUE, sep="", seq)
    }
}

