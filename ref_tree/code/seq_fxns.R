library(ape)
library(seqinr)

# Function for adding TipLabel column to taxonomy file
add_tiplabel <- function(meta) {

    meta$TipLabel <- paste0(meta$Genus,"_",meta$Species)
    meta
}

# Function for renaming sequences with tip labels according to
# metadata file. Matching is based on GenBank record ID
rename_seqs <- function(seqs, meta) {

    x <- names(seqs)
    w <- character()
    
    for (i in 1:length(x)) {
        w[i] <- strsplit(x[i]," ")[[1]][1]
    }

    w <- unlist(w)

    if (sum(is.na(match(w,meta$GenBank)))!=0) {
        cat("ERROR: There are sequences that do not match GenBank column entries in the metadata file\n")
        cat("The non-matching sequences are:\n")
        cat(w[is.na(match(w,meta$GenBank))],sep=",")
        cat("\n")
        return (NULL)
    }

    names(seqs) <- meta$TipLabel[match(w,meta$GenBank)]

    seqs
}

# Function for extracting coding sequences in DNAbin format (ape)
extract_coding <- function(seqs, meta) {

    if (class(seqs)!="DNAbin") {
        cat("ERROR: Expecting sequences in DNAbin format and not",class(seqs),"format\n")
        return (NULL)
    }
    x <- names(seqs)

    # Here are the standard stop codons in binary format. Note that APE uses 18, 28, 48 and 88 for T, C, G and A, respectively.
    stop_codon1 <- c("18","88","88")   # TAA in DNAbin
    stop_codon2 <- c("18","88","48")   # TAG in DNAbin

    for (i in 1:length(x)) {
        
        meta_idx <- match(x[i],meta$TipLabel,nomatch=NA)
        if (is.na(meta_idx)) {
            cat("ERROR: No match for ",x[i],"\n")
            return (NULL)
        }

        len <- meta$Stop[meta_idx] - meta$Start[meta_idx] + 1

        if (meta$Stop[meta_idx] > length(seqs[[i]])) {
            cat ("ERROR: Length of sequence", names(seqs)[i], "shorter than specified stop site. Please correct!\n")
            return (NULL)
        }
        if (len %% 3 != 0) {
            cat ("ERROR: Length for sequence", names(seqs)[i], "not divisible by 3. Please correct!\n")
            return (NULL)
        }
        last_codon <- as.character(seqs[[i]][(meta$Stop[meta_idx]-2):meta$Stop[meta_idx]])
        if (identical(last_codon,stop_codon1) || identical(last_codon,stop_codon2))
            seqs[[i]] <- seqs[[i]][meta$Start[meta_idx]:(meta$Stop[meta_idx]-3)]
        else
            seqs[[i]] <- seqs[[i]][meta$Start[meta_idx]:meta$Stop[meta_idx]]
    }

    seqs
}

# Function for extracting barcode sequences in DNAbin format (ape)
# based on alignment with a target barcode
extract_barcode <- function(seqs, meta) {

    if (class(seqs)!="DNAbin") {
        cat("ERROR: Expecting sequences in DNAbin format and not",class(seqs),"format\n")
        return (NULL)
    }
    x <- names(seqs)

    # Here are the standard stop codons in binary format.
    # Note that APE uses 18, 28, 48 and 88 for T, C, G and A, respectively.
    stop_codon1 <- c("18","88","88")   # TAA in DNAbin
    stop_codon2 <- c("18","88","48")   # TAG in DNAbin

    for (i in 1:length(x)) {
        
        meta_idx <- match(x[i],meta$TipLabel,nomatch=NA)
        if (is.na(meta_idx)) {
            cat("ERROR: No match for ",x[i],"\n")
            return (NULL)
        }

        len <- meta$aln_q_end[meta_idx] - meta$aln_q_start[meta_idx] + 1

        if (len > length(seqs[[i]])) {
            cat ("ERROR: Length of sequence", names(seqs)[i], "shorter than specified alignment. Please correct!\n")
            return (NULL)
        }

        # Find out codon start (assuming the number is the codon position of the first site
        codon_start <- ((meta$aln_t_start[meta_idx]+1) %% 3) + 1    ## target coordinate 652 == codon site 3
        if (codon_start == 2) len <- len - 2
        if (codon_start == 3) len <- len - 1
        len <- floor(len / 3)

        last_codon <- as.character(seqs[[i]][(meta$Stop[meta_idx]-2):meta$Stop[meta_idx]])
        if (identical(last_codon,stop_codon1) || identical(last_codon,stop_codon2))
            seqs[[i]] <- seqs[[i]][meta$Start[meta_idx]:(meta$Stop[meta_idx]-3)]
        else
            seqs[[i]] <- seqs[[i]][codon_start:(codon_start+len)]
    }

    seqs
}

# Function for converting aligned fasta sequences to nexus data block
fasta2nexus <- function(infile, outfile) {

    seqs <- read.alignment(infile, format="fasta")

    ntax <- length(seqs$seq)
    nchar <- nchar(seqs$seq[1])

    cat ("#NEXUS\n\nbegin data;\n\tdimensions ntax =", ntax, " nchar =", nchar, ";\n", sep="", file=outfile)
    cat ("\tformat datatype=dna gap=- interleave=no;\n\tmatrix\n", sep="", file=outfile, append=TRUE)

    for (i in 1:length(seqs$seq)) {
        cat(seqs$nam[i], "\t", toupper(seqs$seq[i]), "\n", sep="", file=outfile, append=TRUE)
    }

    cat("\t;\nend;\n", file=outfile, append=TRUE)
}

# Function for finding stop codons in aligned DNA sequences
# in ape format, using trans fxn in ape
no_stop_codons <- function(seqs, codonstart) {

    aa_seqs <- trans(seqs,code=5, codonstart=codonstart)
    x <- logical()
    for (i in 1:length(aa_seqs))
        x[i] <- is.na(match("*",as.character(aa_seqs[[i]])))

    x
}

