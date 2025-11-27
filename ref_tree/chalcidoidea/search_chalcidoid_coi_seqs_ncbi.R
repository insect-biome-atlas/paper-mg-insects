# Script for retrieving chalcidoid coi sequences from ncbi

library(ape)
library(rentrez)
library(stringr)
library(Biostrings)
source("gb_fxns.R")

# Define function for getting metadata from gb record
get_coi_meta <- function(gb_record) {

    a <- strsplit(gb_record,split="\n")[[1]]

    cds_text <- get_cds_text(a, c("coi","cox1"))
    if (is.na(cds))
        return (NA)
    res <- get_start_stop(cds_text[1])

    # Extract codon start
    codon_start <- extract_elem("codon_start",cds_text)

    # Extract transl_table
    transl_table <- extract_elem("transl_table",cds_text)

    c(res, list(codon_start=codon_start, transl_table=transl_table))
}

# Define function for fetching sequences
fetch_sequence <- function(D, TipLabel, taxa) {

    for (taxon in taxa) {
        cat("Fetching sequence for taxon/a: ",taxa,"\n")
        D <- fetch_mt_sequence(D, TipLabel, taxa)
        if (D$gb_accn[match(TipLabel,D$TipLabel)]=="")
            D <- fetch_coi_sequence(D, TipLabel, taxa)
    }
    return (D)
}

# Define function for fetching mitochondrial genome sequence if any
fetch_mt_sequence <- function(D, TipLabel, taxa) {

    success <- FALSE
    for (taxon in taxa) {
        
        if (success==TRUE)
            break

        cat("Fetching mitochondrion for taxon: ",taxon,"\n")

        res <- entrez_search(db="nucleotide", term=paste0("(mitochondrion[TITL] AND ",taxon,"[ORGN]) NOT UNVERIFIED[TITL]"))

        if (res$count > 0) {
       
            for (i in 1:length(res$ids)) {
                # Get cds info
                gb_record <- entrez_fetch(db="nucleotide", id=res$ids[i], rettype="gb", retmode="text")
                gb_txt <- strsplit(gb_record,split="\n")[[1]]
                cds_txt <- get_cds_text(gb_txt, c("coi","cox1"))
                if (is.na(cds_txt[1])) {
                    cat("Failed on gb_uid:",res$ids[i],"\n")
                    next
                }
                meta <- get_cds_meta(cds_txt)

                # Retrieve sequence and write to file
                seq <- entrez_fetch(db="nucleotide", id=res$ids[1], rettype="fasta")
                cat(file=seq_file, append=TRUE, sep="", seq)

                # Retrieve gb accession and taxonomic information from fasta header and save to data frame
                words <- strsplit(seq," ")[[1]]
                words[1] <- substring(words[1],first=2)
                if (words[3]=="sp.") words[3] <- "sp"

                # Extend data frame
                idx <- match(TipLabel,D$TipLabel)
                D$gb_accn[idx] <- words[1]
                D$gb_genus[idx] <- words[2]
                D$gb_species[idx] <- words[3]
                D$seq_start[idx] <- meta$seq_start
                D$seq_stop[idx] <- meta$seq_stop
                D$seq_len[idx] <- meta$seq_len
                D$codon_start[idx] <- meta$codon_start
                D$transl_table[idx] <- meta$transl_table
                D$mt_genome[idx] <- TRUE
                cat("Successfully fetched mitochondrial genome sequence for TipLabel:",TipLabel,"\n")
                success <- TRUE
                break
            }
        }
    }
    return (D)
}

# Define function for retrieving COI sequence from single sequence record
fetch_coi_sequence <- function(D, TipLabel, taxa) {

    success <- FALSE
    for (taxon in taxa) {

        if (success==TRUE)
            break;

        cat("Fetching coi sequence for taxon: ",taxon,"\n")

        res <- entrez_search(db="nucleotide", term=paste0("((coi[GENE] OR cox1[GENE]) AND ",taxon,"[ORGN]) NOT UNVERIFIED[TITL]"), retmax=50)

        if (res$count > 0) {

            # We cycle through all records returned and store the one with the highest alignment score
            seqs <- list()
            seqs_meta <- data.frame(gb_uid=numeric(),gb_accn=numeric())
            for (i in 1:length(res$ids)) {
                
                # Get cds info
                gb_record <- entrez_fetch(db="nucleotide", id=res$ids[i], rettype="gb", retmode="text")
                gb_txt <- strsplit(gb_record,split="\n")[[1]]
                cds_txt <- get_cds_text(gb_txt, c("coi","cox1"))
                # Retrieve sequence if it has a cds
                if (!is.na(cds_txt[1])) {
                    seq <- entrez_fetch(db="nucleotide", id=res$ids[i], rettype="fasta")
                    words <- strsplit(seq,split=" ")[[1]]
                    gb_accn <- substring(words[1],first=2)
                    seqs_meta <- rbind(seqs_meta,list(gb_uid=res$ids[i],gb_accn=gb_accn))
                    seqs <- c(seqs,seq)
                } else {
                    cat("Failed on gb_uid:",res$ids[i],"\n")
                }

            }
            
            # Score the sequences according to alignment match with COI
            if (length(seqs)==0) {   # No sequences with CDS
                cat("No sequences with CDS\n")
                next
            }
            scores <- get_aln_scores(seqs)
            if (nrow(scores)==0) { # No alignment match
                cat("No alignment match\n")
                next
            }

            # Identify the query sequence with the best score
            scores <- scores[order(-scores$aln_len, -scores$pct_id),]
            gb_uid <- seqs_meta$gb_uid[match(scores$query[1],seqs_meta$gb_accn)]

            # Get cds info (guaranteed to be available by code above)
            gb_record <- entrez_fetch(db="nucleotide", id=gb_uid, rettype="gb", retmode="text")
            gb_txt <- strsplit(gb_record,split="\n")[[1]]
            cds_txt <- get_cds_text(gb_txt, c("coi","cox1"))
            meta <- get_cds_meta(cds_txt)

            # Retrieve sequence and write to file
            seq <- entrez_fetch(db="nucleotide", id=gb_uid, rettype="fasta")
            cat(file=seq_file, append=TRUE, sep="", seq)

            # Retrieve gb accession and taxonomic information from fasta header and save to data frame
            words <- strsplit(seq," ")[[1]]
            words[1] <- substring(words[1],first=2)
            if (words[3]=="sp.") words[3] <- "sp"

            # Extend data frame
            idx <- match(TipLabel,D$TipLabel)
            D$gb_accn[idx] <- words[1]
            D$gb_genus[idx] <- words[2]
            D$gb_species[idx] <- words[3]
            D$seq_start[idx] <- meta$seq_start
            D$seq_stop[idx] <- meta$seq_stop
            D$seq_len[idx] <- meta$seq_len
            D$codon_start[idx] <- meta$codon_start
            D$transl_table[idx] <- meta$transl_table
            D$mt_genome[idx] <- FALSE
            cat("Successfully fetched coi sequence for tip label:",TipLabel,"\n")
            success <- TRUE
        }
    }

    return (D)
}


# Read in Cruaud et al tree
tree <- read.tree("IQ_COMBINED.tre")

# Read in base taxonomy for tips in the Cruaud et al tree
D <- read.delim("cruaud_taxonomy.tsv")

# Add columns for NCBI records
n <- nrow(D)
D$gb_accn <- character(n)
D$gb_genus <- character(n)
D$gb_species <- character(n)
D$seq_start <- numeric(n)
D$seq_stop <- numeric(n)
D$seq_len <- numeric(n)
D$codon_start <- numeric(n)
D$transl_table <- numeric(n)
D$mt_genome <- logical(n)

# Create fresh file for sequences
seq_file <- "cruaud_CO1.fasta"
cat(file=seq_file,"") # Make sure we start from scratch

# 1. Get all species-level sequences
for (i in 1:nrow(D)) {

    next

    if (D$Species_epithet[i] != "sp") {
        cat("Fetching species",D$Species[i],"\n")
        D <- fetch_sequence(D, D$TipLabel[i], D$Species[i])
    } else {
        sp_indet <- c("sp","n.a.","sp01","sp1","sp2","nsp018","nsp035","nsp045")
        AHE_species <- D$AHE_species[i][length(D$AHE_species[i])] # Get rid of prefixes like aff, cf, nr and the like
        UCE_species <- D$UCE_species[i][length(D$UCE_species[i])] # Ditto
        if (!(AHE_species %in% sp_indet)  || !(UCE_species %in% sp_indet)) {
            species <- character()
            if (!(AHE_species %in% sp_indet))
                species <- paste(D$AHE_genus[i],AHE_species)
            if (!(UCE_species %in% sp_indet))
                species <- c(species, paste(D$UCE_genus[i],UCE_species))
            cat("Fetching species",species,"\n")
            D <- fetch_sequence(D, D$TipLabel[i], species)
        }
    }
}

# 2. Get all singleton genus sequences
for (i in 1:nrow(D)) {
    
    next

    genus <- D$Genus[i]
    cat("Processing genus",genus,"\n")
    if (sum(D$Genus==genus)==1 && D$gb_accn[i] != "") {
        cat("Fetching genus",genus,"\n")
        D <- fetch_sequence(D,D$TipLabel[i], genus)
    }
}

# 3. Fill in missing genera if monophyletic
multi_genera <- unique(D$Genus[duplicated(D$Genus)])
for (genus in multi_genera[15:length(multi_genera)]) {

    cat("Processing multigenus",genus,"\n")
    idx <- which(D$Genus==genus)
    if (sum(D$gb_accn[idx]!="")==0 && is.monophyletic(tree,D$TipLabel[idx])) {
        cat("Fetching multigenus",genus,"\n")
        D <- fetch_sequence(D,D$TipLabel[idx[1]], genus)
    }
}

write.table(D,"cruaud_coi_extension_taxonomy.tsv", sep="\t",row.names=FALSE)
