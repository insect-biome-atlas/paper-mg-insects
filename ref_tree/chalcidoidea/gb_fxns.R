# Functions for retrieving sequences from genbank

library(ape)
library(rentrez)
library(stringr)
library(Biostrings)

# Build a reference database for the searches
# Here we are using a 652 bp barcode from Torymus,
# which is representative of the standard barcode
# for most chalcidoids.
system("vsearch --makeudb_usearch reference_barcode_Torymus_652.fasta --output ref.udb")

get_elem_value <- function(elem_name,txt) {
   
    idx <- which(grepl(elem_name,txt))
    if (length(idx)==0) {
        cat("WARNING: Could not find '",elem_name,"' in record\n")
        return (NA)
    }
    b <- strsplit(txt[idx[1]],split=paste0(elem_name,"="))[[1]]
    return (as.numeric(b[length(b)]))
}

get_elem_text <- function(elem_name,txt) {

    idx <- which(grepl(elem_name,txt))
    if (length(idx)==0) {
        cat("WARNING: Could not find '",elem_name,"' in record\n")
        return (NA)
    }
    b <- strsplit(txt[idx[1]],split=paste0(elem_name,"="))[[1]]
    return (b[length(b)])
}

get_cds_text <- function(gb_txt, gene_names) {
   
    for (i in 1:length(gene_names))
        gene_names[i] <- paste0('\"',tolower(gene_names[i]),'\"')

    idx <- which(grepl("/gene=",gb_txt))
    idx1 <- numeric()
    for (i in 1:length(idx)) {
        a <- get_elem_text("/gene", gb_txt[idx[i]])
        if (tolower(a) %in% gene_names)
            idx1 <- c(idx1,idx[i])
    }

    if (length(idx1) == 0) {
        cat("WARNING: No feature for this gene\n")
        return(character())
    }
    if (length(idx1) > 2) {
        cat("WARNING: More than one feature for this gene\n")
        return(character())
    }

    idx2 <- which(grepl("CDS",gb_txt[idx1[1]:length(gb_txt)])) + idx1[1] - 1
    if (length(idx2)==0) {
        cat("WARNING: No CDS for this gene\n")
        return(character())
    }

    idx3 <- which(grepl("/gene",gb_txt[idx2[1]:length(gb_txt)])) + idx2[1] - 1
    if (length(idx3)==0 || idx3[1] != idx1[2]) {
        cat("WARNING: No matching CDS for this gene\n")
        return(character())
    }

    return (gb_txt[idx2[1]:length(gb_txt)])
}

get_start_stop <- function(cds_line) {

    # Find the start and end specification
    a <- strsplit(cds_line,split=" ")[[1]]
    b <- a[length(a)]   # Pick last element

    cds_complement <- FALSE
    # Remove extraneous characters    
    if (grepl("complement",b)) {
        b <- gsub("complement(","",b,fixed=TRUE)
        b <- gsub(")","",b,fixed=TRUE)
        cds_complement <- TRUE
    }
    if (grepl("<",b)) {
        b <- gsub("<","",b,fixed=TRUE)
        b <- gsub(">","",b,fixed=TRUE)
    }

    # Split and save start and stop
    c <- strsplit(b,split="..",fixed=TRUE)[[1]]
    cds_start <- as.numeric(c[1])
    cds_end <- as.numeric(c[2])

    list(cds_complement=cds_complement, cds_start=cds_start, cds_end=cds_end, cds_len=cds_end-cds_start+1)
}

get_cds_meta <- function(gb_txt, gene_names=character()) {

    res <- list(cds_info=FALSE,
                cds_complement=NA,
                cds_start=NA,
                cds_end=NA,
                cds_len=NA,
                cds_codon_start=NA,
                cds_transl_table=NA
                )

    if (length(gb_txt)==0) 
        return (data.frame())

    if (length(gene_names)==0) {
        cat("ERROR: Missing gene names\n")
        return (NULL)
    }

    cds_text <- get_cds_text(gb_txt, gene_names)

    if (length(cds_text) > 0) {

        x <- get_start_stop(cds_text[1])
        codon_start <- get_elem_value("codon_start",cds_text)
        transl_table <- get_elem_value("transl_table",cds_text)
        res <- c(list(cds_info=TRUE), x, list(cds_codon_start=codon_start, cds_transl_table=transl_table))
    }

    return (res)
}

get_aln_scores <- function(seqs) {

    cat(file="query.fasta","")
    for(i in 1:length(seqs)) {
        cat(file="query.fasta", seqs[[i]], sep="\n",append=TRUE)
    }

    system("vsearch --usearch_global query.fasta --db ref.udb --id 0.60 --strand both --blast6out vsearch_res.tsv --alnout vsearch_aln.txt")

    # Retrieve the blast-type results
    res <- read.table("vsearch_res.tsv", header=FALSE,
                        col.names=c("aln_query","aln_target","aln_pct_id","aln_len","aln_mis","aln_gaps",
                        "aln_q_start","aln_q_end","aln_t_start","aln_t_end","aln_evalue","aln_bits"))

    # Get start and end of the local alignment (only global values returned by --usearch_global)
    # but only if some alignments succeeded
    if (nrow(res) > 0) {
        aln_data <- get_aln_info("vsearch_aln.txt")
        res$aln_q_start <- aln_data$q_start
        res$aln_q_end <- aln_data$q_end
        res$aln_t_start <- aln_data$t_start
        res$aln_t_end <- aln_data$t_end
        res$aln_q_neg <- aln_data$q_neg
    } else {
        res <- data.frame(aln_query="",
                          aln_target="",
                          aln_pct_id=0,
                          aln_len=0,
                          aln_mis=0,
                          aln_gaps=0,
                          aln_q_start=0,
                          aln_q_end=0,
                          aln_t_start=0,
                          aln_t_end=0,
                          aln_evalue=0,
                          aln_bits=0,
                          aln_q_neg=NA
        )
    }

    return (res)
}

get_aln_info <- function(aln_file) {
    
    aln_data <- data.frame()

    aln <- readLines(aln_file)
    aln_end_lines <- which(grepl("cols",aln))

    sub_aln_start <- 1
    for (sub_aln_end in aln_end_lines) {

        sub_aln <- aln[sub_aln_start:sub_aln_end]
        
        qry <- sub_aln[grepl("Qry",sub_aln)]
        tgt <- sub_aln[grepl("Tgt",sub_aln)]
        qry <- str_squish(qry)
        tgt <- str_squish(tgt)
        
        words <- strsplit(qry[1]," ")[[1]]
        q_start <- as.numeric(words[2])
        q_neg <- words[3]=="-"
        words <- strsplit(qry[length(qry)]," ")[[1]]
        q_end <- as.numeric(words[length(words)])
    
        words <- strsplit(tgt[1]," ")[[1]]
        t_start <- as.numeric(words[2])
        words <- strsplit(tgt[length(tgt)]," ")[[1]]
        t_end <- as.numeric(words[length(words)])

        aln_data <- rbind(aln_data,list(q_start=q_start,q_end=q_end,q_neg=q_neg,t_start=t_start,t_end=t_end))

        sub_aln_start <- sub_aln_end+1
    }

    return (aln_data)
}

# Define function for fetching coi sequences
fetch_coi_sequences <- function(taxa, seq_file) {

    df <- data.frame()
    if (length(taxa)==0)
        return (df)

    for (i in 1:length(taxa)) {

        taxon <- taxa[i]
        cat("Fetching mt coi sequences for taxon:",taxon,"\n")

        res <- entrez_search(db="nucleotide", term=paste0("(mitochondrion[TITL] AND ",taxon,"[ORGN]) NOT UNVERIFIED[TITL]"))
        df <- rbind(df, fetch_sequences(res$ids, gene_names=c("coi","cox1"), seq_file=seq_file))
    }

    for (i in 1:length(taxa)) {

        taxon <- taxa[i]
        cat("Fetching single coi sequences for taxon:",taxon,"\n")

        res <- entrez_search(db="nucleotide", term=paste0("((coi[GENE] OR cox1[GENE]) AND ",taxon,"[ORGN]) NOT UNVERIFIED[TITL]"), retmax=50)
        df <- rbind(df, fetch_sequences(res$ids, gene_names=c("coi","cox1"), seq_file=seq_file))
    }

    return (df)
}

# Define function for fetching all coi sequences
fetch_all_coi_sequences <- function(taxa, seq_file, P) {

    df <- data.frame()
    if (length(taxa)==0)
        return (df)

    for (i in 1:length(taxa)) {

        taxon <- taxa[i]
        cat("Fetching ALL single coi sequences for taxon:",taxon,"\n")

        res <- entrez_search(db="nucleotide", term=paste0("((coi[GENE] OR cox1[GENE]) AND ",taxon,"[ORGN]) NOT UNVERIFIED[TITL]"), retmax=50)
        retmax <- res$count
        res <- entrez_search(db="nucleotide", term=paste0("((coi[GENE] OR cox1[GENE]) AND ",taxon,"[ORGN]) NOT UNVERIFIED[TITL]"), retmax=retmax)

        gb_uids <- res$ids[!res$ids %in% P$gb_uids]
        df <- rbind(df, fetch_sequences(gb_uids, gene_names=c("coi","cox1"), seq_file=seq_file))
    }

    return (df)
}

# Define function for fetching cytochrome oxidase sequences
fetch_cytochrome_oxidase_sequences <- function(taxa, seq_file, P) {

    df <- data.frame()
    if (length(taxa)==0)
        return (df)

    for (i in 1:length(taxa)) {

        taxon <- taxa[i]
        cat("Fetching cytochrome oxidase sequences for taxon:",taxon,"\n")

        res <- entrez_search(db="nucleotide", term=paste0("(\"cytochrome oxidase\"[TITL] AND ",taxon,"[ORGN]) NOT UNVERIFIED[TITL] NOT (coi[GENE] OR cox1[GENE])"))
        
        gb_uids <- res$ids[!res$ids %in% P$gb_uids]
        df <- rbind(df, fetch_sequences(gb_uids, gene_names=c("coi","cox1"), seq_file=seq_file))
    }

    return (df)
}

# Define function for fetching specific sequence accession numbers
fetch_coi_gb_accsns <- function(gb_accsns, seq_file) {

    df <- data.frame()
    if (length(gb_accsn)==0)
        return (df)

    for (i in 1:length(gb_accsn)) {

        gb_accsn <- gb_accsns[i]
        cat("Fetching genbank accsn:",gb_accsn,"\n")

        res <- entrez_search(db="nucleotide", term=paste0(gb_accsn,"[ACCN]"))
        df <- rbind(df, fetch_sequences(res$ids, gene_names=c("coi","cox1"), seq_file=seq_file))
    }

    return (df)
}

# Define function for fetching sequences
fetch_sequences <- function(gb_uids, gene_names, seq_file) {

    # Create empty data frame for results
    seq_meta <- data.frame()

    # Create empty list for sequences
    seqs <- list()

    if (length(gb_uids) > 0) {

        # Cycle over results
        for (i in 1:length(gb_uids)) {
                
            # Get cds info (if any)
            gb_record <- entrez_fetch(db="nucleotide", id=gb_uids[i], rettype="gb", retmode="text")
            gb_txt <- strsplit(gb_record,split="\n")[[1]]
            cds_meta <- get_cds_meta(gb_txt, gene_names)

            # Retrieve sequence and write to file
            seq <- entrez_fetch(db="nucleotide", id=gb_uids[i], rettype="fasta")
            cat(file=seq_file, append=TRUE, sep="", seq)
            seqs <- c(seqs,seq)

            # Retrieve gb accession and taxonomic information from fasta header 
            words <- strsplit(seq," ")[[1]]
            words[1] <- substring(words[1],first=2)
            if (words[3]=="sp.") words[3] <- "sp"

            # Add to results
            seq_meta <- rbind(seq_meta,
                              c(list(gb_uid=gb_uids[i],
                                     gb_accn=words[1],
                                     gb_genus=words[2],
                                     gb_species=words[3]
                                     ),
                                cds_meta
                               )
                             )
        }

        # Get alignment scores
        aln_scores <- get_aln_scores(seqs)
        if (nrow(aln_scores > 0)) {
            aln_scores$gb_accn <- aln_scores$aln_query
            seq_meta <- merge(seq_meta, aln_scores, by="gb_accn", all.x=TRUE)
        }
    }

    return (seq_meta)
}

