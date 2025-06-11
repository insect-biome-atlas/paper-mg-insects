library(rentrez)

D <- read.delim("Chesters_Braconidae_taxonomy_updated.csv", sep=";")

E <- read.delim("uce_braconidae_genera_taxonomy.tsv")
E <- E[!duplicated(E$Genus),]   # Just in case

seq_file <- "missing_braconidae_genera_CO1.fasta"
file.remove(seq_file)   # Start from scratch

T <- data.frame(GenBank=character(), Subfamily=character(), Genus=character(), Species=character(), Start=numeric(), Stop=numeric())
P <- data.frame(Subfamily=character(), Genus=character())

for (i in 1:nrow(E)) {

    # We only want genera not present in Chesters but belonging to subfamilies represented there.
    # The genera of subfamilies missing in Chesters are added manually because they are mostly rare and not represented
    # by standard barcode sequences in GenBank.
    cat ("Now considering genus ", E$Genus[i], "\n", sep="")
    if (E$Genus[i] %in% D$Genus || !(E$Subfamily[i] %in% D$Subfamily_updated))
        next;

    cat ("Now analyzing genus ", E$Genus[i], "\n", sep="")
    genus <- E$Genus[i]

    res <- entrez_search(db="nucleotide", term=paste0("(COI[GENE] AND ",genus,"[ORGN] AND Braconidae[PORG] AND 658[SLEN]) NOT UNVERIFIED[TITL]"))
    print(res)

    if (res$count > 0 && grepl("codon_start=2", entrez_fetch(db="nucleotide", id=res$ids[1], rettype="gb"))) {

        # retrieve sequence and write to file
        seq <- entrez_fetch(db="nucleotide", id=res$ids[1], rettype="fasta")
        cat(file=seq_file, append=TRUE, sep="", seq)

        # retrieve taxonomic information and save to data frame
        words <- strsplit(seq," ")[[1]]
        words[1] <- substring(words[1],first=2)
        if (words[3]=="sp.") words[3] <- "sp"
        T <- rbind(T, list(GenBank=words[1], Subfamily=E$Subfamily[i], Genus=words[2], Species=words[3], Start=2, Stop=658))

    } else {

        # Check to see if there is potentially a non-standard CO1 sequence worth retrieving
        res <- entrez_search(db="nucleotide", term=paste0("(COI[GENE] AND ",genus,"[ORGN] AND Braconidae[PORG]) NOT UNVERIFIED[TITL]"))
        cat("Relaxed search\n")
        print(res)

        if (res$count > 0)
            P <- rbind(P, list(Subfamily=E$Subfamily[i], Genus=genus))
    }
}

# Fix some erroenous GenBank entries
idx <- which(T$GenBank=="JN288652.1")
T$Genus[idx] <- "Histeromerus"
T$Species[idx] <- "canadensis"
idx <- which(T$GenBank=="JN288574.1")
T$Genus[idx] <- "Neothlipsis"
T$Species[idx] <- "cincta"

write.table(T,"missing_braconidae_genera_taxonomy.tsv", sep="\t",row.names=FALSE)
write.table(P,"missing_braconidae_genera_check_needed.tsv", sep="\t",row.names=FALSE)

