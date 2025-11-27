# Define function for selecting the best coi reference sequence
select_best_ref_seq <- function(df, min_barcode_overlap=412, pct_id_threshold=70) {
    
    df <- df[!is.na(df$aln_query),]
    if (nrow(df)==0)
        return ("")

    df <- df[df$aln_pct_id > pct_id_threshold,]
    if (nrow(df)==0)
        return ("")

    df[order(df$aln_len,-df$aln_gaps),]

    X <- df[df$aln_len==652 & df$aln_gaps==0,]
    if (nrow(X) > 0)
        return (X$gb_accn[1])

    df$metabarcode_aln_end <- 412 - (652 - df$aln_t_end)
    df$metabarcode_aln_start <- numeric(nrow(df))
    for (i in 1:nrow(df)) { df$metabarcode_aln_start[i] <- max(1, df$aln_t_start[i] - 240) }
    df$metabarcode_aln_len <- df$metabarcode_aln_end - df$metabarcode_aln_start + 1

    X <- df[df$metabarcode_aln_len >= min_barcode_overlap,]
    if (nrow(X) > 0)
        return (X$gb_accn[1])

    return ("")
}

