# Functions for getting lca annotations from NCBI taxon IDs

# Function for generating taxon ids
get_taxon_ids <- function(ncbi_taxa) {

    ids <- character()
    elems <- strsplit(ncbi_taxa, split=";")
    for (i in 1:length(elems)) {
        ids <- c(ids,elems[[i]][1])
    }
    unique(ids)
}

# Function for generating lca results
lca <- function(T, taxon_ids) {

    X <- T[T$TaxonID %in% taxon_ids,]
    rank_names <- colnames(T)[10:3]
#    cat("rank_names = ",rank_names,"\n")
    res <- list()
    lca <- "unclassified"
    for (rank_name in rank_names) {
        a <- unique(X[,rank_name])
#        cat("Now looking at rank", rank_name, "and taxon/taxa", a, "\n")
        if (length(a) == 1) {
            if (a!="")
                lca <- a
            res <- c(res,list(a))
        }
        else {
            if (lca=="unclassified")
                res <- c(res,list("unclassified"))
            else
                res <- c(res,list(paste0("unclassified.",lca)))
        }
    }
    names(res) <- rank_names

    # Add name itself if it is ranked
    a <- unique(X[,1])
    if (length(a)==1 && a!="") {
        for (rank_name in rank_names) {
            if (res[rank_name] == "" && a %in% T[,rank_name]) {
                res[rank_name] <- a
                break
            }
        }
    }
    res
}

# Add lca information per row in input data
# from NCBI taxonomy file
add_lca_info <- function(B,T) {
    
    X <- data.frame()

    for (i in 1:nrow(B)) {

#        cat("Now considering this row:\n")
#        print(B[i,])
        taxon_ids <- get_taxon_ids(B$taxon.id[i])
        res <- lca(T, taxon_ids)
#        cat("lca fxn returned\n")
#        print(res)
        X <- rbind(X,res)
    }
    B$hit <- 1:nrow(B)
    X$hit <- 1:nrow(B)

    merge(B,X)
}

