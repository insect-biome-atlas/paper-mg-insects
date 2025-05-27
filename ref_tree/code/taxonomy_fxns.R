library(taxize)

get_subfam_supertribe <- function(sci_id) {

    subfam <- character(length(sci_id))
    supertribe <- character(length(sci_id))

    x <- classification(sci_id=sci_id,db="ncbi")

    for (i in 1:length(sci_id)) {
        idx <- which(x[[i]]$rank=="subfamily")
        subfam[i] <- x[[i]]$name[idx]
        supertribe[i] <- x[[i]]$name[idx+1]
    }

    supertribe[!grepl("idi",supertribe)] <- ""
    
    list(Subfamily=subfam, Supertribe=supertribe)
}

get_subfam <- function(sci_id) {

    subfam <- character(length(sci_id))

    x <- classification(sci_id=sci_id,db="ncbi")

    for (i in 1:length(sci_id)) {
        idx <- which(x[[i]]$rank=="subfamily")
        subfam[i] <- x[[i]]$name[idx]
    }

    list(Subfamily=subfam)
}

