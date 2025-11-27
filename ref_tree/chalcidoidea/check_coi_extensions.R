source("select_ref_seq_fxn.R")

check_coi_extensions  <- function(D,E) {

    misses <- 0
    for (i in 1:nrow(D)) {

        if (D$gb_accn[i]!="")
            next

        F <- E[E$TipLabel==D$TipLabel[i],]

        gb_accn <- select_best_ref_seq(F,300,60)

        if (gb_accn!="") {
            misses <- misses + 1
            cat("For TipLabel", D$TipLabel[i], "there is a gb_accn", gb_accn, "that passes the criteria\n")
        }
    }

    if (misses==0)
        cat("Found", misses, " sequences that could have been included in the coi sequences assembled according to criteria\n")
}

