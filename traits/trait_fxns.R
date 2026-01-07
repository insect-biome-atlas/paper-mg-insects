# Functions for processing trait data

# Functions for checking the coding
# =================================

check_coding <- function(D,key="Clade") {
    
    cat("Checking the coding\n")
    E <- D[,c(key,"Niche","Habitat")]

    if (sum(duplicated(E[,key])) != 0) {
        cat("WARNING: Inconsistent entries:\n")
        dup_key <- E[,key][duplicated(E[,key])]
        print(E[E[,key] %in% dup_key,])
    } else {
        cat("Passed the check\n")
    }
}

check_family <- function(D) {

    cat("Checking the uniqueness of Family entries\n")
    if (sum(!grepl("idae",D$Family))!=0) {
        cat("WARNING: Non-family entries:\n")
        print(D[!grepl("idae",D$Family),])
    } else {
        cat("Passed the check\n")
    }
}

# Functions deriving community and trophic level info
community <- function(x) {
    if (grepl("Predator",x))
        return ("Predator")
    else if (grepl("Saprophage",x))
        return ("Saprophage")
    else if (grepl("Phytophage",x))
        return ("Phytophage")
    else if (grepl("Parasite",x))
        return ("Homeothermic_vertebrate")
    else return ("Unknown")
}

level <- function(x) {
    if (grepl("Predator-parasitoid",x) | grepl("Parasite",x))
        return ("Tertiary")
    else if (grepl("-parasitoid",x))
        return ("Secondary")
    else if (x=="")
        return ("Unknown")
    else
        return ("Primary")
}

