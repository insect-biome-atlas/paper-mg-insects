library(ape)

# Define function that defines a taxset
add_taxset <- function(name, taxa, out_file) {

    # Allow also single-taxon taxset (a synonym of the taxon)
    if (length(taxa) >= 1) {
        cat("   taxset", name, "=", taxa, sep=" ", file=out_file, append=TRUE)
        cat(";\n", file=out_file, append=TRUE)
    }
}

# Define function that adds a hard constraint
add_hard_constraint <- function(name, ingroup, out_file, taxset=FALSE) {

    if (length(ingroup)>=2 || (taxset==TRUE && length(ingroup)>=1)) {
        cat ("   constraint ", name, " = ", sep="", file=out_file, append=TRUE)
        cat (ingroup, sep=" ", file=out_file, append=TRUE)
        cat (";\n", file=out_file, append=TRUE)
    }
}

# Define function that adds a partial constraint
add_soft_constraint <- function(name, ingroup, outgroup, out_file, taxset=FALSE) {

    if ((length(ingroup)>=2 && length(outgroup)!=0) || (taxset==TRUE && length(ingroup)>=1 && length(outgroup)>=1)) {
        cat ("   constraint ", name, " partial = ", sep="", file=out_file, append=TRUE)
        cat (ingroup, sep=" ", file=out_file, append=TRUE)
        cat (" : ", file=out_file, append=TRUE)
        cat (outgroup, sep=" ", file=out_file, append=TRUE)
        cat (";\n", file=out_file, append=TRUE)
    }
}

# Define subset partial constraint function
add_partial_subset_constraint <- function(name, ingroup, alltaxa, out_file) {

    outgroup <- alltaxa[!(alltaxa %in% ingroup)]

    add_soft_constraint(name, ingroup, outgroup, out_file)
}

# Define function that converts a tree to partial constraints
gen_mb_con_file <- function(tree_file=tree_file, out_file=out_file) {

    # Read in the tree
    T <- read.tree(tree_file)

    if (!is.rooted(T)) {
        cat ("ERROR: Input tree is not rooted")
    }

    # Get basic numbers
    num_nodes <- T$Nnode
    num_tips <- length(T$tip.label)
    
    # Print header to output file
    cat("#NEXUS\n\nbegin mrbayes;\n",file=out_file)

    # Loop over clades and extract tips recursively. We do
    # not rely on any particular order of nodes here (O(n*n)
    # complexity insteaqd of O(n))
    for (node in (num_tips+1):(num_tips+num_nodes)) {
        
        tips <- leaves(T, node)
 
        # Output informative constraint partitions
        if (length(tips)<num_tips) {
            cat ("   constraint node", node, " partial = ", sep="", file=out_file, append=TRUE)
            cat (T$tip.label[tips], sep=" ", file=out_file, append=TRUE)
            cat (" : ", file=out_file, append=TRUE)
            cat (T$tip.label[-tips], sep=" ", file=out_file, append=TRUE)
            cat (";\n", file=out_file, append=TRUE)
        }
    }

    # Print tail to output file
    cat ("end;\n", file=out_file, append=TRUE)
}


# Recursive function to get leaf/tip indices of a node
# from a tree in ape format
leaves <- function(T, node) {
    if (node <= length(T$tip.label))
        return (node);
    
    desc_rows <- which(T$edge[,1]==node)
    left  <- T$edge[desc_rows[1],2]
    right <- T$edge[desc_rows[2],2]
    return ( c(leaves(T, left), leaves(T,right)) )
}


