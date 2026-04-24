# Script for adding age estimates to placements

library(ape)

# Read in reference tree
tree <- read.tree("../../ref_tree/expanded_tree/chesters_expanded.nwk")

# Read in placement tree with edge numbers in tip labels and node labels
# NB! The root has been removed in this tree
ptree <- read.tree("../source/mg_placement_tree_edgenum.nwk")

# Read in placements (only includes LWR > 0.5 placements)
P <- readRDS("../source/mg_cluster_placements.rds")

# Read in cluster taxonomy
T <- read.delim("../../iba_data/cluster_taxonomy_mg.tsv")

# Restrict placements to those in cluster taxonomy
P <- P[P$cluster_rep %in% T$ASV,]

# Remove the source column
idx <- which(colnames(P)=="source")
P <- P[,-idx]

# Get the ages of the nodes and leaves in the reference tree (in Ma)
y <- node.depth.edgelength(tree)
ages <- round(y[1] - y, 2)

# Translate from ptree (edge_num) to ptree node index
# We need to traverse the tree to do this safely, as
# the postorder traversal does not happen in the same
# way in ape as in the algorithm used for the placement
# tree, referred to in the epa-ng placement data.
# We use a recursive postorder algorithm that maps
# edge numbers (in node.label or tip.label) to
# ape tree node numbers.
edge_nums2nodes <- function(tree, node, nodes) {
    
    if (node > length(tree$tip.label)) {
        idx <- which(tree$edge[,1]==node)
        for (i in idx)
            nodes <- edge_nums2nodes(tree, tree$edge[i,2], nodes)
        lbl <- tree$node.label[node-length(tree$tip.label)]
        if (lbl!="NA") {
            edge_idx <- as.numeric(lbl) + 1 # Edge numbers are 0 offset
            nodes[edge_idx] <- node
        }
    } else {
        edge_idx <- as.numeric(tree$tip.label[node]) + 1   # Edge numbers are 0 offset
        nodes[edge_idx] <- node
    }
    return (nodes)
}


# Get edge length
get_edge_length <- function(node, tree) {
    idx <- which(tree$edge[,2]==node)
    return (tree$edge.length[idx])
}

# Translate from ptree (unrooted) to tree (rooted)
ptree_node2tree_node <- function(x, tree) {

    if (x <= length(tree$tip.label)) {
        return (x)
    } else { return (x+1) }
}

cat("Computing placement index to node index mapping\n")
nodes <- numeric(length(ptree$edge.length))
nodes <- edge_nums2nodes(ptree, length(ptree$tip.label)+1, nodes)
P$ptree_node <-numeric(nrow(P))
P$ptree_node <- nodes[P$edge_num+1]
cat("Getting edge lengths\n")
P$edge_length <- numeric(nrow(P))
for (i in 1:nrow(P)) { P$edge_length[i] <- get_edge_length(P$ptree_node[i], ptree) }
cat("Getting placement tree to reference tree node mapping\n")
P$tree_node <-numeric(nrow(P))
for (i in 1:nrow(P)) { P$tree_node[i] <- ptree_node2tree_node(P$ptree_node[i], tree) }
cat("Getting placement ages\n")
P$child_age <- ages[P$tree_node]
for (i in 1:nrow(P)) { P$parent_age[i] <- ages[tree$edge[which(tree$edge[,2]==P$tree_node[i]),1]] }
P$placement_age <- (P$distal_length / P$edge_length) * (P$parent_age - P$child_age) + P$child_age
saveRDS(P,"../data/placement_stats.rds")

