# Script for generating corrected Chesters tree

# We start from the Insecta tree presented by Misof
# et al. (2014, Science) and use this to constrain the
# interordinal relationships in Insecta. We extract
# each order subtree from the ultrametric tree in Chesters
# (2017), and glue it onto the order-level backbone.

# Specifically, we start with the version of the Chesters
# tree presented by  Sundh et al (2024), in which tips
# with poor-quality sequences have been removed. This
# version of the tree also has the original outgroups
# removed (since they are sprinkled over the tree), and
# replaced by an expanded outgroup tree including more
# appropriate taxa for phylogenetic placement of hexapod
# sequences.

# We use a modified version of this tree, in which some
# erroneous sequences have been removed using the script
# 'correct_sundh_et_al.R' contained in this directory.
# The corrected tree is in the file named
# 'chesters_new_outgroups_updated.nwk".

# We remove Psocoptera + Phthiraptera (=Psocodea) as one
# clade from the original tree, as Psocoptera is not
# monophyletic with respect to Phthiraptera.

# We then start from the Misof et al (2014) backbone tree
# for Insecta, extract each order from the Sundh et al (2024)
# tree and reinsert it where it belongs in the backbone
# tree.

# We do not add the outgroup part to the resulting tree
# as we do that when expanding the resulting tree.

# Load APE
library(ape)

# Read in Misof et al backbone tree for Insecta
M14 <- read.tree("misof_2014_fig1_insecta_orders_dated.nwk")

# Read in Sundh et al version of Chesters 2017 tree
# with erroneous or problematic sequences removed.
S24 <- read.tree("chesters_new_outgroups_updated.nwk")

# Read in the corrected taxonomy data for Sundh et al
D <- read.delim("chesters_new_outgroups_taxonomy_updated.tsv")

# Psocodea: combine Psocoptera and Phthiraptera
# in the S24 tree.
D$Order[D$Order %in% c("Phthiraptera","Psocoptera")] <- "Psocodea"

# Neuropterida: combine Neuroptera, Raphidioptera and
# Megaloptera in S24 tree.
D$Order[D$Order %in% c("Neuroptera","Megaloptera","Raphidioptera")] <- "Neuropterida"

# Define vector of orders (using Neuropterida and Psocodea)
simple_orders <- c(
    "Archaeognatha",
    "Zygentoma",
    "Ephemeroptera",
    "Odonata",
    "Dermaptera",
    "Plecoptera",
    "Orthoptera",
    "Blattodea",
    "Mantodea",
    "Phasmatodea",
    "Psocodea",
    "Thysanoptera",
    "Hemiptera",
    "Hymenoptera",
    "Lepidoptera",
    "Trichoptera",
    "Diptera",
    "Mecoptera",
    "Siphonaptera",
    "Strepsiptera",
    "Coleoptera",
    "Neuropterida",
    "Zoraptera",
    "Embioptera",
    "Grylloblattodea",
    "Mantophasmatodea")

# Define a function for extracting a clade from the S24
# tree and pasting it into the M14-derived input tree
replace_subtree <- function(order, tree) {

    # Special case for just one taxon
    tips <- D$TipLabel[D$Order==order]
    if (length(tips) == 1) {
        tree$tip.label[tree$tip.label==order] <- tips
        return (tree)
    }

    # Extract donor tree and compute height
    donor_tree <- keep.tip(S24,tips)
    donor_height <- node.depth.edgelength(donor_tree)[1]     # Distance from root node to first tip

    # Get crown node and stalk length from source tree
    crown_node <- getMRCA(S24,donor_tree$tip.label)
    stalk_length <- S24$edge.length[which(S24$edge[,2]==crown_node)]

    # Identify insert position in receptor tree
    insert_tip <- which(tree$tip.label==order)
    insert_height <- tree$edge.length[which(tree$edge[,2]==insert_tip)]

    # Adjust scaling to receptor tree
    scale_factor <- insert_height / (donor_height + stalk_length)
    donor_tree$edge.length <- donor_tree$edge.length * scale_factor
    
    # Add a half-length stalk (the rest will be added when the insert_tip is removed)
    donor_tree$root.edge <- scale_factor * 0.5 * stalk_length
    temp_pos <- scale_factor * (donor_height + 0.5 * stalk_length)

    # Join trees
    new_tree <- bind.tree(tree, donor_tree, where=insert_tip, position=temp_pos)
    
    # Drop the insert_tip and return tree
    drop.tip(new_tree, insert_tip)
}


new_tree <- M14
for (ord in simple_orders) {
    cat("Processing",ord,"\n")
    new_tree <- replace_subtree(ord, new_tree)
}

# Write the resulting tree
write.tree(new_tree,"s24insecta_m14backbone.nwk")

# Write the taxonomic info for this tree
D <- D[D$Class=="Insecta",]
write.table(D,"s24insecta_m14backbone_taxonomy.tsv",sep="\t",row.names=FALSE)

