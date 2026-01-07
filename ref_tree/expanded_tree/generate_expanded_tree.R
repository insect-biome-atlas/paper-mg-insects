# Script for expanding the reference tree for selected
# groups

library(ape)

# Define function for replacing a clade
replace_clade <- function(tree, ingroup, new_tree, new_ingroup) {

#    print(str(tree))
#    print(str(ingroup))
#    print(str(new_tree))
#    print(str(new_ingroup))

#    cat("Getting nodes in tree\n")
    if (!is.monophyletic(tree, ingroup)) {
        cat("ERROR: Old ingroup is not monophyletic, expanded tree is unlikely to be correct\n")
        return (NULL);
    }
    if (length(ingroup)==1)
        crown_node <- which(tree$tip.label==ingroup)
    else
        crown_node <- getMRCA(tree,ingroup)
#    cat("crown_node",crown_node,"\n")
    stem_node <- tree$edge[which(tree$edge[,2]==crown_node),1]
#    cat("stem_node",stem_node,"\n")
    desc <- tree$edge[which(tree$edge[,1]==stem_node),2]
#    cat("desc",desc,"\n")
    sister_node <- desc[which(desc!=crown_node)]
#    cat("sister_node",sister_node,"\n")
    if (sister_node <= length(tree$tip.label))
        sistergroup <- tree$tip.label[sister_node]
    else
        sistergroup <- extract.clade(tree,sister_node)$tip.label
#    print(str(sistergroup))

    depths <- node.depth.edgelength(tree)
    stem_height <- depths[1] - depths[stem_node]
    insert_pos <- tree$edge.length[which(tree$edge[,2]==sister_node)]
#    cat("insert_pos",insert_pos, "\n")

    tree <- drop.tip(tree,ingroup)

#    cat("Getting nodes in new_tree\n")
    if (length(new_ingroup)==1)
        new_crown_node <- which(new_tree$tip.label==new_ingroup)
    else
        new_crown_node <- getMRCA(new_tree,new_ingroup)
#    cat("New crown node:",new_crown_node,"\n")
    new_stem_node <- new_tree$edge[which(new_tree$edge[,2]==new_crown_node),1]
#    cat("New stem node:",new_stem_node,"\n")

    new_depths <- node.depth.edgelength(new_tree)
    new_stem_height <- new_depths[1] - new_depths[new_stem_node]
    new_tree_stalk <- new_depths[new_crown_node] - new_depths[new_stem_node]

#    cat("Rescaling new_tree\n")
    scale_factor <- stem_height / new_stem_height
    new_tree$edge.length <- new_tree$edge.length * scale_factor

#    cat("Dropping tips in new_tree\n")
    new_tree <- keep.tip(new_tree, new_ingroup) # This also drops the "stalk", so the stalk needs to be reinserted
    if (length(new_ingroup) > 1)  
      new_tree$root.edge <- new_tree_stalk * scale_factor
    else
      new_tree$root.edge <- 0.0   # Avoid adding stalk twice if new tree has only one taxon
    
#    cat("Binding trees\n")
    if (length(sistergroup)==1)
        bind.tree(tree, new_tree, where=which(tree$tip.label==sistergroup), position=insert_pos)
    else
        bind.tree(tree, new_tree, where=getMRCA(tree,sistergroup), position=insert_pos)
}

# Define function for updating taxonomy
update_taxonomy <- function(taxonomy, ingroup, new_taxonomy, new_ingroup, cols) {

    rbind(taxonomy[!(taxonomy$TipLabel %in% ingroup),cols], new_taxonomy[new_taxonomy$TipLabel %in% new_ingroup,cols])
}

# Define function for updating sequence file
update_sequences <- function(seqs, ingroup, new_seqs, new_ingroup) {

    c(seqs[!names(seqs) %in% ingroup], new_seqs[names(new_seqs) %in% new_ingroup])
}


# Read in start tree, taxonomy and sequences
cat("Reading in starting tree\n")
expanded_tree <- read.tree("../chesters_2017/s24insecta_m14backbone.nwk")
expanded_taxonomy <- read.delim("../chesters_2017/s24insecta_m14backbone_taxonomy.tsv")
expanded_seqs <- read.FASTA("../chesters_2017/chesters_new_outgroups.fasta")
expanded_seqs <- expanded_seqs[names(expanded_seqs) %in% expanded_tree$tip.label]
if (length(expanded_seqs) != length(expanded_tree$tip.label))
    cat("ERROR: TipLabel mismatch between sequence file and tree\n")

# Add in Clade (defaults to same as Family)
expanded_taxonomy$Clade <- expanded_taxonomy$Family

# Define the columns we wish to keep in the taxonomy
cols <- c("TipLabel","Kingdom","Phylum","Class","Order","Family","Clade","Genus","Species")
expanded_taxonomy <- expanded_taxonomy[,cols]


# Replace Ichneumonidae
# ---------------------

cat("Replacing Ichneumonidae\n")

# Read new ichneumonid tree (pick the last tree as a sample from the posterior)
trees <- read.nexus("../ichneumonidae/mb_runs/run1/tree_sample.tre")
ichneumonidae_tree <- trees[[length(trees)]]

# Read in new ichneumonid taxonomy and sequences
ichneumonidae_taxonomy <- read.delim("../ichneumonidae/expanded_ichneumonidae_taxonomy.tsv")
ichneumonidae_seqs <- read.FASTA("../ichneumonidae/expanded_ichneumonidae.fasta")

# Get respective ingroup taxa
ingroup <- expanded_taxonomy$TipLabel[expanded_taxonomy$Family=="Ichneumonidae"]
new_ingroup <- ichneumonidae_taxonomy$TipLabel[ichneumonidae_taxonomy$Family=="Ichneumonidae"]

# Update tree, taxonomy and sequences
expanded_tree <- replace_clade(expanded_tree, ingroup, ichneumonidae_tree, new_ingroup)
expanded_taxonomy <- update_taxonomy(expanded_taxonomy, ingroup, ichneumonidae_taxonomy, new_ingroup, cols)
expanded_seqs <- update_sequences(expanded_seqs, ingroup, ichneumonidae_seqs, new_ingroup)


# Replace Braconidae
# ------------------

cat("Replacing Braconidae\n")

# Read new braconid tree (pick the last tree as a sample from the posterior)
trees <- read.nexus("../braconidae/mb_runs/local_run1/tree_sample.tre")
braconidae_tree <- trees[[length(trees)]]

# Read in new braconid taxonomy and sequences
braconidae_taxonomy <- read.delim("../braconidae/expanded_braconidae_taxonomy.tsv")
braconidae_seqs <- read.FASTA("../braconidae/expanded_braconidae.fasta")

# Get respective ingroup taxa
ingroup <- expanded_taxonomy$TipLabel[expanded_taxonomy$Family=="Braconidae"]
new_ingroup <- braconidae_taxonomy$TipLabel[braconidae_taxonomy$Family=="Braconidae"]

# Update tree, taxonomy and sequences
expanded_tree <- replace_clade(expanded_tree, ingroup, braconidae_tree, new_ingroup)
expanded_taxonomy <- update_taxonomy(expanded_taxonomy, ingroup, braconidae_taxonomy, new_ingroup, cols)
expanded_seqs <- update_sequences(expanded_seqs, ingroup, braconidae_seqs, new_ingroup)


# Replace Cecidomyiidae
# ---------------------

cat("Replacing Cecidomyiidae\n")

# Read new cecidomyiid tree (pick the last tree as a sample from the posterior)
trees <- read.nexus("../cecidomyiidae/mb_runs/run1/tree_sample.tre")
cecidomyiidae_tree <- trees[[length(trees)]]

# Read in new cecidomyiid taxonomy and sequences
cecidomyiidae_taxonomy <- read.delim("../cecidomyiidae/expanded_cecidomyiidae_taxonomy.tsv")
cecidomyiidae_seqs <- read.FASTA("../cecidomyiidae/expanded_cecidomyiidae.fasta")

# Get respective ingroup taxa
ingroup <- expanded_taxonomy$TipLabel[expanded_taxonomy$Family=="Cecidomyiidae"]
new_ingroup <- cecidomyiidae_taxonomy$TipLabel[cecidomyiidae_taxonomy$Family=="Cecidomyiidae"]

# Update tree, taxonomy and sequences
expanded_tree <- replace_clade(expanded_tree, ingroup, cecidomyiidae_tree, new_ingroup)
expanded_taxonomy <- update_taxonomy(expanded_taxonomy, ingroup, cecidomyiidae_taxonomy, new_ingroup, cols)
expanded_seqs <- update_sequences(expanded_seqs, ingroup, cecidomyiidae_seqs, new_ingroup)


# Replace Zoraptera
# -----------------

cat("Replacing Zoraptera\n")

# Read new zorapteran tree (pick the last tree as a sample from the posterior)
trees <- read.nexus("../zoraptera/mb_runs/run1/tree_sample.tre")
zoraptera_tree <- trees[[length(trees)]]

# Read in new zorapteran taxonomy and sequences
zoraptera_taxonomy <- read.delim("../zoraptera/zoraptera_taxonomy.tsv")
zoraptera_seqs <- read.FASTA("../zoraptera/zoraptera.fasta")

# Get respective ingroup taxa
ingroup <- expanded_taxonomy$TipLabel[expanded_taxonomy$Order=="Zoraptera"]
new_ingroup <- zoraptera_taxonomy$TipLabel[zoraptera_taxonomy$Order=="Zoraptera"]

# Update tree, taxonomy and sequences
expanded_tree <- replace_clade(expanded_tree, ingroup, zoraptera_tree, new_ingroup)
expanded_taxonomy <- update_taxonomy(expanded_taxonomy, ingroup, zoraptera_taxonomy, new_ingroup, cols)
expanded_seqs <- update_sequences(expanded_seqs, ingroup, zoraptera_seqs, new_ingroup)


# Replace Embioptera
# ------------------

cat("Replacing Embioptera\n")

# Read new embiopteran tree (pick the last tree as a sample from the posterior)
trees <- read.nexus("../emb_pha_not/mb_runs/run1/tree_sample.tre")
embioptera_tree <- trees[[length(trees)]]

# Read in new embiopteran taxonomy and sequences
embioptera_taxonomy <- read.delim("../emb_pha_not/emb_pha_not_taxonomy.tsv")
embioptera_seqs <- read.FASTA("../emb_pha_not/emb_pha_not.fasta")

# Get respective ingroup taxa
ingroup <- expanded_taxonomy$TipLabel[expanded_taxonomy$Order=="Embioptera"]
new_ingroup <- embioptera_taxonomy$TipLabel[embioptera_taxonomy$Order=="Embioptera"]

# Update tree, taxonomy and sequences
expanded_tree <- replace_clade(expanded_tree, ingroup, embioptera_tree, new_ingroup)
expanded_taxonomy <- update_taxonomy(expanded_taxonomy, ingroup, embioptera_taxonomy, new_ingroup, cols)
expanded_seqs <- update_sequences(expanded_seqs, ingroup, embioptera_seqs, new_ingroup)


# Replace Grylloblattodea
# -----------------------

cat("Replacing Grylloblattodea\n")

# Read new grylloblattodean tree (pick the last tree as a sample from the posterior)
trees <- read.nexus("../emb_pha_not/mb_runs/run1/tree_sample.tre")
grylloblattodea_tree <- trees[[length(trees)]]

# Read in new grylloblattodean taxonomy and sequences
grylloblattodea_taxonomy <- read.delim("../emb_pha_not/emb_pha_not_taxonomy.tsv")
grylloblattodea_seqs <- read.FASTA("../emb_pha_not/emb_pha_not.fasta")

# Get respective ingroup taxa
ingroup <- expanded_taxonomy$TipLabel[expanded_taxonomy$Order=="Grylloblattodea"]
new_ingroup <- grylloblattodea_taxonomy$TipLabel[grylloblattodea_taxonomy$Order=="Grylloblattodea"]

# Update tree, taxonomy and sequences
expanded_tree <- replace_clade(expanded_tree, ingroup, grylloblattodea_tree, new_ingroup)
expanded_taxonomy <- update_taxonomy(expanded_taxonomy, ingroup, grylloblattodea_taxonomy, new_ingroup, cols)
expanded_seqs <- update_sequences(expanded_seqs, ingroup, grylloblattodea_seqs, new_ingroup)


# Replace Mantophasmatodea
# ------------------------

cat("Replacing Mantophasmatodea\n")

# Read new mantophasmatodean tree (pick the last tree as a sample from the posterior)
trees <- read.nexus("../emb_pha_not/mb_runs/run1/tree_sample.tre")
mantophasmatodea_tree <- trees[[length(trees)]]

# Read in new mantophasmatodean taxonomy and sequences
mantophasmatodea_taxonomy <- read.delim("../emb_pha_not/emb_pha_not_taxonomy.tsv")
mantophasmatodea_seqs <- read.FASTA("../emb_pha_not/emb_pha_not.fasta")

# Get respective ingroup taxa
ingroup <- expanded_taxonomy$TipLabel[expanded_taxonomy$Order=="Mantophasmatodea"]
new_ingroup <- mantophasmatodea_taxonomy$TipLabel[mantophasmatodea_taxonomy$Order=="Mantophasmatodea"]

# Update tree, taxonomy and sequences
expanded_tree <- replace_clade(expanded_tree, ingroup, mantophasmatodea_tree, new_ingroup)
expanded_taxonomy <- update_taxonomy(expanded_taxonomy, ingroup, mantophasmatodea_taxonomy, new_ingroup, cols)
expanded_seqs <- update_sequences(expanded_seqs, ingroup, mantophasmatodea_seqs, new_ingroup)


# Replace Chalcidoidea
# --------------------

cat("Replacing Chalcidoidea\n")

# Set version paths
ver <- 1
tree_path <- paste0("../chalcidoidea/mb_runs/run",ver,"/tree_sample.tre")
tax_path  <- paste0("../chalcidoidea/expanded_chalcidoidea",ver,"_taxonomy.tsv")
seq_path  <- paste0("../chalcidoidea/expanded_chalcidoidea",ver,".fasta")

# Read new chalcid tree (pick the last tree as a sample from the posterior)
trees <- read.nexus(tree_path)
chalcidoidea_tree <- trees[[length(trees)]]

# Read in new chalcidoidea taxonomy and sequences
chalcidoidea_taxonomy <- read.delim(tax_path)
chalcidoidea_seqs <- read.FASTA(seq_path)

# Get respective ingroup taxa
X <- read.delim("../chalcidoidea/chesters_chalcidoidea_taxonomy.tsv")
chesters_chalcidoidea <- unique(X$Family[X$Family!="Mymarommatidae"])
ingroup <- expanded_taxonomy$TipLabel[expanded_taxonomy$Family %in% chesters_chalcidoidea]
new_ingroup <- chalcidoidea_taxonomy$TipLabel[chalcidoidea_taxonomy$Family!="Mymarommatidae"]

# Update tree, taxonomy and sequences
expanded_tree <- replace_clade(expanded_tree, ingroup, chalcidoidea_tree, new_ingroup)
expanded_taxonomy <- update_taxonomy(expanded_taxonomy, ingroup, chalcidoidea_taxonomy, new_ingroup, cols)
expanded_seqs <- update_sequences(expanded_seqs, ingroup, chalcidoidea_seqs, new_ingroup)


# Add in outgroup
# ---------------

cat("Adding in outgroup\n")

# Read entognathan and outgroup tree
trees <- read.nexus("../entognatha_outgroup/mb_runs/run1/tree_sample.tre")
entognatha_outgroup_tree <- trees[[length(trees)]]
entognatha_outgroup_taxonomy <- read.delim("../entognatha_outgroup/entognatha_outgroup_taxonomy.tsv")

# Identify key nodes and heights in root tree
insecta <- entognatha_outgroup_taxonomy$TipLabel[entognatha_outgroup_taxonomy$Class=="Insecta"]
diplura <- entognatha_outgroup_taxonomy$TipLabel[entognatha_outgroup_taxonomy$Class=="Diplura"]
insect_node <- getMRCA(entognatha_outgroup_tree, insecta)
diplura_node <- getMRCA(entognatha_outgroup_tree, diplura)
insert_node <- getMRCA(entognatha_outgroup_tree, c(insecta, diplura))
depths <- node.depth.edgelength(entognatha_outgroup_tree)
entognatha_outgroup_tree_height <- depths[1]     # Total height or depth of tree
insect_height <- depths[insect_node]     # Distance from root to insect node

# Get donor tree and height
donor_tree <- expanded_tree
donor_tree_height <- node.depth.edgelength(donor_tree)[1]

# Rescale the root tree
scale_factor <- donor_tree_height / (entognatha_outgroup_tree_height - insect_height)
entognatha_outgroup_tree$edge.length <- entognatha_outgroup_tree$edge.length * scale_factor
depths <- node.depth.edgelength(entognatha_outgroup_tree)

# Add root edge to donor tree
donor_tree$root.edge <- depths[insect_node] - depths[insert_node]

# Now we are ready to combine the trees
insert_pos <- depths[diplura_node] - depths[insert_node]
receptor_tree <- drop.tip(entognatha_outgroup_tree, insecta)
insert_node <- getMRCA(receptor_tree, diplura)  # Node number of diplura node might have changed
expanded_tree <- bind.tree(receptor_tree, donor_tree, where=insert_node, position=insert_pos)

# Finally add the taxonomy info (and info on translation table)
cols <- c(cols,"TranslationTable")
expanded_taxonomy$TranslationTable <- 5
expanded_taxonomy <- rbind(expanded_taxonomy, entognatha_outgroup_taxonomy[entognatha_outgroup_taxonomy$Class!="Insecta",cols])

# Read root sequences
entognatha_outgroup_seqs <- read.FASTA("../entognatha_outgroup/entognatha_outgroup.fasta")
entognatha_outgroup_seqs <- entognatha_outgroup_seqs[!names(entognatha_outgroup_seqs) %in% insecta]


# Write resulting tree, taxonomy info, and sequences
# --------------------------------------------------

write.tree(expanded_tree,"chesters_expanded.nwk")
write.table(expanded_taxonomy, "chesters_expanded_taxonomy.tsv", sep="\t", row.names=FALSE)
write.FASTA(expanded_seqs, "chesters_expanded.fasta")
write.FASTA(entognatha_outgroup_seqs, "chesters_expanded.fasta", append=TRUE)

