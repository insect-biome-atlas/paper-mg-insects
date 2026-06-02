# Script for generating box plots of dataset contribution to
# microhabitat sampling of the MG fauna

library(ggplot2)
library(patchwork)


# Define plot function
# --------------------

box_plot <- function(D, acc_vals, plot_title, ylab=NULL) {
  ggplot(D, aes(x=dataset, y=otu_prop, fill=dataset)) +
    theme_minimal(base_size=20) +
    geom_boxplot() +
    geom_point(data=acc_vals, aes(x=dataset, y=otu_prop), shape=23, size=7, fill="red") +
    scale_x_discrete(name = NULL, labels = NULL) +  # drop x labels
    scale_y_continuous(name=ylab, limits=c(0.0,0.65)) +
    scale_fill_manual(values = c("deepskyblue","goldenrod4","blue"),
                      labels = c("Malaise", "Litter", "Combined"),
                      name   = "Dataset") +
    ylab("Proportion of OTUs") +
    ggtitle(plot_title) +
    theme(
        legend.position = "bottom",
        legend.key.size=unit(1.5,"cm")
#        legend.title = element_text(size = 10),
#        legend.text = element_text(size = 9)
#        plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
#        axis.text.x = element_blank(),
#        axis.ticks.x = element_blank()
    )
}


# Read in and prepare data
# -------------------------

# Read in sample metadata
# TODO: Set path to IBA metadata on your system
path <- "~/dev/figshare-repos/iba/raw_data/v6/"
meta_mg <- read.delim(paste0(path,"CO1_sequencing_metadata_MG.tsv"))
malaise_samples_mg <- read.delim(paste0(path,"samples_metadata_malaise_MG.tsv"))
litter_samples_mg <- read.delim(paste0(path,"samples_metadata_litter_MG.tsv"))
sites_mg <- read.delim(paste0(path,"sites_metadata_MG.tsv"))

# Read in counts data
lysate <- readRDS("../../iba_data/cluster_counts_malaise_long_mg.rds")
litter <- readRDS("../../iba_data/cluster_counts_litter_long_mg.rds")

# Read in cluster taxonomy info
cluster_taxonomy_mg <- readRDS("../../iba_data/cluster_taxonomy_mg.rds")

# Read in traits data
lht_mg <- read.delim("../../traits/clade_trait_data_mg.tsv")

# Add the required metadata
add_meta <- function(D, meta, samples, sites, taxonomy, lht) {
    D$sampleID_FIELD <- meta$sampleID_FIELD[match(D$sampleID_NGI,meta$sampleID_NGI)]
    idx <- match(D$sampleID_FIELD,samples$sampleID_FIELD)
    D$trapID <- samples$trapID[idx]
    D$Order <- taxonomy$Order[match(D$cluster,taxonomy$cluster)]
    D$Clade <- taxonomy$Clade[match(D$cluster,taxonomy$cluster)]
    D$Niche <- lht$Niche[match(D$Clade,lht$Clade)]
    D$Habitat <- lht$Habitat[match(D$Clade,lht$Clade)]
    return(D)
}
lysate <- add_meta(lysate, meta_mg, malaise_samples_mg, sites_mg, cluster_taxonomy_mg, lht_mg)
litter <- add_meta(litter, meta_mg, litter_samples_mg, sites_mg, cluster_taxonomy_mg, lht_mg)


# Focus on the major habitats
habitats <- c("Plants","Soil","Water","Wood","Temporary habitats","Fungi")
lysate <- lysate[lysate$Habitat %in% habitats,]
litter <- litter[litter$Habitat %in% habitats,]


# Generate box plot data
# ----------------------

# Aggregate reads for each trap and order
lysate_sites <- aggregate(read_count~trapID+cluster+Habitat,data=lysate,FUN=sum)
litter_sites <- aggregate(read_count~trapID+cluster+Habitat,data=litter,FUN=sum)
lysate_litter_sites <- aggregate(read_count~trapID+cluster+Habitat,data=rbind(lysate,litter),FUN=sum)

compute_dataset_prop <- function(df, dataset_name, Habitat) {
    X <- df[df$Habitat==Habitat,]
    Y <- table(df$trapID)
    sum_total <- nrow(X)
    res <- data.frame(table(X$trapID))
    colnames(res) <- c("trapID","OTUs")
    res$otu_prop <- res$OTUs / as.numeric(Y[match(res$trapID,names(Y))])
    res$dataset <- dataset_name
    return (res)
}

generate_plot_data <- function(Habitat) {
    res <- data.frame()
    res <- compute_dataset_prop(lysate_sites, "Lysate", Habitat)
    res <- rbind(res, compute_dataset_prop(litter_sites, "Litter", Habitat))
    res <- rbind(res, compute_dataset_prop(lysate_litter_sites, "Lysate + litter", Habitat))
    res$dataset <- factor(res$dataset, levels=c("Lysate","Litter","Lysate + litter"))
    return(res)
}


# Get accumulated values
# ----------------------

get_accumulated_vals <- function(D, dataset_name) {
    X <- unique(D[,c("cluster","Habitat")])
    df <- data.frame(table(X$Habitat))
    colnames(df) <- c("Habitat","OTUs")
    df$otu_prop <- df$OTUs / sum(df$OTUs)
    df$dataset <- dataset_name
    return(df)
}
acc_vals <- get_accumulated_vals(lysate_sites, "Lysate")
acc_vals <- rbind(acc_vals, get_accumulated_vals(litter_sites, "Litter"))
acc_vals <- rbind(acc_vals, get_accumulated_vals(lysate_litter_sites, "Lysate + litter"))


# Make plots
# ------------

D1 <- generate_plot_data("Plants")
D2 <- generate_plot_data("Soil")
D3 <- generate_plot_data("Water")
D4 <- generate_plot_data("Wood")
D5 <- generate_plot_data("Temporary habitats")
D6 <- generate_plot_data("Fungi")
p1 <- box_plot(D1, acc_vals[acc_vals$Habitat=="Plants",], "Plants", "Proportion of OTUs")
p2 <- box_plot(D2, acc_vals[acc_vals$Habitat=="Soil",], "Soil")
p3 <- box_plot(D3, acc_vals[acc_vals$Habitat=="Water",],  "Water")
p4 <- box_plot(D4, acc_vals[acc_vals$Habitat=="Wood",], "Wood", "Proportion of OTUs")
p5 <- box_plot(D5, acc_vals[acc_vals$Habitat=="Temporary habitats",], "Temporary habitats")
p6 <- box_plot(D6, acc_vals[acc_vals$Habitat=="Fungi",], "Fungi")

# Save plots
ggsave(file = "../figs/Fig_dataset_contribution_habitat_mg.jpg",
       width = 21,
       height = 14,
       plot = p1 + p2 + p3 + p4 + p5 + p6 +
           plot_layout(guides="collect", ncol=3) +
           plot_annotation(tag_levels="A") &
           theme(legend.position="bottom")
       )

