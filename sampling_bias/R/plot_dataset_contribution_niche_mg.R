# Script for generating box plots of dataset contribution to
# niche sampling of the MG fauna

library(ggplot2)
library(patchwork)
source("../../fig_settings/fig_colours.R")


# Define plot function
# --------------------

box_plot <- function(D, acc_vals, plot_title, ylims, ylab=NULL) {
  ggplot(D, aes(x=dataset, y=otu_prop, fill=dataset)) +
    theme_minimal(base_size=20) +
    geom_boxplot() +
    geom_point(data=acc_vals, aes(x=dataset, y=otu_prop), shape=23, size=7, fill=diamond_col) +
    scale_x_discrete(name = NULL, labels = NULL) +  # drop x labels
    scale_fill_manual(values = c("deepskyblue","goldenrod4","blue"),
                      labels = c("Malaise", "Litter", "Combined"),
                      name   = "Dataset") +
    scale_y_continuous(name=ylab, limits=ylims) +
    ggtitle(plot_title) +
    theme(
        legend.position = "bottom",
#        legend.title = element_text(size = 20),
#        legend.text = element_text(size = 18),
        legend.key.size=unit(1.5,"cm")
#        plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
#        axis.text.x = element_blank(),
#        axis.ticks.x = element_blank()
    )
}

box_plot_ratio <- function(D, acc_vals, plot_title, ylims, ylab=NULL) {
    D$otu_prop <- D$ph_ratio
    acc_vals$otu_prop <- acc_vals$ph_ratio
    box_plot(D, acc_vals, plot_title, ylims, ylab)
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

# Focus on the major niches
niches <- c("Saprophage","Saprophage-parasitoid","Phytophage","Phytophage-parasitoid","Predator","Predator-parasitoid")
lysate <- lysate[lysate$Niche %in% niches,]
litter <- litter[litter$Niche %in% niches,]


# Generate box plot data
# ----------------------

# Aggregate reads for each trap and niche
lysate_sites <- aggregate(read_count~trapID+cluster+Niche,data=lysate,FUN=sum)
litter_sites <- aggregate(read_count~trapID+cluster+Niche,data=litter,FUN=sum)
lysate_litter_sites <- aggregate(read_count~trapID+cluster+Niche,data=rbind(lysate,litter),FUN=sum)

compute_dataset_major_prop <- function(df, dataset_name, Niche) {
    major <- c("Saprophage","Phytophage","Predator")
    X <- df[df$Niche==Niche,]
    Y <- table(df$trapID[df$Niche %in% major])
    res <- data.frame(table(X$trapID))
    colnames(res) <- c("trapID","OTUs")
    res$otu_prop <- res$OTUs / as.numeric(Y[match(res$trapID,names(Y))])
    res$dataset <- dataset_name
    return (res)
}

compute_dataset_ph_ratio <- function(df, dataset_name, Niche) {
    hosts <- c("Saprophage","Phytophage","Predator")
    parasitoids <- c("Saprophage-parasitoid","Phytophage-parasitoid","Predator-parasitoid")
    host <- hosts[which(parasitoids==Niche)]
    X <- df[df$Niche==Niche,]
    Y <- table(df$trapID[df$Niche==host])
    res <- data.frame(table(X$trapID))
    colnames(res) <- c("trapID","OTUs")
    res$ph_ratio <- res$OTUs / as.numeric(Y[match(res$trapID,names(Y))])
    res$dataset <- dataset_name
    return (res)
}

generate_major_plot_data <- function(Niche) {
    res <- data.frame()
    res <- compute_dataset_major_prop(lysate_sites, "Lysate", Niche)
    res <- rbind(res, compute_dataset_major_prop(litter_sites, "Litter", Niche))
    res <- rbind(res, compute_dataset_major_prop(lysate_litter_sites, "Lysate + litter", Niche))
    res$dataset <- factor(res$dataset, levels=c("Lysate","Litter","Lysate + litter"))
    return(res)
}

generate_parasitoid_plot_data <- function(Niche) {
    res <- data.frame()
    res <- compute_dataset_ph_ratio(lysate_sites, "Lysate", Niche)
    res <- rbind(res, compute_dataset_ph_ratio(litter_sites, "Litter", Niche))
    res <- rbind(res, compute_dataset_ph_ratio(lysate_litter_sites, "Lysate + litter", Niche))
    res$dataset <- factor(res$dataset, levels=c("Lysate","Litter","Lysate + litter"))
    return(res)
}


# Get accumulated values
# ----------------------

get_accumulated_major_props <- function(D, dataset_name) {
    major <- c("Saprophage","Phytophage","Predator")
    X <- unique(D[,c("cluster","Niche")])
    X <- X[X$Niche %in% major,]
    df <- data.frame(table(X$Niche))
    colnames(df) <- c("Niche","OTUs")
    df$otu_prop <- df$OTUs / sum(df$OTUs)
    df$dataset <- dataset_name
    return(df)
}

get_accumulated_ph_ratios <- function(D, dataset_name) {
    hosts <- c("Saprophage","Phytophage","Predator")
    parasitoids <- c("Saprophage-parasitoid","Phytophage-parasitoid","Predator-parasitoid")
    X <- unique(D[,c("cluster","Niche")])
    df <- data.frame(table(X$Niche))
    colnames(df) <- c("Niche","OTUs")
    res <- data.frame()
    for (i in 1:length(hosts)) {
        ph_ratio <- df$OTUs[df$Niche==parasitoids[i]] / df$OTUs[df$Niche==hosts[i]]
        res <- rbind(res, data.frame(list(Niche=parasitoids[i],ph_ratio=ph_ratio)))
    }
    res$dataset <- dataset_name
    return(res)
}

acc_props <- get_accumulated_major_props(lysate_sites, "Lysate")
acc_props <- rbind(acc_props, get_accumulated_major_props(litter_sites, "Litter"))
acc_props <- rbind(acc_props, get_accumulated_major_props(lysate_litter_sites, "Lysate + litter"))

acc_ratios <- get_accumulated_ph_ratios(lysate_sites, "Lysate")
acc_ratios <- rbind(acc_ratios, get_accumulated_ph_ratios(litter_sites, "Litter"))
acc_ratios <- rbind(acc_ratios, get_accumulated_ph_ratios(lysate_litter_sites, "Lysate + litter"))


# Make plots
# ------------

D1 <- generate_major_plot_data("Saprophage")
D2 <- generate_major_plot_data("Phytophage")
D3 <- generate_major_plot_data("Predator")
D4 <- generate_parasitoid_plot_data("Saprophage-parasitoid")
D5 <- generate_parasitoid_plot_data("Phytophage-parasitoid")
D6 <- generate_parasitoid_plot_data("Predator-parasitoid")
p1 <- box_plot(D1, acc_props[acc_props$Niche=="Saprophage",],  "Saprophage", c(0.0,0.8), "Proportion of OTUs")
p2 <- box_plot(D2, acc_props[acc_props$Niche=="Phytophage",], "Phytophage", c(0.0,0.8))
p3 <- box_plot(D3, acc_props[acc_props$Niche=="Predator",], "Predator", c(0.0,0.8))
p4 <- box_plot_ratio(D4, acc_ratios[acc_ratios$Niche=="Saprophage-parasitoid",], "Saprophage-parasitoid", c(0.0,1.05), "Parasitoid:host OTU ratio")
p5 <- box_plot_ratio(D5, acc_ratios[acc_ratios$Niche=="Phytophage-parasitoid",], "Phytophage-parasitoid", c(0.0,1.05))
p6 <- box_plot_ratio(D6, acc_ratios[acc_ratios$Niche=="Predator-parasitoid",], "Predator-parasitoid", c(0.0,1.05))

# Save plots
ggsave(file = "../figs/Fig_dataset_contribution_niche_mg.jpg",
       width = 21,
       height = 14,
       plot = p1 + p2 + p3 + p4 + p5 + p6 +
           plot_layout(guides="collect", ncol=3) +
           plot_annotation(tag_levels="A") &
           theme(legend.position="bottom")
       )

