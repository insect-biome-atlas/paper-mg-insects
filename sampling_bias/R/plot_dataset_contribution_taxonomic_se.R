# Script for generating box plots of dataset contribution
# to taxonomic sampling

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
    scale_y_continuous(name = ylab, limits=c(0.0,0.65)) + 
    scale_fill_manual(values = c("deepskyblue","goldenrod4","blue"),
                      labels = c("Malaise", "Litter", "Combined"),
                      name   = "Dataset") +
    ylab(ylab) +
    ggtitle(plot_title) +
    theme(
        legend.position = "bottom",
        legend.key.size = unit(1.5,"cm")
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
meta_se <- read.delim(paste0(path,"CO1_sequencing_metadata_SE.tsv"))
malaise_samples_se <- read.delim(paste0(path,"samples_metadata_malaise_SE.tsv"))
litter_samples_se <- read.delim(paste0(path,"samples_metadata_soil_litter_SE.tsv"))
sites_se <- read.delim(paste0(path,"sites_metadata_SE.tsv"))

# Read in counts data
lysate <- readRDS("../../iba_data/cluster_counts_malaise_long_se.rds")
litter <- readRDS("../../iba_data/cluster_counts_litter_long_se.rds")

# Read in cluster taxonomy info
cluster_taxonomy_se <- readRDS("../../iba_data/cluster_taxonomy_se.rds")

# Add the required metadata
add_meta <- function(D, meta, samples, sites, taxonomy) {
    D$sampleID_FIELD <- meta$sampleID_FIELD[match(D$sampleID_NGI,meta$sampleID_NGI)]
    idx <- match(D$sampleID_FIELD,samples$sampleID_FIELD)
    D$trapID <- samples$trapID[idx]
    D$Order <- taxonomy$Order[match(D$cluster,taxonomy$cluster)]
    return(D)
}
lysate <- add_meta(lysate, meta_se, malaise_samples_se, sites_se, cluster_taxonomy_se)
litter <- add_meta(litter, meta_se, litter_samples_se, sites_se, cluster_taxonomy_se)

# Group other orders
big_five <- c("Diptera","Hymenoptera","Coleoptera","Lepidoptera","Hemiptera")
lysate$Order[!(lysate$Order %in% big_five)] <- "Other"
litter$Order[!(litter$Order %in% big_five)] <- "Other"

# Filter traps based on minimal data requirement
filter_dataset <- function(D, min_count) {
    X <- unique(lysate[,c("trapID","sampleID_NGI")])
    trap_samples <- table(X$trapID)
    traps_to_keep <- names(trap_samples[trap_samples >= min_count])
    D[D$trapID %in% traps_to_keep,]
}
lysate <- filter_dataset(lysate,10)     # Require at least 10 Malaise trap samples
litter <- filter_dataset(litter,10)


# Generate box plot data
# ----------------------

# Aggregate reads for each trap and order
lysate_sites <- aggregate(read_count~trapID+cluster+Order,data=lysate,FUN=sum)
litter_sites <- aggregate(read_count~trapID+cluster+Order,data=litter,FUN=sum)
lysate_litter_sites <- aggregate(read_count~trapID+cluster+Order,data=rbind(lysate,litter),FUN=sum)

compute_dataset_prop <- function(df, dataset_name, Order) {
    X <- df[df$Order==Order,]
    Y <- table(df$trapID)
    sum_total <- nrow(X)
    res <- data.frame(table(X$trapID))
    colnames(res) <- c("trapID","OTUs")
    res$otu_prop <- res$OTUs / as.numeric(Y[match(res$trapID,names(Y))])
    res$dataset <- dataset_name
    return (res)
}

generate_plot_data <- function(Order) {
    res <- data.frame()
    res <- compute_dataset_prop(lysate_sites, "Lysate", Order)
    res <- rbind(res, compute_dataset_prop(litter_sites, "Litter", Order))
    res <- rbind(res, compute_dataset_prop(lysate_litter_sites, "Lysate + litter", Order))
    res$dataset <- factor(res$dataset, levels=c("Lysate","Litter","Lysate + litter"))
    return(res)
}


# Get accumulated values
# ----------------------

get_accumulated_vals <- function(D, dataset_name) {
    X <- unique(D[,c("cluster","Order")])
    df <- data.frame(table(X$Order))
    colnames(df) <- c("Order","OTUs")
    df$otu_prop <- df$OTUs / sum(df$OTUs)
    df$dataset <- dataset_name
    return(df)
}
acc_vals <- get_accumulated_vals(lysate_sites, "Lysate")
acc_vals <- rbind(acc_vals, get_accumulated_vals(litter_sites, "Litter"))
acc_vals <- rbind(acc_vals, get_accumulated_vals(lysate_litter_sites, "Lysate + litter"))


# Make plots
# ------------

D1 <- generate_plot_data("Diptera")
D2 <- generate_plot_data("Hymenoptera")
D3 <- generate_plot_data("Coleoptera")
D4 <- generate_plot_data("Lepidoptera")
D5 <- generate_plot_data("Hemiptera")
D6 <- generate_plot_data("Other")
p1 <- box_plot(D1, acc_vals[acc_vals$Order=="Diptera",], "Diptera", "Proportion of OTUs")
p2 <- box_plot(D2, acc_vals[acc_vals$Order=="Hymenoptera",], "Hymenoptera")
p3 <- box_plot(D3, acc_vals[acc_vals$Order=="Coleoptera",], "Coleptera")
p4 <- box_plot(D4, acc_vals[acc_vals$Order=="Lepidoptera",], "Lepidoptera", "Proportion of OTUs")
p5 <- box_plot(D5, acc_vals[acc_vals$Order=="Hemiptera",], "Hemiptera")
p6 <- box_plot(D6, acc_vals[acc_vals$Order=="Other",], "Other")

# Save plots
ggsave(file = "../figs/Fig_dataset_contribution_taxonomic_se.jpg",
       width = 21,
       height = 14,
       plot = p1 + p2 + p3 + p4 + p5 + p6 +
           plot_layout(guides="collect", ncol=3) +
           plot_annotation(tag_levels="A") &
           theme(legend.position="bottom")
       )

