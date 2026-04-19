# Script for generating box plots of sampling method bias
# in habitat sampling

library(ggplot2)
library(patchwork)


# Define plot function
# --------------------

box_plot <- function(D, estimated_val, acc_vals, plot_title) {
  ggplot(D, aes(x=dataset, y=otu_prop, fill=dataset)) +
    theme_minimal(base_size=20) +
    geom_boxplot() +
    geom_point(data=acc_vals, aes(x=dataset, y=otu_prop), shape=23, size=3, fill="red") +
    geom_hline(yintercept=estimated_val, colour="red", linetype="dotted", linewidth=1) +
    scale_x_discrete(name = NULL, labels = NULL) +  # drop x labels
    scale_fill_manual(values = c("deepskyblue","darkolivegreen","linen","aquamarine","blue","brown"),
#                      labels = c("Lysate", "Homogenate", "Litter", "Lysate + homogenate", "Lysate + litter", "All"),
                      name   = "Dataset") +
    ylab("Proportion of OTUs") +
    ggtitle(plot_title) +
    theme(
        legend.position = "bottom",
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
homogenate <- readRDS("../../iba_data/cluster_counts_homogenate_long_se.rds")
litter <- readRDS("../../iba_data/cluster_counts_litter_long_se.rds")

# Read in cluster taxonomy info
cluster_taxonomy_se <- read.delim("../../iba_data/cluster_taxonomy_se.tsv")

# Read in traits data
lht_se <- read.delim("../../traits/clade_trait_data_se.tsv")

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
lysate <- add_meta(lysate, meta_se, malaise_samples_se, sites_se, cluster_taxonomy_se, lht_se)
homogenate <- add_meta(homogenate, meta_se, malaise_samples_se, sites_se, cluster_taxonomy_se, lht_se)
litter <- add_meta(litter, meta_se, litter_samples_se, sites_se, cluster_taxonomy_se, lht_se)

# Focus on major niches
habitats <- c("Plants","Soil","Water","Wood","Temporary habitats","Fungi")
lysate <- lysate[lysate$Habitat %in% habitats,]
homogenate <- homogenate[homogenate$Habitat %in% habitats,]
litter <- litter[litter$Habitat %in% habitats,]

# Filter traps based on minimal data requirement for homogenates
filter_dataset <- function(D, min_homogenate_count) {
    X <- unique(homogenate[,c("trapID","sampleID_NGI")])
    trap_samples <- table(X$trapID)
    traps_to_keep <- names(trap_samples[trap_samples >= min_homogenate_count])
    D[D$trapID %in% traps_to_keep,]
}
lysate <- filter_dataset(lysate,3)          # 3 homogenate samples = 12 lysate samples
homogenate <- filter_dataset(homogenate,3)
litter <- filter_dataset(litter,3)


# Get expected values from 2020 estimate
# --------------------------------------

habitat_comp <- read.delim("../data/habitat_comp.tsv")
sum_total <- sum(habitat_comp$species)
prop_estimate <- function(habitat) { habitat_comp$species[habitat_comp$habitat==habitat] / sum_total }
plants_estimate <- prop_estimate("Plants")
soil_estimate <- prop_estimate("Soil")
water_estimate <- prop_estimate("Water")
wood_estimate <- prop_estimate("Wood")
temporary_estimate <- prop_estimate("Temporary habitats")
fungi_estimate <- prop_estimate("Fungi")


# Generate box plot data
# ----------------------

# Aggregate reads for each trap and order
lysate_sites <- aggregate(read_count~trapID+cluster+Habitat,data=lysate,FUN=sum)
homogenate_sites <- aggregate(read_count~trapID+cluster+Habitat,data=homogenate,FUN=sum)
litter_sites <- aggregate(read_count~trapID+cluster+Habitat,data=litter,FUN=sum)
lysate_homogenate_sites <- aggregate(read_count~trapID+cluster+Habitat,data=rbind(lysate,homogenate),FUN=sum)
lysate_litter_sites <- aggregate(read_count~trapID+cluster+Habitat,data=rbind(lysate,homogenate),FUN=sum)
all_sites <- aggregate(read_count~trapID+cluster+Habitat,data=rbind(lysate,homogenate,litter),FUN=sum)

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
    res <- rbind(res, compute_dataset_prop(homogenate_sites, "Homogenate", Habitat))
    res <- rbind(res, compute_dataset_prop(litter_sites, "Litter", Habitat))
    res <- rbind(res, compute_dataset_prop(lysate_homogenate_sites, "Lysate + homogenate", Habitat))
    res <- rbind(res, compute_dataset_prop(lysate_litter_sites, "Lysate + litter", Habitat))
    res <- rbind(res, compute_dataset_prop(all_sites, "All", Habitat))
    res$dataset <- factor(res$dataset, levels=c("Lysate","Homogenate","Litter","Lysate + homogenate","Lysate + litter", "All"))
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
acc_vals <- rbind(acc_vals, get_accumulated_vals(homogenate_sites, "Homogenate"))
acc_vals <- rbind(acc_vals, get_accumulated_vals(litter_sites, "Litter"))
acc_vals <- rbind(acc_vals, get_accumulated_vals(lysate_homogenate_sites, "Lysate + homogenate"))
acc_vals <- rbind(acc_vals, get_accumulated_vals(lysate_litter_sites, "Lysate + litter"))
acc_vals <- rbind(acc_vals, get_accumulated_vals(all_sites, "All"))


# Make plots
# ------------

D1 <- generate_plot_data("Plants")
D2 <- generate_plot_data("Soil")
D3 <- generate_plot_data("Water")
D4 <- generate_plot_data("Wood")
D5 <- generate_plot_data("Temporary habitats")
D6 <- generate_plot_data("Fungi")
p1 <- box_plot(D1, dipt_estimate, acc_vals[acc_vals$Habitat=="Plants",], "Plants")
p2 <- box_plot(D2, hyme_estimate, acc_vals[acc_vals$Habitat=="Soil",], "Soil")
p3 <- box_plot(D3, cole_estimate, acc_vals[acc_vals$Habitat=="Water",], "Water")
p4 <- box_plot(D4, lepi_estimate, acc_vals[acc_vals$Habitat=="Wood",], "Wood")
p5 <- box_plot(D5, hemi_estimate, acc_vals[acc_vals$Habitat=="Temporary habitats",], "Temporary habitats")
p6 <- box_plot(D6, other_estimate, acc_vals[acc_vals$Habitat=="Fungi",], "Fungi")

# Save plots
ggsave(file = "../figs/Fig_sampling_bias_habitat.jpg",
       width = 14,
       height = 21,
       plot = p1 + p2 + p3 + p4 + p5 + p6 +
           plot_layout(axis_titles="collect_y", guides="collect", ncol=2) +
           plot_annotation(tag_levels="A") &
           theme(legend.position="bottom")
       )

