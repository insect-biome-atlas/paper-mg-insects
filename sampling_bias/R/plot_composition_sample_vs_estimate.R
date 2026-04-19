# Script generating bar charts showing taxonomic and ecological composition
# of the sampled and estimated SE faunas

library(ggplot2)
library(patchwork)

# Define niches and habitats of interest
# (negligible fraction of species outside)
# --------------------------------------

niches <- c("Saprophage","Saprophage-parasitoid","Phytophage","Phytophage-parasitoid","Predator","Predator-parasitoid")
habitats <- c("Plants","Soil","Water","Wood","Temporary habitats","Fungi")


# Get expected values from 2020 estimate
# --------------------------------------

D1 <- read.delim("../data/niche_comp_by_order.tsv")
D2 <- read.delim("../data/habitat_comp_by_order.tsv")
D1 <- D1[D1$niche %in% niches,]
D2 <- D2[D2$habitat %in% habitats,]
D1$study <- "2020 estimate"
D2$study <- "2020 estimate"
D1$otu_prop <- numeric(nrow(D1))
D2$otu_prop <- numeric(nrow(D2))
sum_total1 <- aggregate(species~order, data=D1, FUN=sum)
sum_total2 <- aggregate(species~order, data=D2, FUN=sum)
D1$otu_prop <- D1$species / sum_total1$species[match(D1$order,sum_total1$order)]
D2$otu_prop <- D2$species / sum_total2$species[match(D2$order,sum_total2$order)]


# Get values from current study
# -----------------------------

# Read in taxonomy data
T <- read.delim("../../iba_data/cluster_taxonomy_se.tsv")

# Read in life history trait data
lht <- read.delim("../../traits/clade_trait_data_se.tsv") 
remove_annot <- function(D) {
    ranks <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
    idx <- which(colnames(D) %in% ranks)
    return (D[,-idx])
}  
lht <- remove_annot(lht)

# Merge data frames
T <- merge(T,lht,by="Clade")

# Get counts data
malaise <- readRDS("../../iba_data/cluster_counts_malaise_long_se.rds")
litter <- readRDS("../../iba_data/cluster_counts_litter_long_se.rds")
all_clusters <- unique(c(malaise$cluster,litter$cluster))

# Restrict taxonomy to encountered clusters. Add order group
T <- T[T$cluster %in% all_clusters,]
big_five <- c("Diptera","Hymenoptera","Coleoptera","Lepidoptera","Hemiptera")
T$Order[!(T$Order %in% big_five)] <- "Other"
orders <- c(big_five,"Other")


# Compute current study composition by order
# ------------------------------------------

compute_niche_composition <- function(df, orders) {
    X <- data.frame(table(df$Niche,by=df$Order))
    colnames(X) <- c("niche","order","species")
    Y <- aggregate(species~order, data=X, FUN=sum)
    X$otu_prop <- X$species / Y$species[match(X$order,Y$order)]
    return (X)
}

compute_habitat_composition <- function(df, orders) {
    X <- data.frame(table(df$Habitat,by=df$Order))
    colnames(X) <- c("habitat","order","species")
    Y <- aggregate(species~order, data=X, FUN=sum)
    X$otu_prop <- X$species / Y$species[match(X$order,Y$order)]
    return (X)
}

T1 <- T[T$Niche %in% niches,]
E1 <- compute_niche_composition(T1, orders)

T2 <- T[T$Habitat %in% habitats,]
E2 <- compute_habitat_composition(T2, orders)

E1$study <- "Current study"
E2$study <- "Current study"


# Merge data series
# -----------------

niche_comp <- rbind(D1,E1)
habitat_comp <- rbind(D2,E2)
niche_comp$study <- factor(niche_comp$study, levels=c("Current study","2020 estimate"))
habitat_comp$study <- factor(habitat_comp$study, levels=c("Current study","2020 estimate"))


# Generate plots
# --------------

# Plot functions
plot_Niche <- function(df, plot_title) {
    ggplot(data=df, aes(x=niche, y=otu_prop, group=study, fill=study)) +
        geom_col(position="dodge") +
        coord_flip() +
        scale_fill_manual(values = c("blue","linen")) +
        labs(title = plot_title,
             x = NULL,
             y = "Proportion of OTUs",
             fill = "Study") +
        theme_minimal(base_size=17)
}

plot_Habitat <- function(df, plot_title) {
    ggplot(data=df, aes(x=habitat, y=otu_prop, group=study, fill=study)) +
        geom_col(position="dodge") +
        coord_flip() +
        scale_fill_manual(values = c("blue","linen")) +
        labs(title = plot_title,
             x = NULL,
             y = "Proportion of OTUs",
             fill = "Study") +
        theme_minimal(base_size=17)
}

# Render plots
p1 <- plot_Niche(niche_comp[niche_comp$order=="Coleoptera",], "Coleoptera niches")
p2 <- plot_Habitat(habitat_comp[habitat_comp$order=="Coleoptera",], "Coleoptera microhabitats")
p3 <- plot_Niche(niche_comp[niche_comp$order=="Diptera",], "Diptera niches")
p4 <- plot_Habitat(habitat_comp[habitat_comp$order=="Diptera",], "Diptera microhabitats")
p5 <- plot_Niche(niche_comp[niche_comp$order=="Hymenoptera",], "Hymenoptera niches")
p6 <- plot_Habitat(habitat_comp[habitat_comp$order=="Hymenoptera",], "Hymenoptera microhabitats")

# Put plots together
ggsave(
       file = "../figs/Fig_composition_sample_vs_estimate.jpg",
       width = 14.0,
       height = 18.0,
       plot = p1 + p2 + p3 + p4 + p5 + p6 +
              plot_layout(ncol=2, axis_titles="collect", guides="collect") +
              plot_annotation(tag_levels="A") & theme(legend.position="bottom")
       )

