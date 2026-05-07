# Plot site composition data

library(ggplot2)
library(patchwork)


# Read in data
# ------------

otu_site_meta_mg <- read.delim("../data/otu_site_meta_mg.tsv")
otu_site_meta_se <- read.delim("../data/otu_site_meta_se.tsv")
source("filter_se_traps.R")
traps_to_keep <- filter_se_traps(20)    # Only keep se traps with 20 or more samples
otu_site_meta_se <- otu_site_meta_se[otu_site_meta_se$trapID %in% traps_to_keep,]
otu_site_meta_rf <- otu_site_meta_mg[otu_site_meta_mg$trap_habitat %in% c("Rainforest","Montane_Rainforest"),]
otu_site_meta_df <- otu_site_meta_mg[otu_site_meta_mg$trap_habitat == "Dry_Forest",]
otu_site_meta_rf$trap_habitat <- "Rainforest"
otu_site_meta_df$trap_habitat <- "Dry forest"
otu_site_meta_se$trap_habitat <- "Temperate forest"
D <- rbind(otu_site_meta_rf,otu_site_meta_df,otu_site_meta_se)
D$trap_habitat <- factor(D$trap_habitat,levels=c("Rainforest","Dry forest","Temperate forest"))
idx <- which(D$Habitat=="Temporary habitats")
D$Habitat[idx] <- "Temporary"

# Define plot functions
# ---------------------

box_plot <- function(X, col, plot_title=NULL, y_label=NULL, y_limits=numeric()) {
    idx <- which(colnames(X)==col)
    colnames(X)[idx] <- "count"
    if (length(y_limits)==0)
        y_limits <- c(min(X$count),max(X$count))
    ggplot(X, aes(x=trap_habitat, y=count, fill=trap_habitat)) +
        theme_minimal() +
        geom_boxplot() +
        scale_x_discrete(name = NULL, labels = NULL) +  # drop x labels
        scale_fill_manual(values = c("darkgreen","lightgreen","blue"), labels = c("Rainforest", "Tropical dry forest","Temperate forest"), name = "Habitat") +
        ylab(y_label) +
        ylim(y_limits) +
        ggtitle(plot_title) +
        theme(
            legend.position = "bottom",
            legend.title = element_text(size = 10),
            legend.text = element_text(size = 9)
#      plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
#      axis.text.x = element_blank(),
#      axis.ticks.x = element_blank()
    )
}


# Compute and plot niche composition
# ----------------------------------

traps <- unique(D$trapID)
trap_habitat <- D$trap_habitat[match(traps,D$trapID)]
res <- data.frame(list(trapID=traps,trap_habitat=trap_habitat))
hosts <- c("Saprophage","Phytophage","Predator")
parasitoids <- c("Saprophage-parasitoid","Phytophage-parasitoid","Predator-parasitoid")
for (niche in c(hosts,parasitoids)) {
    X <- data.frame(table(D$trapID[D$Niche==niche]))
    colnames(X) <- c("trapID",niche)
    res <- merge(res,X)
}
for (i in 1:length(hosts)) {
    host <- hosts[i]
    parasitoid <- parasitoids[i]
    res[,parasitoid] <- res[,parasitoid] / res[,host]
}
sum_major <- rowSums(res[,hosts])
for (niche in hosts)
    res[,niche] <- res[,niche]/sum_major

y_limits1 <- c(0.0,0.8)
y_limits2 <- c(0.0,1.05)
p1 <- box_plot(res, "Saprophage", "Saprophage", "Proportion of major-niche OTUs", y_limits1)
p2 <- box_plot(res, "Saprophage-parasitoid", "Saprophage parasitoid", "Parasitoid:host OTU ratio", y_limits2)
p3 <- box_plot(res, "Phytophage", "Phytophage", "Proportion of major-niche OTUs", y_limits1)
p4 <- box_plot(res, "Phytophage-parasitoid", "Phytophage parasitoid", "Parasitoid:host OTU ratio", y_limits2)
p5 <- box_plot(res, "Predator", "Predator", "Proportion of major-niche OTUs", y_limits1)
p6 <- box_plot(res, "Predator-parasitoid", "Predator parasitoid", "Parasitoid:host OTU ratio", y_limits2)

ggsave(file = "../figs/Fig_composition_by_forest_type_niche.jpg",
       width = 7,
       height = 10.5,
       plot = p1 + p2 + p3 + p4 + p5 + p6 +
           plot_layout(axis_titles="collect_y", guides="collect", ncol=2) +
           plot_annotation(tag_levels="A") &
           theme(legend.position="bottom")
       )

# Compute and plot habitat composition
# ----------------------------------

res <- data.frame(list(trapID=traps,trap_habitat=trap_habitat))
habitats <- c("Plants","Soil","Water","Wood","Temporary","Fungi")
for (habitat in habitats) {
    X <- data.frame(table(D$trapID[D$Habitat==habitat]))
    colnames(X) <- c("trapID",habitat)
    res <- merge(res,X)
}
sum_habitats <- rowSums(res[,habitats])
for (habitat in habitats)
    res[,habitat] <- res[,habitat]/sum_habitats

y_limits <- c(0.0,0.7)
p1 <- box_plot(res, "Plants", "Plants", "Proportion of OTUs", y_limits)
p2 <- box_plot(res, "Soil", "Soil", "Proportion of OTUs", y_limits)
p3 <- box_plot(res, "Water", "Water", "Proportion of OTUs", y_limits)
p4 <- box_plot(res, "Wood", "Wood", "Proportion of OTUs", y_limits)
p5 <- box_plot(res, "Temporary", "Temporary", "Proportion of OTUs", y_limits)
p6 <- box_plot(res, "Fungi", "Fungi", "Proportion of OTUs", y_limits)

ggsave(file = "../figs/Fig_composition_by_forest_type_habitat.jpg",
       width = 7,
       height = 10.5,
       plot = p1 + p2 + p3 + p4 + p5 + p6 +
           plot_layout(axis_titles="collect_y", guides="collect", ncol=2) +
           plot_annotation(tag_levels="A") &
           theme(legend.position="bottom")
       )

# Compute and plot taxonomic composition
# --------------------------------------

res <- data.frame(list(trapID=traps,trap_habitat=trap_habitat))
big_five <- c("Diptera","Hymenoptera","Coleoptera","Lepidoptera","Hemiptera")
for (order in big_five) {
    X <- data.frame(table(D$trapID[D$Order==order]))
    colnames(X) <- c("trapID",order)
    res <- merge(res,X)
}
X <- data.frame(table(D$trapID[!(D$Order %in% big_five)]))
colnames(X) <- c("trapID","Other")
res <- merge(res,X)
cols <- c(big_five,"Other")
sum_orders <- rowSums(res[,cols])
for (col in cols)
    res[,col] <- res[,col]/sum_orders

y_limits <- c(0.0,0.65)
p1 <- box_plot(res, "Diptera", "Diptera", "Proportion of OTUs", y_limits)
p2 <- box_plot(res, "Hymenoptera", "Hymenoptera", "Proportion of OTUs", y_limits)
p3 <- box_plot(res, "Coleoptera", "Coleoptera", "Proportion of OTUs", y_limits)
p4 <- box_plot(res, "Lepidoptera", "Lepidoptera", "Proportion of OTUs", y_limits)
p5 <- box_plot(res, "Hemiptera", "Hemiptera", "Proportion of OTUs", y_limits)
p6 <- box_plot(res, "Other", "Other", "Proportion of OTUs", y_limits)

ggsave(file = "../figs/Fig_composition_by_forest_type_taxonomic.jpg",
       width = 7,
       height = 10.5,
       plot = p1 + p2 + p3 + p4 + p5 + p6 +
           plot_layout(axis_titles="collect_y", guides="collect", ncol=2) +
           plot_annotation(tag_levels="A") &
           theme(legend.position="bottom")
       )

