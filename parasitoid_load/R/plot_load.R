# Script for generating box plots of parasitoid loads per trap site

library(ggplot2)
library(patchwork)


# Define plot function
box_plot <- function(D, plot_title = NULL) {
  ggplot(D, aes(x=trap_habitat, y=ph_ratio, fill=trap_habitat)) +
    theme_minimal() +
    geom_boxplot() +
    scale_x_discrete(name = NULL, labels = NULL) +  # drop x labels
    scale_fill_manual(values = c("darkgreen","lightgreen","blue"), labels = c("Rainforest", "Tropical dry forest","Temperate forest"), name = "Habitat") +
    ylab("Parasite:host read ratio") +
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


# Read in and aggregate data
D1 <- readRDS("../../iba_data/cluster_counts_malaise_long_mg.rds")
D2 <- readRDS("../../iba_data/cluster_counts_malaise_long_se.rds")

meta1 <- read.delim("../../iba_data/malaise_sample_meta_mg.tsv")
meta2 <- read.delim("../../iba_data/malaise_sample_meta_se.tsv")

D1$trapID <- meta1$trapID[match(D1$sampleID_NGI,meta1$sampleID_NGI)]
D2$trapID <- meta2$trapID[match(D2$sampleID_NGI,meta2$sampleID_NGI)]

X <- unique(D2[,c("trapID","sampleID_NGI")])
trap_samples <- table(X$trapID)
traps_to_keep <- names(trap_samples[trap_samples >= 20])    # drop Swedish traps with few samples
D2 <- D2[D2$trapID %in% traps_to_keep,]

D1$trap_habitat <- meta1$trap_habitat[match(D1$sampleID_NGI,meta1$sampleID_NGI)]
D2$trap_habitat <- meta2$trap_habitat[match(D2$sampleID_NGI,meta2$sampleID_NGI)]
D2 <- D2[D2$trap_habitat=="forest",]
D2$trap_habitat <- "Temperate forest"
D1$trap_habitat[D1$trap_habitat=="Dry_Forest"] <- "Tropical dry forest"
D1$trap_habitat[D1$trap_habitat=="Montane_Rainforest"] <- "Rainforest"

T1 <- readRDS("../../iba_data/cluster_taxonomy_mg.rds")
T2 <- readRDS("../../iba_data/cluster_taxonomy_se.rds")

D1$Clade <- T1$Clade[match(D1$cluster,T1$cluster)]
D2$Clade <- T2$Clade[match(D2$cluster,T2$cluster)]

lht1 <- read.delim("../../traits/clade_trait_data_mg.tsv")
lht2 <- read.delim("../../traits/clade_trait_data_se.tsv")

D1$Niche <- lht1$Niche[match(D1$Clade,lht1$Clade)]
D2$Niche <- lht2$Niche[match(D2$Clade,lht2$Clade)]
D1$Habitat <- lht1$Habitat[match(D1$Clade,lht1$Clade)]
D2$Habitat <- lht2$Habitat[match(D2$Clade,lht2$Clade)]

# Aggregate reads for each trap and niche and trap habitat
X1 <- aggregate(read_count~trapID+Niche+trap_habitat,data=D1,FUN=sum)
X2 <- aggregate(read_count~trapID+Niche+trap_habitat,data=D2,FUN=sum)
X <- rbind(X1,X2)

# Generate box plot data for saprophage communities
B1 <- X[X$Niche=="Saprophage",]
B2 <- X[X$Niche=="Saprophage-parasitoid",]
sapro <- data.frame(trapID=B1$trapID)
sapro$trap_habitat <- B1$trap_habitat[match(sapro$trapID,B1$trapID)]
sapro$trap_habitat <- factor(sapro$trap_habitat, levels=c("Rainforest","Tropical dry forest","Temperate forest"))
sapro$h_reads <- B1$read_count[match(sapro$trapID,B1$trapID)]
sapro$p_reads <- B2$read_count[match(sapro$trapID,B2$trapID)]
sapro$p_reads[is.na(sapro$p_reads)] <- 0
sapro$ph_ratio <- sapro$p_reads/sapro$h_reads

# Generate box plot data for phytophage communities
B1 <- X[X$Niche=="Phytophage",]
B2 <- X[X$Niche=="Phytophage-parasitoid",]
phyto <- data.frame(trapID=B1$trapID)
phyto$trap_habitat <- B1$trap_habitat[match(phyto$trapID,B1$trapID)]
phyto$trap_habitat <- factor(phyto$trap_habitat, levels=c("Rainforest","Tropical dry forest","Temperate forest"))
phyto$h_reads <- B1$read_count[match(phyto$trapID,B1$trapID)]
phyto$p_reads <- B2$read_count[match(phyto$trapID,B2$trapID)]
phyto$p_reads[is.na(phyto$p_reads)] <- 0
phyto$ph_ratio <- phyto$p_reads/phyto$h_reads

# Generate box plot data for predator communities
B1 <- X[X$Niche=="Predator",]
B2 <- X[X$Niche=="Predator-parasitoid",]
pred <- data.frame(trapID=B1$trapID)
pred$trap_habitat <- B1$trap_habitat[match(pred$trapID,B1$trapID)]
pred$trap_habitat <- factor(pred$trap_habitat, levels=c("Rainforest","Tropical dry forest","Temperate forest"))
pred$h_reads <- B1$read_count[match(pred$trapID,B1$trapID)]
pred$p_reads <- B2$read_count[match(pred$trapID,B2$trapID)]
pred$p_reads[is.na(pred$p_reads)] <- 0
pred$ph_ratio <- pred$p_reads/pred$h_reads

p1 <- box_plot(sapro, "Saprophage community")
p2 <- box_plot(phyto, "Phytophage community")
p3 <- box_plot(pred, "Predator community")

ggsave(file = "../figs/Fig_parasitoid_load.jpg",
       width = 12,
       height = 4,
       plot = p1 + p2 + p3 +
           plot_layout(axis_titles="collect_y", guides="collect", ncol=3) +
           plot_annotation(tag_levels="A") &
           theme(legend.position="bottom")
       )
