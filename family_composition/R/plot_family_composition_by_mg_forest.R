library(ggplot2)
library(patchwork)
#library(ggpubr)
#library(grid)
library(tidyverse)
library(shadowtext)
library(scales)
library(cowplot)

source("plot_family_fxn.R")


# Prepare stacked bar chart for totals
# ------------------------------------

# Get counts data
malaise_mg <- readRDS("../../iba_data/cluster_counts_malaise_long_mg.rds")
litter_mg <- readRDS("../../iba_data/cluster_counts_litter_long_mg.rds")
combined_mg <- rbind(malaise_mg,litter_mg)

# Add trap habitat
meta_mg <- read.delim("../../iba_data/malaise_litter_sample_meta_mg.tsv")
combined_mg$trap_habitat <- meta_mg$trap_habitat[match(combined_mg$sampleID_NGI,meta_mg$sampleID_NGI)]
idx <- which(combined_mg$trap_habitat=="Montane_Rainforest")
combined_mg$trap_habitat[idx] <- "Rainforest"

# Get forest contributions
get_otus_by_habitat <- function(df) {
    x <- numeric()
    set1 <- unique(df$cluster[df$trap_habitat=="Rainforest"])
    set2 <- unique(df$cluster[df$trap_habitat=="Dry_Forest"])
    x[3] <- sum(set1 %in% set2)
    x[1] <- length(set1) - x[3]
    x[2] <- length(set2) - x[3]
    return (x)
}
D <- data.frame(list(habitat=c("Rainforest","Dry forest","Shared"),
                     OTUs=get_otus_by_habitat(combined_mg),
                     dataset=rep("Combined",times=3)
                    )
                )
D$habitat <- factor(D$habitat,levels=c("Dry forest","Shared","Rainforest"))

# Render the plot
p1 <- ggplot(data=D, aes(x=dataset, y=OTUs, group=habitat, fill=habitat)) +
        geom_col(width=1.0) +
        coord_flip() +
        scale_fill_manual(name = "Habitat",
                          values = c("Dry forest"="lightgreen","Shared"="gray","Rainforest"="forestgreen"),
                          breaks = c("Rainforest","Shared","Dry forest"),
                          labels = c("Rainforest"="Rainforest","Shared"="Shared","Dry forest"="Dry forest")) +
        
        # Set basic theme
        theme_minimal(base_size=17) +

        # Get rid of title and category label on y axis (this refers to original x axis, go figure...)
        theme(axis.text.y = element_blank()) +
        theme(axis.title.y = element_blank()) +
        theme(legend.position = "top") +

        # Customize y axis (now displayed as x axis because of coord_flip)
        scale_y_continuous(
            name = "Number of OTUs",
            limits = c(0, sum(D$OTUs)),
            breaks = seq(0, 50000, by = 10000)
            ) +
 
        # Set the titles
        labs(title = NULL,
             fill = "Habitat",
             tag = "A") +
       
        # This refers to original y axis...
        theme(axis.title.x = element_text(size=17, margin=margin(t=15, r=0, b=0, l=0)))


# Prepare stacked bar chart for Madagascan radiations
# ---------------------------------------------------

# Read in and prepare data
P <- readRDS("../../placements/data/placement_stats.rds")
X <- table(P$edge_num)
P <- P[P$edge_num %in% names(X[X>4]),]
T <- readRDS("../../iba_data/cluster_taxonomy_mg.rds")
P$cluster <- T$cluster[match(P$cluster_rep,T$ASV)]

# Get forest contributions
get_radiations_by_habitat <- function(P, C) {
    x <- numeric()
    set1 <- unique(P$edge_num[P$cluster %in% C$cluster[C$trap_habitat=="Rainforest"]])
    set2 <- unique(P$edge_num[P$cluster %in% C$cluster[C$trap_habitat=="Dry_Forest"]])
    x[3] <- sum(set1 %in% set2)
    x[1] <- length(set1) - x[3]
    x[2] <- length(set2) - x[3]
    return (x)
}
D <- data.frame(list(habitat=c("Rainforest","Dry forest","Shared"),
                     Radiations=get_radiations_by_habitat(P, combined_mg),
                     dataset=rep("Combined",times=3)
                    )
                )
D$habitat <- factor(D$habitat,levels=c("Dry forest","Shared","Rainforest"))

# Render the plot
p2 <- ggplot(data=D, aes(x=dataset, y=Radiations, group=habitat, fill=habitat)) +
        geom_col(width=1.0) +
        coord_flip() +
        scale_fill_manual(name = "Habitat",
                          values = c("Dry forest"="lightgreen","Shared"="gray","Rainforest"="forestgreen"),
                          breaks = c("Rainforest","Shared","Dry forest"),
                          labels = c("Rainforest"="Rainforest","Shared"="Shared","Dry forest"="Dry forest")) +

        # Set basic theme
        theme_minimal(base_size=17) +

        # Get rid of title and category label on y axis (this refers to original x axis, go figure...)
        theme(axis.text.y = element_blank()) +
        theme(axis.title.y = element_blank()) +
        theme(legend.position = "none") +

        # Customize y axis (now displayed as x axis because of coord_flip)
        scale_y_continuous(
            name = "Number of radiations",
            limits = c(0, sum(D$Radiations)),
            breaks = seq(0, 1250, by = 250)
            ) +

        # Set the titles
        labs(title = NULL,
             fill = NULL,
             tag = "B") +

        # This refers to original y axis...
        theme(axis.title.x = element_text(size=17, margin=margin(t=15, r=0, b=0, l=0)))


# Prepare family composition bar charts
# -------------------------------------

# Add family to dataset
tax_mg <- readRDS("../../iba_data/cluster_taxonomy_mg.rds")
combined_mg$Family <- tax_mg$Family[match(combined_mg$cluster,tax_mg$cluster)]

# Get number of OTUs
mg_rainforest_otus <- unique(combined_mg[combined_mg$trap_habitat=="Rainforest",c("cluster","Family")])
mg_dryforest_otus <- unique(combined_mg[combined_mg$trap_habitat=="Dry_Forest",c("cluster","Family")])

family_stats <- function(df) {
    res <- data.frame(table(df$Family))
    colnames(res) <- c("Family","OTUs")
    res <- res[order(res$OTUs,decreasing=TRUE),]
    return (res)
}

mg_rainforest_families <- family_stats(mg_rainforest_otus)
mg_dryforest_families <- family_stats(mg_dryforest_otus)

# Add order info
mg_fams <- unique(tax_mg[,c("Family","Order")])
mg_rainforest_families$Order <- mg_fams$Order[match(mg_rainforest_families$Family,mg_fams$Family)]
mg_dryforest_families$Order <- mg_fams$Order[match(mg_dryforest_families$Family,mg_fams$Family)]
big_five <- c("Diptera","Hymenoptera","Coleoptera","Lepidoptera","Hemiptera")
idx <- which(!(mg_rainforest_families$Order %in% big_five))
mg_rainforest_families$Order[idx] <- "Other"
idx <- which(!(mg_dryforest_families$Order %in% big_five))
mg_dryforest_families$Order[idx] <- "Other"

# Order families
mg_rainforest_families <- mg_rainforest_families[order(mg_rainforest_families$OTUs,decreasing=TRUE),]
mg_dryforest_families <- mg_dryforest_families[order(mg_dryforest_families$OTUs,decreasing=TRUE),]

# Generate label with rainforest order before family name
mg_rainforest_families$order_num <- 1:nrow(mg_rainforest_families)
mg_rainforest_families$family_label <- paste0(mg_rainforest_families$order_num,". ",mg_rainforest_families$Family)
mg_dryforest_families$family_label <- mg_rainforest_families$family_label[match(mg_dryforest_families$Family,mg_rainforest_families$Family)]

# Use augmented labels as family names
mg_rainforest_families$Family <- mg_rainforest_families$family_label
mg_dryforest_families$Family <- mg_dryforest_families$family_label

# Order as factor with levels in order
orders <- c(big_five, "Other")
mg_rainforest_families$Order <- factor(mg_rainforest_families$Order, levels=orders)
mg_dryforest_families$Order <- factor(mg_dryforest_families$Order, levels=orders)


# Generate plots
# --------------

# Make plots

p3 <- plot_family(mg_rainforest_families, 20, 9000, 6, "Rainforest OTUs", 0.3) + guides(fill=guide_legend(nrow=1))
p4 <- plot_family(mg_dryforest_families, 20, 3500, 7, "Dry forest OTUs", 0.3) + guides(fill=guide_legend(nrow=1))
p3 <- p3 + labs(tag = "C")
p4 <- p4 + labs(tag = "D")

# Extract and remove legend to display one shared legend
# Show "Other" as "Orthoptera", as only Orthoptera included in top hits
pDummy <- plot_family(mg_rainforest_families, 20, 9000, 6, "Rainforest", 0.3, other_orthoptera=TRUE) + guides(fill=guide_legend(nrow=1))
p5 <- cowplot::get_legend(pDummy)   # Extract legend

p3 <- p3 + theme(legend.position="none") # Remove legend
p4 <- p4 + theme(legend.position="none") # Remove legend


# Put plots together
ggsave(
       file = "../figs/Fig_family_composition_mg_forests.jpg",
       width = 14.0,
       height = 20.0,
       plot = p1 / p2 / (p3 + p4) / p5 +
              plot_layout(axis_titles="collect", heights=c(0.05,0.05,0.85,0.05))
       )

