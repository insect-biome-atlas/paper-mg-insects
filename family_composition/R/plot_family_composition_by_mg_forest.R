library(ggplot2)
library(patchwork)
#library(ggpubr)
#library(grid)
library(tidyverse)
library(shadowtext)
library(scales)

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


# Make the plot
# -------------

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
             fill = "Habitat") +
       
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

# Select the top 20 families
mg_rainforest_families <- mg_rainforest_families[order(mg_rainforest_families$OTUs,decreasing=TRUE),]
mg_dryforest_families <- mg_dryforest_families[order(mg_dryforest_families$OTUs,decreasing=TRUE),]
mg_rainforest_families <- mg_rainforest_families[1:20,]
mg_dryforest_families <- mg_dryforest_families[1:20,]

# Order as factor with levels in order
orders <- c(big_five, "Other")
mg_rainforest_families$Order <- factor(mg_rainforest_families$Order, levels=orders)
mg_dryforest_families$Order <- factor(mg_dryforest_families$Order, levels=orders)


# Generate plots
# --------------

my_colors <- viridis_pal()(6)

# Simple plot rendering function
plot_family <- function(D, num_top, max_tick, div, plot_title) {
   
    D <- D[1:num_top,]
    D$Family <- factor(D$Family,levels=c(D$Family[num_top:1]))
    max_x <- max(D$OTUs)
    D$label_cutoff <- max_x * as.numeric(sapply(as.character(D$Family),nchar)) / 40

    ggplot(data=D, aes(x=OTUs, y=Family, fill=Order)) +
        geom_col(width=0.8) +
        theme_minimal(base_size=17) +
        labs(title = plot_title,
             x = "Number of OTUs") +
        scale_fill_manual(
            name = "Order",
            values = c("Diptera" = my_colors[1], "Hymenoptera" = my_colors[2], "Coleoptera" = my_colors[3],
                       "Lepidoptera" = my_colors[4], "Hemiptera" = my_colors[5], "Other" = my_colors[6])) +
        theme(legend.position="bottom") +
        scale_x_continuous(
#            limits = c(0, max_x),
            breaks = seq(0, max_tick, by=max_tick/div),
            position = "top"
        ) +
        # The vertical axis extends upwards and downwards
        scale_y_discrete(expand = expansion(add = c(0.5, 0.5))) +
        theme(
            # Set the color and the width of the grid lines for the horizontal axis
            panel.grid.major.x = element_line(color = "#A8BAC4", linewidth = 0.3),
            panel.grid.minor.x = element_line(color = "white", linewidth = 0.3),
            panel.grid.major.y = element_line(color = "white", linewidth = 0.3),
            # Remove tick marks by setting their length to 0
            axis.ticks.length = unit(0, "mm"),
            # Remove the title for the y axis, keep the x axis label
            axis.title.y = element_text(color="white"),
            axis.title.x = element_text(color="white"),
            # Remove labels from the vertical axis
            axis.text.y = element_blank(),
            # But customize labels for the horizontal axis
            axis.text.x = element_text(family = "Helvetica", size = 14)
        ) +
        # Add labels back in, into the plot
        geom_shadowtext(
            data = subset(D, OTUs < label_cutoff),
            aes(OTUs, y = Family, label = Family),
            hjust = 0,
            nudge_x = max_x/100,
            colour = "black",
            bg.colour = "white",
            bg.r = 0.2,
            family = "Helvetica",
            size = 6) +
        geom_text(
            data = subset(D, OTUs >= label_cutoff),
            aes(0.3*max_x, y = Family, label = Family),
            hjust = 0,
            nudge_x = max_x/100,
            colour = "white",
            family = "Helvetica",
            size = 6)
}


# Make plots
p2 <- plot_family(mg_rainforest_families, 20, 9000, 6, "Rainforest")
p3 <- plot_family(mg_dryforest_families, 20, 3500, 7, "Dry forest")


# Put plots together
ggsave(
       file = "../figs/Fig_family_composition_mg_forests.jpg",
       width = 14.0,
       height = 21.0,
       plot = p1 / (p2 + p3) +
              plot_layout(axis_titles="collect", heights=c(0.05,1.0)) +
              plot_annotation(tag_levels="A")
       )

