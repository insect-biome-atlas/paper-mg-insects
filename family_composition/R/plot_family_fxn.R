library(ggplot2)
library(patchwork)
#library(ggpubr)
#library(grid)
library(tidyverse)
library(shadowtext)
library(scales)


# Family plot rendering function
plot_family <- function(D, num_top, max_tick, div, plot_title, offset_frac=0.0, other_orthoptera=FALSE) {

    my_colors <- viridis_pal()(6)
    
    D <- D[1:num_top,]
    D$Family <- factor(D$Family,levels=c(D$Family[num_top:1]))
    max_x <- max(D$OTUs)

    other_label <- "Other"
    if (other_orthoptera) {
        levels(D$Order)[levels(D$Order)=="Other"] <- "Orthoptera"
        other_label <- "Orthoptera"
    }

    D$label_cutoff <- max_x * as.numeric(sapply(as.character(D$Family),nchar)) / 40

    plt <- ggplot(data=D, aes(x=OTUs, y=Family, fill=Order)) +
        geom_col(width=0.8) +
        theme_minimal(base_size=17) +
        labs(title = plot_title,
            x = "Number of OTUs") +
#        scale_fill_manual(
#            name = "Order",
#            values = c("Diptera" = my_colors[1], "Hymenoptera" = my_colors[2], "Coleoptera" = my_colors[3],
#                       "Lepidoptera" = my_colors[4], "Hemiptera" = my_colors[5], "Other" = my_colors[6]))  +
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
            aes(0 + offset_frac*max_x, y = Family, label = Family),
            hjust = 0,
            nudge_x = max_x/100,
            colour = "white",
            family = "Helvetica",
            size = 6)

    if (other_orthoptera) {
        plt <- plt +
            scale_fill_manual(
            name = "Order",
            values = c("Diptera" = my_colors[1], "Hymenoptera" = my_colors[2], "Coleoptera" = my_colors[3],
                       "Lepidoptera" = my_colors[4], "Hemiptera" = my_colors[5], "Orthoptera" = my_colors[6]))
    } else {
        plt <- plt +
            scale_fill_manual(
                name = "Order",
                values = c("Diptera" = my_colors[1], "Hymenoptera" = my_colors[2], "Coleoptera" = my_colors[3],
                           "Lepidoptera" = my_colors[4], "Hemiptera" = my_colors[5], "Other" = my_colors[6]))
    }
 
    plt

}

# Clade plot rendering function
plot_clade <- function(D, num_top, max_tick, div, plot_title, offset_frac=0.0, other_orthoptera=FALSE) {

    D$Family <- D$Clade
    plot_family(D, num_top, max_tick, div, plot_title, offset_frac, other_orthoptera)
}

