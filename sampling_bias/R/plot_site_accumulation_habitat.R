# Plot site accumulation data for habitats

library(ggplot2)
library(patchwork)


# Read in data
plants <- read.delim("../data/site_acc_plants_prop.tsv")
soil <- read.delim("../data/site_acc_soil_prop.tsv")
water <- read.delim("../data/site_acc_water_prop.tsv")
wood <- read.delim("../data/site_acc_wood_prop.tsv")
temporary <- read.delim("../data/site_acc_temporary_prop.tsv")
fungi <- read.delim("../data/site_acc_fungi_prop.tsv")

# Compute estimates
habitat_comp <- read.delim("../data/habitat_comp.tsv")
habitats <- c("Plants","Soil","Water","Wood","Temporary habitats","Fungi")
sum_major <- sum(habitat_comp$species[habitat_comp$habitat %in% habitats])
prop_estimate <- function(habitat) { habitat_comp$species[habitat_comp$habitat==habitat] / sum_major }
plants_estimate <- prop_estimate("Plants")
soil_estimate <- prop_estimate("Soil")
water_estimate <- prop_estimate("Water")
wood_estimate <- prop_estimate("Wood")
temporary_estimate <- prop_estimate("Temporary habitats")
fungi_estimate <- prop_estimate("Fungi")

# Render plots
plot_prop <- function(D, estimated_val, plot_title) {
    ggplot(data=D, aes(x=samples, y=habitat_otu_prop)) +
        ggtitle(plot_title) +
        geom_point(alpha=.02 , size=2 , show.legend=FALSE, colour="blue") +
        geom_smooth(method="loess", se=TRUE, lwd=1) +
        geom_hline(yintercept=estimated_val, colour="red", linetype="dotted", linewidth=1) +
        theme_linedraw(base_size=20) +
        ylim(c(0.0,0.8)) +
        labs(x = "No. sites" , y = "Proportion of OTUs")
}

p1 <- plot_prop(plants, plants_estimate, "Plants")
p2 <- plot_prop(soil, soil_estimate, "Soil")
p3 <- plot_prop(water, water_estimate, "Water")
p4 <- plot_prop(wood, wood_estimate, "Wood")
p5 <- plot_prop(temporary, temporary_estimate, "Temporary habitats")
p6 <- plot_prop(fungi, fungi_estimate, "Fungi")

# Save plots
ggsave( "../figs/Fig_site_acc_habitat_complete_se.jpg",
        width=14,
        height=21,
        plot=p1 + p2 + p3 + p4 + p5 + p6 +
            plot_layout(axis_titles="collect", ncol=2) +
            plot_annotation(tag_levels="A")
      )

