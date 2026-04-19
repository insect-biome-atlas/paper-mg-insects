# Plot site accumulation data for niches

library(ggplot2)
library(patchwork)


# Read in data
sapro <- read.delim("../data/site_acc_saprophage_prop.tsv")
phyto <- read.delim("../data/site_acc_phytophage_prop.tsv")
pred <- read.delim("../data/site_acc_predator_prop.tsv")
sapro_ph <- read.delim("../data/site_acc_saprophage_ph_ratio.tsv")
phyto_ph <- read.delim("../data/site_acc_phytophage_ph_ratio.tsv")
pred_ph <- read.delim("../data/site_acc_predator_ph_ratio.tsv")

# Compute estimates
niche_comp <- read.delim("../data/niche_comp.tsv")
major_niches <- c("Saprophage","Phytophage","Predator")
sum_major <- sum(niche_comp$species[niche_comp$niche %in% major_niches])
prop_estimate <- function(niche) { niche_comp$species[niche_comp$niche==niche] / sum_major }
sapro_estimate <- prop_estimate("Saprophage")
phyto_estimate <- prop_estimate("Phytophage")
pred_estimate <- prop_estimate("Predator")
ph_ratio_estimate <- function(parasitoid, host) { niche_comp$species[niche_comp$niche==parasitoid] / niche_comp$species[niche_comp$niche==host] }
sapro_ph_estimate <- ph_ratio_estimate("Saprophage-parasitoid","Saprophage")
phyto_ph_estimate <- ph_ratio_estimate("Phytophage-parasitoid","Phytophage")
pred_ph_estimate <- ph_ratio_estimate("Predator-parasitoid","Predator")

# Render plots
plot_prop <- function(D, estimated_val, plot_title) {
    ggplot(data=D, aes(x=samples, y=niche_otu_prop)) +
        ggtitle(plot_title) +
        geom_point(alpha=.02 , size=2 , show.legend=FALSE, colour="blue") +
        geom_smooth(method="loess", se=TRUE, lwd=1) +
        geom_hline(yintercept=estimated_val, colour="red", linetype="dotted", linewidth=1) +
        theme_linedraw(base_size=20) +
        ylim(c(0.0,0.8)) +
        labs(x = "No. sites" , y = "Proportion of major-niche OTUs")
}

plot_ph_ratio <- function(D, estimated_val, plot_title) {
    ggplot(data=D, aes(x=samples, y=ph_otu_ratio)) +
        ggtitle(plot_title) +
        geom_point(alpha=.02 , size=2 , show.legend=FALSE, colour="blue") +
        geom_smooth(method="loess", se=TRUE, lwd=1) +
        geom_hline(yintercept=estimated_val, colour="red", linetype="dotted", linewidth=1) +
        theme_linedraw(base_size=20) +
        ylim(c(0.0,1.05)) +
        labs(x = "No. sites" , y = "Parasitoid:host OTU ratio")
}

p1 <- plot_prop(sapro, sapro_estimate, "Saprophages")
p2 <- plot_ph_ratio(sapro_ph, sapro_ph_estimate, "Saprophage parasitoids")
p3 <- plot_prop(phyto, phyto_estimate, "Phytophages")
p4 <- plot_ph_ratio(phyto_ph, phyto_ph_estimate, "Phytophage parasitoids")
p5 <- plot_prop(pred, pred_estimate, "Predators")
p6 <- plot_ph_ratio(pred_ph, pred_ph_estimate, "Predator parasitoids")


# Save plots
ggsave( "../figs/Fig_site_acc_niche_complete_se.jpg",
        width=14,
        height=21,
        plot=p1 + p2 + p3 + p4 + p5 + p6 +
            plot_layout(axis_titles="collect", ncol=2) +
            plot_annotation(tag_levels="A")
      )

