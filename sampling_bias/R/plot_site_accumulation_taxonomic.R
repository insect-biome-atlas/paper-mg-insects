# Plot site accumulation data for taxonomic groups

library(ggplot2)
library(patchwork)


# Read in data
dipt <- read.delim("../data/site_acc_diptera_prop.tsv")
hyme <- read.delim("../data/site_acc_hymenoptera_prop.tsv")
cole <- read.delim("../data/site_acc_coleoptera_prop.tsv")
lepi <- read.delim("../data/site_acc_lepidoptera_prop.tsv")
hemi <- read.delim("../data/site_acc_hemiptera_prop.tsv")
other <- read.delim("../data/site_acc_other_orders_prop.tsv")

# Compute estimates
taxon_comp <- read.delim("../data/taxon_comp.tsv")
sum_total <- sum(taxon_comp$species)
prop_estimate <- function(taxon) { taxon_comp$species[taxon_comp$order==taxon] / sum_total }
dipt_estimate <- prop_estimate("Diptera")
hyme_estimate <- prop_estimate("Hymenoptera")
cole_estimate <- prop_estimate("Coleoptera")
lepi_estimate <- prop_estimate("Lepidoptera")
hemi_estimate <- prop_estimate("Hemiptera")
other_estimate <- prop_estimate("Other")

# Render plots
plot_prop <- function(D, estimated_val, plot_title) {
    ggplot(data=D, aes(x=samples, y=taxon_otu_prop)) +
        ggtitle(plot_title) +
        geom_point(alpha=.02 , size=2 , show.legend=FALSE, colour="blue") +
        geom_smooth(method="loess", se=TRUE, lwd=1) +
        geom_hline(yintercept=estimated_val, colour="red", linetype="dotted", linewidth=1) +
        theme_linedraw(base_size=20) +
        ylim(c(0.0,0.8)) +
        labs(x = "No. sites" , y = "Proportion of OTUs")
}

p1 <- plot_prop(dipt, dipt_estimate, "Diptera")
p2 <- plot_prop(hyme, hyme_estimate, "Hymenoptera")
p3 <- plot_prop(cole, cole_estimate, "Coleoptera")
p4 <- plot_prop(lepi, lepi_estimate, "Lepidoptera")
p5 <- plot_prop(hemi, hemi_estimate, "Hemiptera")
p6 <- plot_prop(other, other_estimate, "Other")

# Save plots
ggsave( "../figs/Fig_site_acc_taxonomic_complete_se.jpg",
        width=21,
        height=14,
        plot=p1 + p2 + p3 + p4 + p5 + p6 +
            plot_layout(axis_titles="collect", ncol=3) +
            plot_annotation(tag_levels="A")
      )

