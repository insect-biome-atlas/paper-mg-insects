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
prop_estimate <- function(D, taxon, sum_total) { D$species[D$order==taxon] / sum_total }

taxon_comp <- read.delim("../data/taxon_comp.tsv")
sum_total <- sum(taxon_comp$species)
dipt_estimate <- prop_estimate(taxon_comp, "Diptera", sum_total)
hyme_estimate <- prop_estimate(taxon_comp, "Hymenoptera", sum_total)
cole_estimate <- prop_estimate(taxon_comp, "Coleoptera", sum_total)
lepi_estimate <- prop_estimate(taxon_comp, "Lepidoptera", sum_total)
hemi_estimate <- prop_estimate(taxon_comp, "Hemiptera", sum_total)
other_estimate <- prop_estimate(taxon_comp, "Other", sum_total)

taxon_comp <- read.delim("../data/new_taxon_comp.tsv")
sum_total <- sum(taxon_comp$species)
new_dipt_estimate <- prop_estimate(taxon_comp, "Diptera", sum_total)
new_hyme_estimate <- prop_estimate(taxon_comp, "Hymenoptera", sum_total)
new_cole_estimate <- prop_estimate(taxon_comp, "Coleoptera", sum_total)
new_lepi_estimate <- prop_estimate(taxon_comp, "Lepidoptera", sum_total)
new_hemi_estimate <- prop_estimate(taxon_comp, "Hemiptera", sum_total)
new_other_estimate <- prop_estimate(taxon_comp, "Other", sum_total)

# Render plots
plot_prop <- function(D, estimated_val1, estimated_val2, plot_title) {
    ggplot(data=D, aes(x=samples, y=taxon_otu_prop)) +
        ggtitle(plot_title) +
        geom_point(alpha=.02 , size=2 , show.legend=FALSE, colour="blue") +
        geom_smooth(method="loess", se=TRUE, lwd=1) +
        geom_hline(yintercept=estimated_val1, colour="red", linetype="dotted", linewidth=1) +
        geom_hline(yintercept=estimated_val2, colour="red", linetype="dashed", linewidth=1) +
        theme_linedraw(base_size=20) +
        ylim(c(0.0,0.8)) +
        labs(x = "No. sites" , y = "Proportion of OTUs")
}

p1 <- plot_prop(dipt, dipt_estimate, new_dipt_estimate, "Diptera")
p2 <- plot_prop(hyme, hyme_estimate, new_hyme_estimate, "Hymenoptera")
p3 <- plot_prop(cole, cole_estimate, new_cole_estimate, "Coleoptera")
p4 <- plot_prop(lepi, lepi_estimate, new_lepi_estimate, "Lepidoptera")
p5 <- plot_prop(hemi, hemi_estimate, new_hemi_estimate, "Hemiptera")
p6 <- plot_prop(other, other_estimate, new_other_estimate, "Other")

# Save plots
ggsave( "../figs/Fig_site_acc_taxonomic_complete_se.jpg",
        width=21,
        height=14,
        plot=p1 + p2 + p3 + p4 + p5 + p6 +
            plot_layout(axis_titles="collect", ncol=3) +
            plot_annotation(tag_levels="A")
      )

