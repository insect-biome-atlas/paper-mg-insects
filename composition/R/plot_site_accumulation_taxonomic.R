# Plot site accumulation data for taxonomic groups

library(ggplot2)
library(patchwork)
source("../../fig_settings/fig_colours.R")


# Read in data
dipt <- read.delim("../data/site_acc_diptera_prop.tsv")
hyme <- read.delim("../data/site_acc_hymenoptera_prop.tsv")
cole <- read.delim("../data/site_acc_coleoptera_prop.tsv")
lepi <- read.delim("../data/site_acc_lepidoptera_prop.tsv")
hemi <- read.delim("../data/site_acc_hemiptera_prop.tsv")
other <- read.delim("../data/site_acc_other_orders_prop.tsv")


# Render plots
p1 <- ggplot(data=dipt, aes(x=samples, y=taxon_otu_prop, group=country)) +
  ggtitle("Diptera") +
  geom_point(alpha=.05 , size=2 , show.legend=FALSE, aes(colour=country)) +
  theme_linedraw(base_size=20) +
  geom_smooth(method="loess", se=TRUE, lwd=1, show.legend=FALSE, aes(colour=country)) +
  xlim(c(1,50)) +
  ylim(c(0.0,0.65)) +
  scale_color_manual(values=c(mg_col,se_col)) +
  labs(x = "No. sites" , y = "Proportion of OTUs")

p2 <- ggplot(data=hyme, aes(x=samples, y=taxon_otu_prop, group=country)) +
  ggtitle("Hymenoptera") +
  geom_point(alpha=.05 , size=2 , show.legend=FALSE, aes(colour=country)) +
  theme_linedraw(base_size=20) +
  geom_smooth(method="loess", se=TRUE, lwd=1, show.legend=FALSE, aes(colour=country)) +
  xlim(c(1,50)) +
  ylim(c(0.0,0.65)) +
  scale_color_manual(values=c(mg_col,se_col)) +
  labs(x = "No. sites" , y = "Proportion of OTUs")

p3 <- ggplot(data=cole, aes(x=samples, y=taxon_otu_prop, group=country)) +
  ggtitle("Coleoptera") +
  geom_point(alpha=.05 , size=2 , show.legend=FALSE, aes(colour=country)) +
  theme_linedraw(base_size=20) +
  geom_smooth(method="loess", se=TRUE, lwd=1, show.legend=FALSE, aes(colour=country)) +
  xlim(c(1,50)) +
  ylim(c(0.0,0.65)) +
  scale_color_manual(values=c(mg_col,se_col)) +
  labs(x = "No. sites" , y = "Proportion of OTUs")

p4 <- ggplot(data=lepi, aes(x=samples, y=taxon_otu_prop, group=country)) +
  ggtitle("Lepidoptera") +
  geom_point(alpha=.05 , size=2 , show.legend=FALSE, aes(colour=country)) +
  theme_linedraw(base_size=20) +
  geom_smooth(method="loess", se=TRUE, lwd=1, show.legend=FALSE, aes(colour=country)) +
  xlim(c(1,50)) +
  ylim(c(0.0,0.65)) +
  scale_color_manual(values=c(mg_col,se_col)) +
  labs(x = "No. sites" , y = "Proportion of OTUs")

p5 <- ggplot(data=hemi, aes(x=samples, y=taxon_otu_prop, group=country)) +
  ggtitle("Hemiptera") +
  geom_point(alpha=.05 , size=2 , show.legend=FALSE, aes(colour=country)) +
  theme_linedraw(base_size=20) +
  geom_smooth(method="loess", se=TRUE, lwd=1, show.legend=FALSE, aes(colour=country)) +
  xlim(c(1,50)) +
  ylim(c(0.0,0.65)) +
  scale_color_manual(values=c(mg_col,se_col)) +
  labs(x = "No. sites" , y = "Proportion of OTUs")

p6 <- ggplot(data=other, aes(x=samples, y=taxon_otu_prop, group=country)) +
  ggtitle("Other") +
  geom_point(alpha=.05 , size=2 , show.legend=FALSE, aes(colour=country)) +
  theme_linedraw(base_size=20) +
  geom_smooth(method="loess", se=TRUE, lwd=1, show.legend=FALSE, aes(colour=country)) +
  xlim(c(1,50)) +
  ylim(c(0.0,0.65)) +
  scale_color_manual(values=c(mg_col,se_col)) +
  labs(x = "No. sites" , y = "Proportion of OTUs")

# Print Plots
ggsave("../figs/Fig_site_acc_taxonomic.jpg",
       width = 21,
       height = 14,
       plot = p1 + p2 + p3 + p4 + p5 + p6 +
              plot_layout(axis_titles="collect", ncol=3) +
              plot_annotation(tag_levels="A")
       )

