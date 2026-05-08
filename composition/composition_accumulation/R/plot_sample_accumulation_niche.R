# Plot sample accumulation data for niches

library(ggplot2)
library(patchwork)


# Read in data
sapro <- read.delim("../data/sample_acc_saprophage_prop.tsv")
phyto <- read.delim("../data/sample_acc_phytophage_prop.tsv")
pred <- read.delim("../data/sample_acc_predator_prop.tsv")
sapro_ph <- read.delim("../data/sample_acc_saprophage_ph_ratio.tsv")
phyto_ph <- read.delim("../data/sample_acc_phytophage_ph_ratio.tsv")
pred_ph <- read.delim("../data/sample_acc_predator_ph_ratio.tsv")


# Render plots
p1 <- ggplot(data=sapro, aes(x=samples, y=niche_otu_prop, group=country)) +
  ggtitle("Saprophages") +
  geom_point(alpha=.05 , size=2 , show.legend=FALSE, aes(colour=country)) +
  theme_linedraw(base_size=20) +
  geom_smooth(method="loess", se=TRUE, lwd=1) +
  xlim(c(1,40)) +
  scale_color_manual(values=c("Green","Blue")) +
  labs(x = "No. samples" , y = "Proportion of major-niche OTUs")

p2 <- ggplot(data=phyto, aes(x=samples, y=niche_otu_prop, group=country)) +
  ggtitle("Phytophages") +
  geom_point(alpha=.05 , size=2 , show.legend=FALSE, aes(colour=country)) +
  theme_linedraw(base_size=20) +
  geom_smooth(method="loess", se=TRUE, lwd=1) +
  xlim(c(1,40)) +
  scale_color_manual(values=c("Green","Blue")) +
  labs(x = "No. samples" , y = "Proportion of major-niche OTUs")

p3 <- ggplot(data=pred, aes(x=samples, y=niche_otu_prop, group=country)) +
  ggtitle("Predators") +
  geom_point(alpha=.05 , size=2 , show.legend=FALSE, aes(colour=country)) +
  theme_linedraw(base_size=20) +
  geom_smooth(method="loess", se=TRUE, lwd=1) +
  xlim(c(1,40)) +
  scale_color_manual(values=c("Green","Blue")) +
  labs(x = "No. samples" , y = "Proportion of major-niche OTUs")

p4 <- ggplot(data=sapro_ph, aes(x=samples, y=ph_otu_ratio, group=country)) +
  ggtitle("Saprophage parasitoids") +
  geom_point(alpha=.05 , size=2 , show.legend=FALSE, aes(colour=country)) +
  theme_linedraw(base_size=20) +
  geom_smooth(method="loess", se=TRUE, lwd=1) +
  xlim(c(1,40)) +
  scale_color_manual(values=c("Green","Blue")) +
  labs(x = "No. samples" , y = "Parasitoid:host OTU ratio")

p5 <- ggplot(data=phyto_ph, aes(x=samples, y=ph_otu_ratio, group=country)) +
  ggtitle("Phytophage parasitoids") +
  geom_point(alpha=.05 , size=2 , show.legend=FALSE, aes(colour=country)) +
  theme_linedraw(base_size=20) +
  geom_smooth(method="loess", se=TRUE, lwd=1) +
  xlim(c(1,40)) +
  scale_color_manual(values=c("Green","Blue")) +
  labs(x = "No. samples" , y = "Parasitoid:host OTU ratio")

p6 <- ggplot(data=pred_ph, aes(x=samples, y=ph_otu_ratio, group=country)) +
  ggtitle("Predator parasitoids") +
  geom_point(alpha=.05 , size=2 , show.legend=FALSE, aes(colour=country)) +
  theme_linedraw(base_size=20) +
  geom_smooth(method="loess", se=TRUE, lwd=1) +
  xlim(c(1,40)) +
  scale_color_manual(values=c("Green","Blue")) +
  labs(x = "No. samples" , y = "Parasitoid:host OTU ratio")

# Print Plots
ggsave("../figs/Fig_sample_acc_niche.jpg",
       width=21,
       height=14,
       plot=p1 + p2 + p3 + p4 + p5 + p6 +
       plot_layout(axis_titles="collect",ncol=3) +
       plot_annotation(tag_levels="A")
       )

