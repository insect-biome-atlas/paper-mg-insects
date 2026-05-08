# Plot sample accumulation data for habitats

library(ggplot2)
library(patchwork)

# Read in data
plants <- read.delim("../data/sample_acc_plants_prop.tsv")
soil <- read.delim("../data/sample_acc_soil_prop.tsv")
water <- read.delim("../data/sample_acc_water_prop.tsv")
wood <- read.delim("../data/sample_acc_wood_prop.tsv")
temporary <- read.delim("../data/sample_acc_temporary_prop.tsv")
fungi <- read.delim("../data/sample_acc_fungi_prop.tsv")

# Render plots
p1 <- ggplot(data=plants, aes(x=samples, y=habitat_otu_prop, group=country)) +
  ggtitle("Plants") +
  geom_point(alpha=.05 , size=2 , show.legend=FALSE, aes(colour=country)) +
  theme_linedraw(base_size=20) +
  geom_smooth(method="loess", se=TRUE, lwd=1) +
  xlim(c(1,40)) +
  scale_color_manual(values=c("Green","Blue")) +
  labs(x = "No. samples" , y = "Proportion of OTUs")

p2 <- ggplot(data=soil, aes(x=samples, y=habitat_otu_prop, group=country)) +
  ggtitle("Soil") +
  geom_point(alpha=.05 , size=2 , show.legend=FALSE, aes(colour=country)) +
  theme_linedraw(base_size=20) +
  geom_smooth(method="loess", se=TRUE, lwd=1) +
  xlim(c(1,40)) +
  scale_color_manual(values=c("Green","Blue")) +
  labs(x = "No. samples" , y = "Proportion of OTUs")

p3 <- ggplot(data=water, aes(x=samples, y=habitat_otu_prop, group=country)) +
  ggtitle("Water") +
  geom_point(alpha=.05 , size=2 , show.legend=FALSE, aes(colour=country)) +
  theme_linedraw(base_size=20) +
  geom_smooth(method="loess", se=TRUE, lwd=1) +
  xlim(c(1,40)) +
  scale_color_manual(values=c("Green","Blue")) +
  labs(x = "No. samples" , y = "Proportion of OTUs")

p4 <- ggplot(data=wood, aes(x=samples, y=habitat_otu_prop, group=country)) +
  ggtitle("Wood") +
  geom_point(alpha=.05 , size=2 , show.legend=FALSE, aes(colour=country)) +
  theme_linedraw(base_size=20) +
  geom_smooth(method="loess", se=TRUE, lwd=1) +
  xlim(c(1,40)) +
  scale_color_manual(values=c("Green","Blue")) +
  labs(x = "No. samples" , y = "Proportion of OTUs")

p5 <- ggplot(data=temporary, aes(x=samples, y=habitat_otu_prop, group=country)) +
  ggtitle("Temporary habitats") +
  geom_point(alpha=.05 , size=2 , show.legend=FALSE, aes(colour=country)) +
  theme_linedraw(base_size=20) +
  geom_smooth(method="loess", se=TRUE, lwd=1) +
  xlim(c(1,40)) +
  scale_color_manual(values=c("Green","Blue")) +
  labs(x = "No. samples" , y = "Proportion of OTUs")

p6 <- ggplot(data=fungi, aes(x=samples, y=habitat_otu_prop, group=country)) +
  ggtitle("Fungi") +
  geom_point(alpha=.05 , size=2 , show.legend=FALSE, aes(colour=country)) +
  theme_linedraw(base_size=20) +
  geom_smooth(method="loess", se=TRUE, lwd=1) +
  xlim(c(1,40)) +
  scale_color_manual(values=c("Green","Blue")) +
  labs(x = "No. samples" , y = "Proportion of OTUs")

# Print Plots
ggsave("../figs/Fig_sample_acc_habitat.jpg",
       width=21,
       height=14,
       plot=p1 + p2 + p3 + p4 + p5 + p6 +
       plot_layout(axis_titles="collect",ncol=3) +
       plot_annotation(tag_levels="A")
       )

