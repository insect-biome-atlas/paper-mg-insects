# Compute and plot predictive diversification tests

library(ggplot2)
library(patchwork)


# Get MG data, and samples from the predictive distribution
D <- read.delim("../data/diversification_data.tsv")
P <- read.delim("../data/predictive_samples.tsv")

# Get SE data and complement with clade info
S <- read.delim("../../composition/data/otu_site_meta_se.tsv")
T <- readRDS("../../iba_data/cluster_taxonomy_SE.rds")
S$Clade <- T$Clade[match(S$cluster,T$cluster)]

# KL divergence
kldiv <- function(p,q) {
    sum(p*log(p/q))
}

plot_fxn <- function(pred_vals, obs_val, title, xlab) {

    D <- data.frame(list(pred_vals=pred_vals))
    E <- data.frame(list(obs_val=obs_val, y=0.0))
    ggplot(D, aes(x=pred_vals)) +
        theme_minimal() + 
        geom_density(fill="steelblue", alpha=0.5) +
        geom_point(data=E, aes(x=obs_val, y=y), shape=23, size=3, fill="red") +
        ylab(NULL) +
        xlab(xlab) +
        ggtitle(title)
}

cat("\n")   # Make output more clearly separated from command line

# 1. Major niche composition
# --------------------------

# Reference distribution
niches <- c("Saprophage","Phytophage","Predator")
p <- numeric(length(niches))
for (i in 1:length(niches)) {
    p[i] <- length(unique(S$cluster[S$Niche==niches[i]])) / length(unique(S$cluster[S$Niche %in% niches]))
}

# Observed kl_divergence
q <- numeric(length(niches))
for (i in 1:length(niches)) {
    q[i] <- sum(D$OTUs[D$Niche==niches[i]]) / sum(D$OTUs[D$Niche %in% niches])
}
obs_kl_div <- kldiv(p,q)

pred_kl_div <- numeric(ncol(P))
for (i in 1:ncol(P)) {
    otus <- P[,i]
    q <- numeric(length(niches))
    for (j in 1:length(niches)) {
        q[j] <- sum(otus[D$Niche==niches[j]]) / sum(otus[D$Niche %in% niches])
    }
    pred_kl_div[i] <- kldiv(p,q)
}
tail_prob <- sum(pred_kl_div <= obs_kl_div)

cat("Observed KL divergence in major-niche proportions is:", obs_kl_div, "\n")
cat("Tail probability:", tail_prob, "\n")
cat("\n")

p1 <- plot_fxn(pred_kl_div, obs_kl_div, "Major-niche composition", "KL divergence")


# 2. Habitat composition
# ----------------------

# Reference distribution
habitats <- c("Plants","Soil","Water","Wood","Temporary habitats","Fungi")
p <- numeric(length(habitats))
for (i in 1:length(habitats)) {
    p[i] <- length(unique(S$cluster[S$Habitat==habitats[i]])) / length(unique(S$cluster[S$Habitat %in% habitats]))
}

# Observed kl_divergence
q <- numeric(length(habitats))
for (i in 1:length(habitats)) {
    q[i] <- sum(D$OTUs[D$Habitat==habitats[i]]) / sum(D$OTUs[D$Habitat %in% habitats])
}
obs_kl_div <- kldiv(p,q)

pred_kl_div <- numeric(ncol(P))
for (i in 1:ncol(P)) {
    otus <- P[,i]
    q <- numeric(length(habitats))
    for (j in 1:length(habitats)) {
        q[j] <- sum(otus[D$Habitat==habitats[j]]) / sum(otus[D$Habitat %in% habitats])
    }
    pred_kl_div[i] <- kldiv(p,q)
}
tail_prob <- sum(pred_kl_div <= obs_kl_div) / length(pred_kl_div)

cat("Observed KL divergence in habitat proportions is:", obs_kl_div, "\n")
cat("Tail probability (kl_div <= observed):", tail_prob, "\n")
cat("\n")

p2 <- plot_fxn(pred_kl_div, obs_kl_div, "Habitat composition", "KL divergence")


# 3. Order composition
# --------------------

# Reference distribution
orders <- unique(S$Order)
p <- numeric(length(orders))
for (i in 1:length(orders)) {
    p[i] <- (length(unique(S$cluster[S$Order==orders[i]]))+1) / (length(unique(S$cluster))+length(orders))
}

# Observed kl_divergence
q <- numeric(length(orders))
for (i in 1:length(orders)) {
    q[i] <- (sum(D$OTUs[D$Order==orders[i]])+1) / (sum(D$OTUs[D$Order %in% orders])+length(orders))
}
obs_kl_div <- kldiv(p,q)

pred_kl_div <- numeric(ncol(P))
for (i in 1:ncol(P)) {
    otus <- P[,i]
    q <- numeric(length(orders))
    for (j in 1:length(orders)) {
        q[j] <- (sum(otus[D$Order==orders[j]])+1) / (sum(otus[D$Order %in% orders])+length(orders))
    }
    pred_kl_div[i] <- kldiv(p,q)
}
tail_prob <- sum(pred_kl_div <= obs_kl_div)

cat("Observed KL divergence in order proportions is:", obs_kl_div, "\n")
cat("Tail probability (kl_div <= observed):", tail_prob, "\n")
cat("\n")

p3 <- plot_fxn(pred_kl_div, obs_kl_div, "Order composition", "KL divergence")


# 4. Clade composition
# --------------------

# Reference distribution
clades <- unique(S$Clade)
p <- numeric(length(clades))
for (i in 1:length(clades)) {
    p[i] <- (length(unique(S$cluster[S$Clade==clades[i]]))+1) / (length(unique(S$cluster))+length(clades))
}

# Observed kl_divergence
q <- numeric(length(clades))
for (i in 1:length(clades)) {
    q[i] <- (sum(D$OTUs[D$Clade==clades[i]])+1) / (sum(D$OTUs[D$Clade %in% clades])+length(clades))
}
obs_kl_div <- kldiv(p,q)

pred_kl_div <- numeric(ncol(P))
for (i in 1:ncol(P)) {
    otus <- P[,i]
    q <- numeric(length(clades))
    for (j in 1:length(clades)) {
        q[j] <- (sum(otus[D$Clade==clades[j]])+1) / (sum(otus[D$Clade %in% clades])+length(clades))
    }
    pred_kl_div[i] <- kldiv(p,q)
}
tail_prob <- sum(pred_kl_div <= obs_kl_div)

cat("Observed KL divergence in clade proportions is:", obs_kl_div, "\n")
cat("Tail probability (kl_div <= observed):", tail_prob, "\n")
cat("\n")

p4 <- plot_fxn(pred_kl_div, obs_kl_div, "Family-clade composition", "KL divergence")


# 5. Scarcity of parasitoids
# --------------------------

parasitoids <- c("Predator-parasitoid","Phytophage-parasitoid","Saprophage-parasitoid")
parasitoid_otus <- sum(D$OTUs[D$Niche %in% parasitoids])
total_otus <- sum(D$OTUs)
obs_frac <- parasitoid_otus / total_otus

pred_frac <- numeric(ncol(P))
for (i in 1:ncol(P)) {
    otus <- P[,i]
    pred_frac[i] <- sum(otus[D$Niche %in% parasitoids]) / sum(otus)
}
tail_prob <- sum(pred_frac<=obs_frac) / length(pred_frac)

cat("Observed fraction of parasitoids is:", obs_frac, "\n")
cat("Tail probability (frac <= observed):", tail_prob, "\n")
cat("\n")

p5 <- plot_fxn(pred_frac, obs_frac, "Parasitoid fraction", "Fraction")


# 6. Abundance of Coleoptera
# --------------------------

obs_frac <- sum(D$OTUs[D$Order=="Coleoptera"]) / sum(D$OTUs)

pred_frac <- numeric(ncol(P))
for (i in 1:ncol(P)) {
    otus <- P[,i]
    pred_frac[i] <- sum(otus[D$Order=="Coleoptera"]) / sum(otus)
}
tail_prob <- sum(pred_frac>=obs_frac) / length(pred_frac)

cat("Observed fraction of Coleoptera:", obs_frac, "\n")
cat("Tail probability (frac >= observed):", tail_prob, "\n")
cat("\n")

p6 <- plot_fxn(pred_frac, obs_frac, "Coleoptera fraction", "Fraction")

# Save plots
ggsave(file = "../figs/Fig_diversification_tests.jpg",
       width = 7.5,
       height = 10,
       plot = p1 + p2 + p3 + p4 + p5 + p6 +
           plot_layout(ncol=2) +
           plot_annotation(tag_levels="A")
       )

