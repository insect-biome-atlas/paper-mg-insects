# Compute and plot predictive diversification tests
# Use predictive samples corrected for the Coleoptera sampling bias

library(ggplot2)
library(patchwork)
source("../../fig_settings/fig_colours.R")


# Get MG data, and samples from the predictive distribution
D <- read.delim("../data/diversification_data.tsv")
P <- read.delim("../data/predictive_samples_bias.tsv")

# Get SE data and complement with clade info
S <- read.delim("../../composition/data/otu_site_meta_se.tsv")
T <- readRDS("../../iba_data/cluster_taxonomy_SE.rds")
S$Clade <- T$Clade[match(S$cluster,T$cluster)]

# KL divergence
kldiv <- function(p,q) {
    sum(p*log(p/q))
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

# Accumulate results
lbl <- "Major-niche composition"
max_scale_val <- max(pred_kl_div,obs_kl_div)
res1 <- data.frame(label=lbl,vals=pred_kl_div/max_scale_val)
res2 <- data.frame(label=lbl,obs_val=obs_kl_div/max_scale_val,ref_val=0.0)


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

# Accumulate results
lbl <- "Microhabitat use"
max_scale_val <- max(pred_kl_div,obs_kl_div)
res1 <- rbind(res1, data.frame(label=lbl,vals=pred_kl_div/max_scale_val))
res2 <- rbind(res2, data.frame(label=lbl,obs_val=obs_kl_div/max_scale_val,ref_val=0.0))


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

# Accumulate results
lbl <- "Order composition"
max_scale_val <- max(pred_kl_div,obs_kl_div)
res1 <- rbind(res1, data.frame(label=lbl,vals=pred_kl_div/max_scale_val))
res2 <- rbind(res2, data.frame(label=lbl,obs_val=obs_kl_div/max_scale_val,ref_val=0.0))

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

# Accumulate results
lbl <- "Family-clade composition"
max_scale_val <- max(pred_kl_div,obs_kl_div)
res1 <- rbind(res1, data.frame(label=lbl,vals=pred_kl_div/max_scale_val))
res2 <- rbind(res2, data.frame(label=lbl,obs_val=obs_kl_div/max_scale_val,ref_val=0.0))

# Prepare for plotting
res1$label <- factor(res1$label,levels=res2$label[4:1])


# Generate first plot (aspects 1-4)
# ---------------------------------

p1 <- ggplot(data=res1, aes(x=label, y=vals)) +
        geom_violin(width=1.0, linewidth=0.2, fill="steelblue", alpha=0.5) +
        geom_point(data=res2, aes(x=label,y=obs_val), shape=23, size=3, fill=mg_col) +
        geom_point(data=res2, aes(x=label,y=ref_val), shape=23, size=3, fill=se_col) +
#        theme_ipsum() +
        theme_minimal() +
        theme(
            legend.position="none"
        ) +
        coord_flip() + # This switch X and Y axis and allows to get the horizontal version
        xlab("") +
        ylab("Distance (Kullback-Leibler, scaled)")


# 5. Scarcity of parasitoids
# --------------------------

parasitoids <- c("Predator-parasitoid","Phytophage-parasitoid","Saprophage-parasitoid")
all_niches <- c(niches, parasitoids)
parasitoid_otus <- sum(D$OTUs[D$Niche %in% parasitoids])
total_otus <- sum(D$OTUs[D$Niche %in% all_niches])
obs_frac <- parasitoid_otus / total_otus

# SE reference value
parasitoid_otus <- length(unique(S$cluster[S$Niche %in% parasitoids]))
total_otus <- length(unique(S$cluster[S$Niche %in% all_niches]))
se_frac <- parasitoid_otus / total_otus

pred_frac <- numeric(ncol(P))
for (i in 1:ncol(P)) {
    otus <- P[,i]
    pred_frac[i] <- sum(otus[D$Niche %in% parasitoids]) / sum(otus)
}
tail_prob <- sum(pred_frac>=obs_frac) / length(pred_frac)

cat("Observed fraction of parasitoids is:", obs_frac, "\n")
cat("Tail probability (frac >= observed):", tail_prob, "\n")
cat("\n")

# Accumulate results
lbl <- "Parasitoid fraction"
res3 <- data.frame(label=lbl,vals=pred_frac)
res4 <- data.frame(label=lbl,obs_val=obs_frac,ref_val=se_frac)


# 6. Abundance of Coleoptera among non-parasitoids
# ------------------------------------------------

# Observed MG value
D1 <- D[!(D$Niche %in% parasitoids),]
obs_frac <- sum(D1$OTUs[D1$Order=="Coleoptera"]) / sum(D1$OTUs)

# SE reference value
S1 <- S[!(S$Niche %in% parasitoids),]
se_frac <- length(unique(S1$cluster[S1$Order=="Coleoptera"])) / length(unique(S1$cluster))

pred_frac <- numeric(ncol(P))
for (i in 1:ncol(P)) {
    otus <- P[,i]
    otus1 <- otus[!(D$Niche %in% parasitoids)]                          # Subset otus to non-parasitoids
    pred_frac[i] <- sum(otus1[D1$Order=="Coleoptera"]) / sum(otus1)     # Use 1 to 1 mapping between otus1 and D1
}
tail_prob <- sum(pred_frac<=obs_frac) / length(pred_frac)

cat("Observed fraction of Coleoptera among non-parasitoids:", obs_frac, "\n")
cat("Tail probability (frac <= observed):", tail_prob, "\n")
cat("\n")

# Accumulate results
lbl <- "Coleoptera fraction"
res3 <- rbind(res3,data.frame(label=lbl,vals=pred_frac))
res4 <- rbind(res4,data.frame(label=lbl,obs_val=obs_frac,ref_val=se_frac))

# Prepare for plotting
res3$label <- factor(res3$label,levels=res4$label[2:1])


# Generate plots
# --------------

p1 <- ggplot(data=res1, aes(x=label, y=vals)) +
        geom_violin(width=0.95, linewidth=0.2, fill="steelblue", alpha=0.5) +
        geom_point(data=res2, aes(x=label,y=obs_val), shape=23, size=3, fill=mg_col) +
        geom_point(data=res2, aes(x=label,y=ref_val), shape=23, size=3, fill=se_col) +
        theme_linedraw() +
        theme(
            legend.position="none",
            axis.title.x=element_text(margin=margin(t=5)),
            axis.text.y=element_text(size=12)
            ) +
        coord_flip() + # This switch X and Y axis and allows to get the horizontal version
        xlab("") +
        ylab("Distance (Kullback-Leibler, scaled)")

p2 <- ggplot(data=res3, aes(x=label, y=vals)) +
        geom_violin(width=0.95, linewidth=0.2, fill="steelblue", alpha=0.5) +
        geom_point(data=res4, aes(x=label,y=obs_val), shape=23, size=3, fill=mg_col) +
        geom_point(data=res4, aes(x=label,y=ref_val), shape=23, size=3, fill=se_col) +
        theme_linedraw() +
        theme(
            legend.position="none",
            axis.title.x=element_text(margin=margin(t=5)),
            axis.text.y=element_text(size=12)
            ) +
        coord_flip() + # This switch X and Y axis and allows to get the horizontal version
        xlab("") +
        ylab("Proportion of species")

# Save plots
ggsave(file = "../figs/Fig_diversification_tests_bias_violin.jpg",
       width = 7.5,
       height = 9.0,
       plot = p1 + p2 +
           plot_layout(ncol=1, heights=c(4,2)) +
           plot_annotation(tag_levels="A")
       )

