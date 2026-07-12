# Plot shallow-sample family stats

source("plot_family_fxn.R")


# Compute composition stats by sample
# -----------------------------------

# Prepare data
D1 <- data.frame(readRDS("../../iba_data/cluster_counts_malaise_long_mg.rds"))
D2 <- data.frame(readRDS("../../iba_data/cluster_counts_malaise_long_se.rds"))
T1 <- readRDS("../../iba_data/cluster_taxonomy_mg.rds")
T2 <- readRDS("../../iba_data/cluster_taxonomy_se.rds")

D1$Order <- T1$Order[match(D1$cluster,T1$cluster)]
D1$Family <- T1$Family[match(D1$cluster,T1$cluster)]
D1$Clade <- T1$Clade[match(D1$cluster,T1$cluster)]
D1$OTUs <- 1
big_five <- c("Diptera","Hymenoptera","Coleoptera","Lepidoptera","Hemiptera")
idx <- which(!(D1$Order %in% big_five))
D1$Order[idx] <- "Other"
D1$Order <- factor(D1$Order, levels=c(big_five,"Other"))

D2$Order <- T2$Order[match(D2$cluster,T2$cluster)]
D2$Family <- T2$Family[match(D2$cluster,T2$cluster)]
D2$Clade <- T2$Clade[match(D2$cluster,T2$cluster)]
D2$OTUs <- 1
idx <- which(!(D2$Order %in% big_five))
D2$Order[idx] <- "Other"
D2$Order <- factor(D2$Order, levels=c(big_five,"Other"))

# A function for printing stats if we are interested
print_stats <- function(X, title) {
    D <- colSums(X)
    cat(title,"\n")
    print(head(D[order(D,decreasing=TRUE)], n=20))
}

# Get the data we need
X <- xtabs(D1$OTUs~D1$sampleID_NGI+D1$Family)
mg_families <- data.frame(colSums(X))
colnames(mg_families) <- "OTUs"
mg_families$Family <- row.names(mg_families)
mg_families$Order <- D1$Order[match(mg_families$Family,D1$Family)]
mg_families <- mg_families[order(mg_families$OTUs,decreasing=TRUE),]
mg_families$OTUs <- 100*mg_families$OTUs/max(mg_families$OTUs)

X <- xtabs(D2$OTUs~D2$sampleID_NGI+D2$Family)
se_families <- data.frame(colSums(X))
colnames(se_families) <- "OTUs"
se_families$Family <- row.names(se_families)
se_families$Order <- D2$Order[match(se_families$Family,D2$Family)]
se_families <- se_families[order(se_families$OTUs,decreasing=TRUE),]
se_families$OTUs <- 100*se_families$OTUs/max(se_families$OTUs)

# Render plots
num_top <- 20
p1 <- plot_family(mg_families, num_top, 100, 5, "Tropical forest (Madagascar)")
p2 <- plot_family(se_families, num_top, 100, 5, "Temperate forest (Sweden)")

# Put plots together
ggsave(
       file = "../figs/Fig_shallow_family_composition.jpg",
       width = 14.0,
       height = 10.5,
       plot = p1 + p2  +
              plot_layout(ncol=2, axis_titles="collect", guides="collect") +
              plot_annotation(tag_levels="A") & theme(legend.position="bottom")
       )

