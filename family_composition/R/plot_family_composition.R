# Script generating stacked bar charts showing
# catch by sample type and forest type

library(ggplot2)
library(patchwork)

source("plot_family_fxn.R")


# Get values we need
# ------------------

# Get counts data
malaise_mg <- readRDS("../../iba_data/cluster_counts_malaise_long_mg.rds")
malaise_se <- readRDS("../../iba_data/cluster_counts_malaise_long_se.rds")
litter_mg <- readRDS("../../iba_data/cluster_counts_litter_long_mg.rds")
litter_se <- readRDS("../../iba_data/cluster_counts_litter_long_se.rds")
combined_mg <- rbind(malaise_mg,litter_mg)
combined_se <- rbind(malaise_se,litter_se)

# Add trap habitat
meta_mg <- read.delim("../../iba_data/malaise_litter_sample_meta_mg.tsv")
meta_se <- read.delim("../../iba_data/malaise_litter_sample_meta_se.tsv")
combined_mg$trap_habitat <- meta_mg$trap_habitat[match(combined_mg$sampleID_NGI,meta_mg$sampleID_NGI)]
combined_se$trap_habitat <- meta_se$trap_habitat[match(combined_se$sampleID_NGI,meta_se$sampleID_NGI)]
idx <- which(combined_mg$trap_habitat=="Montane_Rainforest")
combined_mg$trap_habitat[idx] <- "Rainforest"

# Add order, family and clade
tax_mg <- readRDS("../../iba_data/cluster_taxonomy_mg.rds")
tax_se <- readRDS("../../iba_data/cluster_taxonomy_se.rds")
combined_mg$Order <- tax_mg$Order[match(combined_mg$cluster,tax_mg$cluster)]
combined_se$Order <- tax_se$Order[match(combined_se$cluster,tax_se$cluster)]
combined_mg$Family <- tax_mg$Family[match(combined_mg$cluster,tax_mg$cluster)]
combined_se$Family <- tax_se$Family[match(combined_se$cluster,tax_se$cluster)]
combined_mg$Clade <- tax_mg$Clade[match(combined_mg$cluster,tax_mg$cluster)]
combined_se$Clade <- tax_se$Clade[match(combined_se$cluster,tax_se$cluster)]

# Group minor orders
big_five <- c("Diptera","Hymenoptera","Coleoptera","Lepidoptera","Hemiptera")
idx1 <- which(!combined_mg$Order %in% big_five)
idx2 <- which(!combined_se$Order %in% big_five)
combined_mg$Order[idx1] <- "Other"
combined_se$Order[idx2] <- "Other"

mg_otus <- unique(combined_mg[,c("cluster","Order","Family","Clade")])
se_otus <- unique(combined_se[combined_se$trap_habitat=="forest",c("cluster","Order","Family","Clade")])
mg_rainforest_otus <- unique(combined_mg[combined_mg$trap_habitat=="Rainforest",c("cluster","Order","Family","Clade")])
mg_dryforest_otus <- unique(combined_mg[combined_mg$trap_habitat=="Dry_Forest",c("cluster","Order","Family","Clade")])

family_stats <- function(df) {
    res <- data.frame(table(df$Family))
    colnames(res) <- c("Family","OTUs")
    res$Order <- df$Order[match(res$Family,df$Family)]
    res$Order <- factor(res$Order, levels=c(big_five,"Other"))
    res <- res[order(res$OTUs,decreasing=TRUE),]
    return (res)
}
clade_stats <- function(df) {
    res <- data.frame(table(df$Clade))
    colnames(res) <- c("Clade","OTUs")
    res$Order <- df$Order[match(res$Clade,df$Clade)]
    res$Order <- factor(res$Order, levels=c(big_five,"Other"))
    res <- res[order(res$OTUs,decreasing=TRUE),]
    return (res)
}

mg_families <- family_stats(mg_otus)
mg_clades <- clade_stats(mg_otus)
se_families <- family_stats(se_otus)
se_clades <- clade_stats(se_otus)
mg_rainforest_families <- family_stats(mg_rainforest_otus)
mg_dryforest_families <- family_stats(mg_dryforest_otus)
mg_rainforest_clades <- clade_stats(mg_rainforest_otus)
mg_dryforest_clades <- clade_stats(mg_dryforest_otus)


# Plot family composition
# -----------------------

# Render plots
num_top <- 20
p1 <- plot_family(mg_families, num_top, 15000, 5, "Tropical forest (Madagascar)")
p2 <- plot_family(se_families, num_top, 3500, 7, "Temperate forest (Sweden)")

# Put plots together
ggsave(
       file = "../figs/Fig_family_composition.jpg",
       width = 14.0,
       height = 10.5,
       plot = p1 + p2  +
              plot_layout(ncol=2, axis_titles="collect") +
              plot_annotation(tag_levels="A") & theme(legend.position="bottom")
       )


# Plot clade composition
# ----------------------

# Render plots
p1 <- plot_clade(mg_clades, num_top, 10000, 5, "Tropical forest (Madagascar)")
p2 <- plot_clade(se_clades, num_top, 3000, 6, "Temperate forest (Sweden)")
p3 <- plot_clade(mg_rainforest_clades, num_top, 10000, 5, "Rainforest (Madagascar)")
p4 <- plot_clade(mg_dryforest_clades, num_top, 3500, 7, "Dry forest (Madagascar)")

# Put plots together
ggsave(
       file = "../figs/Fig_clade_composition.jpg",
       width = 14.0,
       height = 21.0,
       plot = p1 + p2 + p3 + p4 +
              plot_layout(ncol=2, axis_titles="collect") +
              plot_annotation(tag_levels="A") & theme(legend.position="bottom")
       )

