# Script generating stacked bar charts showing
# catch by sample type and forest type

library(ggplot2)
library(patchwork)


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

# Add family and clade
tax_mg <- read.delim("../../iba_data/cluster_taxonomy_mg.tsv")
tax_se <- read.delim("../../iba_data/cluster_taxonomy_se.tsv")
combined_mg$Family <- tax_mg$Family[match(combined_mg$cluster,tax_mg$cluster)]
combined_se$Family <- tax_se$Family[match(combined_se$cluster,tax_se$cluster)]
combined_mg$Clade <- tax_mg$Clade[match(combined_mg$cluster,tax_mg$cluster)]
combined_se$Clade <- tax_se$Clade[match(combined_se$cluster,tax_se$cluster)]

mg_otus <- unique(combined_mg[,c("cluster","Family","Clade")])
se_otus <- unique(combined_se[combined_se$trap_habitat=="forest",c("cluster","Family","Clade")])
mg_rainforest_otus <- unique(combined_mg[combined_mg$trap_habitat=="Rainforest",c("cluster","Family","Clade")])
mg_dryforest_otus <- unique(combined_mg[combined_mg$trap_habitat=="Dry_Forest",c("cluster","Family","Clade")])

family_stats <- function(df) {
    res <- data.frame(table(df$Family))
    colnames(res) <- c("Family","OTUs")
    res <- res[order(res$OTUs,decreasing=TRUE),]
    return (res)
}
clade_stats <- function(df) {
    res <- data.frame(table(df$Clade))
    colnames(res) <- c("Family","OTUs")
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


# Generate plots
# --------------

# Plot rendering functions
plot_family <- function(D, num_top, plot_title) {
    D <- D[1:num_top,]
    D$Family <- factor(D$Family,levels=c(D$Family[num_top:1]))
    ggplot(data=D, aes(x=Family, y=OTUs)) +
        geom_col(width=0.8) +
        coord_flip() +
        labs(title = plot_title,
             x = NULL,
             y = "Number of OTUs") +
        theme_minimal(base_size=17)
}
plot_clade <- function(D, num_top, plot_title) {
    D <- D[1:num_top,]
    D$Clade <- factor(D$Clade,levels=c(D$Family[num_top:1]))
    ggplot(data=D, aes(x=Family, y=OTUs)) +
        geom_col(width=0.8) +
        coord_flip() +
        labs(title = plot_title,
             x = NULL,
             y = "Number of OTUs") +
        theme_minimal(base_size=17)
}

# Family composition
# ------------------

num_top <- 20
mg_families$OTUs <- -mg_families$OTUs
mg_rainforest_families$OTUs <- -mg_rainforest_families$OTUs
p1 <- plot_family(mg_families, num_top, "Tropical forest (Madagascar)")
p2 <- plot_family(se_families, num_top, "Temperate forest (Sweden)")
p3 <- plot_family(mg_rainforest_families, num_top, "Rainforest (Madagascar)")
p4 <- plot_family(mg_dryforest_families, num_top, "Dry forest (Madagascar)")

# Put plots together
ggsave(
       file = "../figs/Fig_family_composition.jpg",
       width = 14.0,
       height = 21.0,
       plot = p1 + p2 + p3 + p4 +
              plot_layout(ncol=2, axis_titles="collect") +
              plot_annotation(tag_levels="A") & theme(legend.position="bottom")
       )


# Clade composition
# -----------------

mg_clades$OTUs <- -mg_clades$OTUs
mg_rainforest_clades$OTUs <- -mg_rainforest_clades$OTUs
p1 <- plot_family(mg_clades, num_top, "Tropical forest (Madagascar)")
p2 <- plot_family(se_clades, num_top, "Temperate forest (Sweden)")
p3 <- plot_family(mg_rainforest_clades, num_top, "Rainforest (Madagascar)")
p4 <- plot_family(mg_dryforest_clades, num_top, "Dry forest (Madagascar)")

# Put plots together
ggsave(
       file = "../figs/Fig_clade_composition.jpg",
       width = 14.0,
       height = 21.0,
       plot = p1 + p2 + p3 + p4 +
              plot_layout(ncol=2, axis_titles="collect") +
              plot_annotation(tag_levels="A") & theme(legend.position="bottom")
       )

