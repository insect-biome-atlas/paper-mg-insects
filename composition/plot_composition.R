# Script generating bar charts showing taxonomic and ecological composition
# of the MG and SE insect faunas

library(ggplot2)
library(patchwork)

# Read in data
# ------------

D1 <- read.delim("../iba_data/cluster_taxonomy_MG.tsv")
D2 <- read.delim("../iba_data/cluster_taxonomy_SE.tsv")

# Read in life history trait data
lht1 <- read.delim("../traits/clade_trait_data_mg.tsv") 
lht2 <- read.delim("../traits/clade_trait_data_se.tsv")
remove_annot <- function(D) {
    ranks <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
    idx <- which(colnames(D) %in% ranks)
    return (D[,-idx])
}  
lht1 <- remove_annot(lht1)
lht2 <- remove_annot(lht2)

# Merge data frames
D1 <- merge(D1,lht1,by="Clade")
D2 <- merge(D2,lht2,by="Clade")

# Write data frames
write.tsv <- function(D,file) { write.table(D,file,row.names=FALSE,sep="\t") }
write.tsv(D1,"total_mg.tsv")
write.tsv(D2,"total_se.tsv")

# Define a function adding a factor that groups all orders except the five largest
add_order_groups <- function(D) {
    big_five <- c("Hymenoptera", "Diptera", "Coleoptera", "Lepidoptera", "Hemiptera")
    x <- ""
    for ( i in 1:length(D$Order) ) {
        if ( D$Order[i] %in% big_five )
            x[i] <- as.character( D$Order[i] )
        else
            x[i] <- "Other"
    }
    D$Order_group <- factor(x)
    return (D)
}
D1 <- add_order_groups(D1)
D2 <- add_order_groups(D2)

# Compute subset data
M1 <- read.delim("../iba_data/cluster_counts_malaise_long_mg.tsv")
M2 <- read.delim("../iba_data/cluster_counts_malaise_long_se.tsv")
L1 <- read.delim("../iba_data/cluster_counts_litter_long_mg.tsv")
L2 <- read.delim("../iba_data/cluster_counts_litter_long_se.tsv")

malaise_samples_mg <- read.delim("../iba_data/malaise_sample_meta_mg.tsv")
malaise_samples_se <- read.delim("../iba_data/malaise_sample_meta_se.tsv")
litter_samples_mg <- read.delim("../iba_data/litter_sample_meta_mg.tsv")
litter_samples_se <- read.delim("../iba_data/litter_sample_meta_se.tsv")

# Habitat cluster sets for subsetting
se_forest_mt <- malaise_samples_se$sampleID_NGI[malaise_samples_se$trap_habitat=="forest"]
se_forest_lit <- litter_samples_se$sampleID_NGI[litter_samples_se$trap_habitat=="forest"]

X <- c(M2$cluster[M2$sampleID_NGI %in% se_forest_mt],L2$cluster[L2$sampleID_NGI %in% se_forest_lit])
D3 <- D2[D2$cluster %in% X,]


# Taxonomic composition
# ---------------------

composition_df <- function(D, col, name=col) {
    x <- data.frame(table(D[,col]))
    colnames(x) <- c(name,"OTUs")
    x$Fraction <- x$OTUs/sum(x$OTUs)
    return (x)
}
tax_comp_mg <- composition_df(D1,"Order_group","Order")
tax_comp_se <- composition_df(D2,"Order_group","Order")
# write.tsv(tax_comp_mg,"taxonomic_composition_mg.tsv")
# write.tsv(tax_comp_se,"taxonomic_composition_se.tsv")
tax_comp_se_f <- composition_df(D3,"Order_group","Order")

# Niche composition
# -----------------

niche_comp_mg <- composition_df(D1,"Niche")
niche_comp_se <- composition_df(D2,"Niche")
# write.tsv(niche_comp_mg,"niche_composition_mg.tsv")
# write.tsv(niche_comp_se,"niche_composition_se.tsv")
niche_comp_se_f <- composition_df(D3,"Niche")

# Habitat composition
# -------------------

habitat_comp_mg <- composition_df(D1,"Habitat")
habitat_comp_se <- composition_df(D2,"Habitat")
# write.tsv(habitat_comp_mg,"habitat_composition_mg.tsv")
# write.tsv(habitat_comp_se,"habitat_composition_se.tsv")
habitat_comp_se_f <- composition_df(D3,"Habitat")


# Horizontal bar charts
plot_Tax <- function(df) {
    p <- ggplot(df, aes(x = Order, y = OTUs)) +
        geom_col() +
        coord_flip() +
        theme_minimal()
    return (p)
}

plot_Niche <- function(df) {
    p <- ggplot(df, aes(x = Niche, y = OTUs)) +
        geom_col() +
        coord_flip() +
        theme_minimal()
    return (p)
}

plot_Habitat <- function(df) {
    p <- ggplot(df, aes(x = Habitat, y = OTUs)) +
        geom_col() +
        coord_flip() +
        theme_minimal()
    return (p)
}

p1 <- plot_Tax(tax_comp_mg)
p2 <- plot_Tax(tax_comp_se)
p2A <- plot_Tax(tax_comp_se_f)

p3 <- plot_Niche(niche_comp_mg)
p4 <- plot_Niche(niche_comp_se)
p4A <- plot_Niche(niche_comp_se_f)

p5 <- plot_Habitat(habitat_comp_mg)
p6 <- plot_Habitat(habitat_comp_se)
p6A <- plot_Habitat(habitat_comp_se_f)

# Put plots together into figures

# Plot total data
total_plot <- (p1 + p2 + p3 + p4 + p5 + p6) +
  plot_layout(ncol = 2) +
  plot_annotation(tag_levels = "A")

ggsave(file = "Fig_composition_total.jpg",
       plot = total_plot,  height = 10.5,  width = 7,  dpi = 300)

# Focus on forest data
forest_plot <- (p1 + p2A + p3 + p4A + p5 + p6A) +
  plot_layout(ncol = 2) +
  plot_annotation(tag_levels = "A")

ggsave(file = "Fig_composition_forest.jpg",
       plot = total_plot,  height = 10.5,  width = 7,  dpi = 300)

# Compare dry forest and rainforest in MG


