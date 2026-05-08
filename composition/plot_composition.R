# Script generating bar charts showing taxonomic and ecological composition
# of the MG and SE (forest) insect faunas

library(ggplot2)
library(patchwork)

# Read in and prepare data
# ------------------------

# Read in taxonomy data
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

# Only keep clusters present in forest sites in Sweden
se_forest_mt <- malaise_samples_se$sampleID_NGI[malaise_samples_se$trap_habitat=="forest"]
se_forest_lit <- litter_samples_se$sampleID_NGI[litter_samples_se$trap_habitat=="forest"]
X <- c(M2$cluster[M2$sampleID_NGI %in% se_forest_mt],L2$cluster[L2$sampleID_NGI %in% se_forest_lit])
D2 <- D2[D2$cluster %in% unique(X),]


# Function computing composition
# ------------------------------

composition_df <- function(D, col, name=col) {
    x <- data.frame(table(D[,col]))
    colnames(x) <- c(name,"OTUs")
    x$Fraction <- x$OTUs/sum(x$OTUs)
    x$kOTUs <- x$OTUs/1000
    return (x)
}


# Taxonomic composition
# ---------------------

orders <- c("Diptera","Hymenoptera","Coleoptera","Lepidoptera","Hemiptera","Other")
orders <- orders[6:1]
tax_comp_mg <- composition_df(D1,"Order_group","Order")
tax_comp_se <- composition_df(D2,"Order_group","Order")
tax_comp_mg$Order <- factor(tax_comp_mg$Order,levels=orders)
tax_comp_se$Order <- factor(tax_comp_se$Order,levels=orders)

# Niche composition
# -----------------

niches <- c("Saprophage","Saprophage-parasitoid","Phytophage","Phytophage-parasitoid","Predator","Predator-parasitoid")
niches <- niches[6:1]
niche_comp_mg <- composition_df(D1,"Niche")
niche_comp_se <- composition_df(D2,"Niche")
niche_comp_mg <- niche_comp_mg[niche_comp_mg$Niche %in% niches,]
niche_comp_se <- niche_comp_se[niche_comp_se$Niche %in% niches,]
niche_comp_mg$Niche <- factor(niche_comp_mg$Niche,levels=niches)
niche_comp_se$Niche <- factor(niche_comp_se$Niche,levels=niches)

# Habitat composition
# -------------------

habitats <- c("Plants","Soil","Water","Wood","Temporary habitats","Fungi")
habitats <- habitats[6:1]
habitat_comp_mg <- composition_df(D1,"Habitat")
habitat_comp_se <- composition_df(D2,"Habitat")
habitat_comp_mg <- habitat_comp_mg[habitat_comp_mg$Habitat %in% habitats,]
habitat_comp_se <- habitat_comp_se[habitat_comp_se$Habitat %in% habitats,]
habitat_comp_mg$Habitat <- factor(habitat_comp_mg$Habitat,levels=habitats)
habitat_comp_se$Habitat <- factor(habitat_comp_se$Habitat,levels=habitats)


# Generate plots
# --------------

# Horizontal bar chart fxns
plot_Tax <- function(df) {
    ggplot(df, aes(x = Order, y = kOTUs)) +
        geom_col() +
        coord_flip() +
        theme_minimal()
}

plot_Niche <- function(df) {
    ggplot(df, aes(x = Niche, y = kOTUs)) +
        geom_col() +
        coord_flip() +
        theme_minimal()
}

plot_Habitat <- function(df) {
    ggplot(df, aes(x = Habitat, y = kOTUs)) +
        geom_col() +
        coord_flip() +
        theme_minimal()
}

# Render plots
tax_comp_mg$kOTUs <- -tax_comp_mg$kOTUs
p1 <- plot_Tax(tax_comp_mg)
p2 <- plot_Tax(tax_comp_se)

niche_comp_mg$kOTUs <- -niche_comp_mg$kOTUs
p3 <- plot_Niche(niche_comp_mg)
p4 <- plot_Niche(niche_comp_se)

habitat_comp_mg$kOTUs <- -habitat_comp_mg$kOTUs
p5 <- plot_Habitat(habitat_comp_mg)
p6 <- plot_Habitat(habitat_comp_se)

# Put plots together
ggsave(file = "Fig_composition.jpg",
       width = 7.0,
       height = 10.5,
       plot = p1 + p2 + p3 + p4 + p5 + p6 +
              plot_layout(ncol=2, axis_titles="collect") +
              plot_annotation(tag_levels="A")
       )

