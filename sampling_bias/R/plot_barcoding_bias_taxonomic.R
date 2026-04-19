# Plot barcoding bias

library(ggplot2)
library(patchwork)

# Read in barcoding data (15 samples)
D <- read.delim("../Barcoding_cleaned_matched_corrected.csv",sep=";")

# Restrict to hexapods
hexapods <- c("Insecta","Protura","Diplura","Collembola")
D <- D[D$Class %in% hexapods,]

# Add life history annotations
# The family->clade matching will create problems for a few specimens,
# but the annotations for them are too unprecise to allow Clade matching
# anyway, so it makes sense to disregard them
lht <- read.delim("../../traits/clade_trait_data_se.tsv")
D$Niche <- lht$Niche[match(D$Family,lht$Clade)]
D$Habitat <- lht$Habitat[match(D$Family,lht$Clade)]

# Compute the success rates
big_five <- c("Diptera","Hymenoptera","Coleoptera","Lepidoptera","Hemiptera")
D$Order_group <- D$Order
D$Order_group[!(D$Order_group %in% big_five)] <- "Other"
D$matched <- grepl("YES",D$Matched_to_Cluster)

generate_plot_data <- function(D, col, new_name) {
    X <- data.frame(table(D[,col]))
    colnames(X) <- c(new_name,"total_specimens")
    Y <- data.frame(table(D[D$matched,col]))
    colnames(Y) <- c(new_name,"matched_specimens")
    E <- merge(X,Y)
    E$success_rate <- E$matched_specimens / E$total_specimens
    return(E)
}

E1 <- generate_plot_data(D, "Order_group", "Order")
E1$Order <- factor(E1$Order,levels=c("Other","Hemiptera","Lepidoptera","Coleoptera","Hymenoptera","Diptera"))
E2 <- generate_plot_data(D, "Niche", "Niche")
E2$Niche <- factor(E2$Niche,levels=c("Predator-parasitoid","Predator","Saprophage-parasitoid","Saprophage","Phytophage-parasitoid","Phytophage"))
E3 <- generate_plot_data(D, "Habitat", "Habitat")
E3$Habitat <- factor(E3$Habitat,levels=c("Fungi","Temporary habitats","Wood","Water","Soil","Plants"))

p1 <- ggplot(E1, aes(x=Order, y=success_rate)) +
        geom_col(width=0.7) +
        coord_flip() +
        xlab("Order") +
        ylab("Success rate") +
        theme_minimal(base_size=15)

p2 <- ggplot(E2, aes(x=Niche, y=success_rate)) +
        geom_col(width=0.7) +
        coord_flip() +
        xlab("Niche") +
        ylab("Success rate") +
        theme_minimal(base_size=15)

p3 <- ggplot(E3, aes(x=Habitat, y=success_rate)) +
        geom_col(width=0.7) +
        coord_flip() +
        xlab("Microhabitat") +
        ylab("Success rate") +
        theme_minimal(base_size=15)

ggsave( 
        file = "../figs/Fig_barcoding_bias.jpg",
        width = 21,
        height = 7,
        plot = p1 + p2 + p3 +
            plot_layout(axis_titles="collect", ncol=3) +
            plot_annotation(tag_levels="A")
      )

