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
malaise_mg$trap_habitat <- meta_mg$trap_habitat[match(malaise_mg$sampleID_NGI,meta_mg$sampleID_NGI)]
malaise_se$trap_habitat <- meta_se$trap_habitat[match(malaise_se$sampleID_NGI,meta_se$sampleID_NGI)]
litter_mg$trap_habitat <- meta_mg$trap_habitat[match(litter_mg$sampleID_NGI,meta_mg$sampleID_NGI)]
litter_se$trap_habitat <- meta_se$trap_habitat[match(litter_se$sampleID_NGI,meta_se$sampleID_NGI)]
combined_mg$trap_habitat <- meta_mg$trap_habitat[match(combined_mg$sampleID_NGI,meta_mg$sampleID_NGI)]
combined_se$trap_habitat <- meta_se$trap_habitat[match(combined_se$sampleID_NGI,meta_se$sampleID_NGI)]
idx <- which(malaise_mg$trap_habitat=="Montane_Rainforest")
malaise_mg$trap_habitat[idx] <- "Rainforest"
idx <- which(litter_mg$trap_habitat=="Montane_Rainforest")
litter_mg$trap_habitat[idx] <- "Rainforest"
idx <- which(combined_mg$trap_habitat=="Montane_Rainforest")
combined_mg$trap_habitat[idx] <- "Rainforest"

# Get forest contributions per sample type
get_otus_by_habitat <- function(df) {
    x <- numeric()
    set1 <- unique(df$cluster[df$trap_habitat=="Rainforest"])
    set2 <- unique(df$cluster[df$trap_habitat=="Dry_Forest"])
    x[3] <- sum(set1 %in% set2)
    x[1] <- length(set1) - x[3]
    x[2] <- length(set2) - x[3]
    return (x)
}
D1 <- data.frame(list(habitat=c("Rainforest","Dry forest","Shared"),
                      OTUs=get_otus_by_habitat(malaise_mg),
                      dataset=rep("Malaise",times=3)
                      )
                )
D2 <- data.frame(list(habitat=c("Rainforest","Dry forest","Shared"),
                      OTUs=get_otus_by_habitat(litter_mg),
                      dataset=rep("Litter",times=3)
                      )
                )
D3 <- data.frame(list(habitat=c("Rainforest","Dry forest","Shared"),
                      OTUs=get_otus_by_habitat(combined_mg),
                      dataset=rep("Combined",times=3)
                      )
                )
D <- rbind(D1,D2)
D <- rbind(D,D3)
D$habitat <- factor(D$habitat,levels=c("Dry forest","Shared","Rainforest"))
D$dataset <- factor(D$dataset,levels=c("Combined","Litter","Malaise"))

# Get sample (dataset) composition by forest type
get_otus_by_dataset <- function(df_malaise, df_litter, habitat) {
    x <- numeric()
    set1 <- unique(df_malaise$cluster[df_malaise$trap_habitat==habitat])
    set2 <- unique(df_litter$cluster[df_litter$trap_habitat==habitat])
    x[3] <- sum(set1 %in% set2)
    x[1] <- length(set1) - x[3]
    x[2] <- length(set2) - x[3]
    return (x)
}
E1 <- data.frame(list(habitat="Rainforest",
                      OTUs=get_otus_by_dataset(malaise_mg, litter_mg, "Rainforest"),
                      dataset=c("Malaise","Litter","Shared")
                      )
                )
E2 <- data.frame(list(habitat="Dry forest",
                      OTUs=get_otus_by_dataset(malaise_mg, litter_mg, "Dry_Forest"),
                      dataset=c("Malaise","Litter","Shared")
                      )
                )
E3 <- data.frame(list(habitat="Temperate forest",
                      OTUs=get_otus_by_dataset(malaise_se, litter_se, "forest"),
                      dataset=c("Malaise","Litter","Shared")
                      )
                )
E <- rbind(E1,E2)
E <- rbind(E,E3)
E$dataset <- factor(E$dataset,levels=c("Litter","Shared","Malaise"))
E$habitat <- factor(E$habitat,levels=c("Temperate forest","Dry forest","Rainforest"))

# Generate plots
# --------------

# Render the plots
p1 <- ggplot(data=D, aes(x=dataset, y=OTUs, group=habitat, fill=habitat)) +
        geom_col(width=0.8) +
        coord_flip() +
        scale_fill_manual(values = c("Dry forest"="aquamarine","Shared"="green","Rainforest"="forestgreen"),
                          breaks = c("Rainforest","Shared","Dry forest"),
                          labels = c("Rainforest"="Rainforest","Shared"="Shared","Dry forest"="Dry forest")) +
        labs(title = NULL,
             x = "Sample",
             y = "Number of OTUs",
             fill = "Habitat") +
        theme_minimal(base_size=17)

p2 <- ggplot(data=E, aes(x=habitat, y=OTUs, group=dataset, fill=dataset)) +
        geom_col(width=0.8) +
        coord_flip() +
        scale_fill_manual(values = c("Litter"="lightgray","Shared"="blue","Malaise"="deepskyblue"), 
                          breaks = c("Malaise","Shared","Litter"),
                          labels = c("Malaise"="Malaise","Shared"="Shared","Litter"="Litter")) +
        labs(title = NULL,
             x = "Habitat",
             y = "Number of OTUs",
             fill = "Dataset") +
        theme_minimal(base_size=17)

# Put plots together
ggsave(
       file = "../figs/Fig_forest_sample_breakdown.jpg",
       width = 14.0,
       height = 5.0,
       plot = p1 + p2 + 
              plot_layout(ncol=2) +
              plot_annotation(tag_levels="A") & theme(legend.position="bottom")
       )

