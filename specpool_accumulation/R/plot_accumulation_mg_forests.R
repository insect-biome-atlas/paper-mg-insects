# Generate accumulation plots

library(vegan)
library(ggplot2)
library(patchwork)

# Set basic plot params
bs    <- 20
lsize <- 2

# Read data
# ---------

cat("Reading data\n")

# Read species * site matrix
sp_matrix_mg <- readRDS("../../iba_data/site_otu_occurrence_combined_mg.rds")

# Read in sample metadata
site_meta_mg <- read.delim("../../iba_data/malaise_litter_sample_meta_mg.tsv")

# Generate trap habitat subsets
site_meta_df <- site_meta_mg[site_meta_mg$trap_habitat=="Dry_Forest",]
site_meta_rf <- site_meta_mg[site_meta_mg$trap_habitat %in% c("Montane_Rainforest","Rainforest"),]

sp_matrix_df <- sp_matrix_mg[rownames(sp_matrix_mg) %in% unique(site_meta_df$trapID),]
sp_matrix_rf <- sp_matrix_mg[rownames(sp_matrix_mg) %in% unique(site_meta_rf$trapID),]


# Generate total species accumulation curves
# ------------------------------------------

cat("Generating plot data\n")

spDF <- specaccum(sp_matrix_df, method="random", permutations=100)
spRF <- specaccum(sp_matrix_rf, method="random", permutations=100)
poolDF <- poolaccum(sp_matrix_df, permutations=1000)
poolRF <- poolaccum(sp_matrix_rf, permutations=1000)

cat("Formatting plot data\n")

dfDF <- data.frame(species=spDF$richness, sd=spDF$sd, sites=spDF$sites)
rfDF <- data.frame(species=spRF$richness, sd=spRF$sd, sites=spRF$sites)
dfPoolDF <- data.frame(poolDF$means)
rfPoolDF <- data.frame(poolRF$means)

# Rain forest
p1 <- ggplot(data=rfDF, aes(x=sites, y=species)) +
        ggtitle("Madagascar (all forests)") +
        geom_line(linewidth=lsize) +
        geom_ribbon(aes(ymin=species-sd, ymax=species+sd), alpha=0.2, fill="blue") +
        geom_line(data=rfPoolDF, linewidth=lsize, aes(x=N, y=Chao), linetype="dashed") +
        theme_linedraw(base_size=bs) +
        labs(x = "Number of sites" , y = "Number of OTUs")

# Sweden
p2 <- ggplot(data=dfDF, aes(x=sites, y=species)) +
        gg_title("Temperate forest (Sweden)") + 
        geom_line(linewidth=lsize) +
        geom_ribbon(aes(ymin=species-sd, ymax=species+sd) , alpha=0.2, fill="blue") +
        geom_line(data=dfPoolDF, linewidth=lsize, aes(x=N, y=Chao), linetype="dashed") +
        theme_linedraw(base_size=bs) +
        labs(x = "Number of sites" , y = "Number of OTUs")

# Save figure
ggsave( file = "../figs/Fig_accumulation_mg_forests.jpg",
        width = 14,
        height = 7,
        plot = p1 + p2 +
            plot_layout(axis_titles="collect", ncol=2) +
            plot_annotation(tag_levels="A")
        )

