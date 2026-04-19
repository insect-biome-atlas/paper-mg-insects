# Generate accumulation plots

library(vegan)
library(ggplot2)
library(patchwork)

# Set basic plot params
bs    <- 35
lsize <- 2

# Read data
# ---------

# Read OTU site occurrence tables
cat("Reading data\n")
site_otu_occurrence_combined_mg <- readRDS("../../iba_data/site_otu_occurrence_combined_mg.rds")
site_otu_occurrence_combined_se <- readRDS("../../iba_data/site_otu_occurrence_combined_se.rds")


# Generate total species accumulation curves
# ------------------------------------------

cat("Generating plot data\n")
spMG <- specaccum(site_otu_occurrence_combined_mg, method="random", permutations=100)
spSE <- specaccum(site_otu_occurrence_combined_se, method="random", permutations=100)
poolMG <- poolaccum(site_otu_occurrence_combined_mg, permutations=1000)
poolSE <- poolaccum(site_otu_occurrence_combined_se, permutations=1000)

cat("Formatting plot data\n")
mgDF <- data.frame(species=spMG$richness, sd=spMG$sd, sites=spMG$sites)
seDF <- data.frame(species=spSE$richness, sd=spSE$sd, sites=spSE$sites)
mgPoolDF <- data.frame(poolMG$means)
sePoolDF <- data.frame(poolSE$means)
mgPoolDF <- mgPoolDF[mgPoolDF$N >=5,]

# Madagascar
p1 <- ggplot(data=mgDF, aes(x=sites, y=species)) +
        geom_line(linewidth=lsize) +
        geom_ribbon(aes(ymin=species-sd, ymax=species+sd), alpha=0.2, fill="blue") +
        geom_line(data=mgPoolDF, linewidth=lsize, aes(x=N, y=Chao), linetype="dashed") +
        theme_linedraw(base_size=bs) +
        labs(x = "Number of sites" , y = "Number of OTUs")

# Sweden
p2 <- ggplot(data=seDF, aes(x=sites, y=species)) +
        geom_line(linewidth=lsize) +
        geom_ribbon(aes(ymin=species-sd, ymax=species+sd) , alpha=0.2, fill="blue") +
        geom_line(data=sePoolDF, linewidth=lsize, aes(x=N, y=Chao), linetype="dashed") +
        theme_linedraw(base_size=bs) +
        labs(x = "Number of sites" , y = "Number of OTUs")

# Save figure
ggsave( file = "../figs/Fig_accumulation_total.jpg",
        width = 28,
        height = 14,
        plot = p1 + p2 +
            plot_layout(axis_titles="collect", ncol=2) +
            plot_annotation(tag_levels="A")
        )

