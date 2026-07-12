# Plot validation of filtered data against known fauna
# and against sintax matching to BOLD


# Required packages: 
library(gridExtra)
library(ggplot2)
library(dplyr)
library(viridis)
library(RColorBrewer)
library(scales)
library(ggnewscale)
library(patchwork)


########### Fxns for plots ######################

# Dummy data for plot
dummy_data <- data.frame(
    x = c(1, 10),
    y1 = c(1, 10),
    y2 = c(0.5, 5)
)

# Match plot
plot_match <- function(df, lbl_pos) {

    ggplot(df, aes(x=Known_2017, y=OTUs)) +
        geom_point(alpha = 0.8, size = 3, aes(colour=Known_2017/Estimated)) +
        geom_text(data = subset(df, log(OTUs) > log(Known_2017) + 1.5),
                  aes(label = Family), color="black", alpha = 1.0, hjust = 0.3, vjust = -0.7) +
#       expand_limits(x = 0, y = 0) +
        scale_x_continuous(trans="log10",breaks=c(1,10,100,1000,10000,100000,1000000),
                           labels=c("1","10","100","1000","10000","100000","1000000")) +  
        scale_y_continuous(trans="log10",breaks=c(1,10,100,1000,10000,100000,1000000),
                           labels=c("1","10","100","1000","10000","100000","1000000")) +  
        scale_colour_viridis_c(option = "rocket", direction = -1, name = "Est. prop. known") +
        theme_minimal() +
        theme(legend.position = lbl_pos) +
        geom_abline(intercept = 0, slope = 1, linetype=1) +
#        geom_abline(intercept = log10(0.5), slope = 1, linetype=3) +
#        geom_line(data = dummy_data, aes(x = x, y = y1, linetype = "100% (y = x)")) +
#        geom_line(data = dummy_data, aes(x = x, y = y2, linetype = "50% (y = 0.5x)")) +
#        scale_linetype_manual(name = "Lines", values = c("100% (y = x)" = 1, "50% (y = 0.5x)" = 3)) +
#        guides(linetype = guide_legend(override.aes = list(colour = c("black", "black")))) +
        coord_cartesian(clip = "off") +
        ylab("No. OTUs found") +
        xlab("No. of known species") +
        theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), 
                           "inches"),
            axis.title.x = element_text(vjust=-1.0),
            axis.title.y = element_text(vjust=0.0),
            axis.text=element_text(size=16),
            axis.title=element_text(size=16),
            legend.text=element_text(size=12),
            legend.title=element_text(size=14),
            plot.tag = element_text(size=22)
        )
}


########### Read in data ####################
ref <- read.delim("../../traits/ronquist_2020_SE_traits_ncbi_taxonomy.tsv")
D1 <- aggregate(ref$Sweden.estimated.total, by=list(Family=ref$NCBI_Family), FUN=sum)
colnames(D1) <- c("Family","Estimated")
D2 <- aggregate(ref$Corrected.values.2017, by=list(Family=ref$NCBI_Family), FUN=sum)
colnames(D2) <- c("Family","Known_2017")
D <- merge(D1,D2)
D <- D[D$Known_2017>0,] # Only include those known in 2017

T <- readRDS("../../iba_data/cluster_taxonomy_se.rds")
TX <- read.delim("~/dev/figshare-repos/iba/processed_data/v4-prel/cleaned_noise_filtered_cluster_taxonomy_SE.tsv")
TX <- TX[TX$representative==1,]
C1 <- readRDS("../../iba_data/cluster_counts_malaise_long_se.rds")
C2 <- readRDS("../../iba_data/cluster_counts_litter_long_se.rds")
C <- rbind(C1,C2)
# C <- C[C$read_count > 10,]
T <- T[T$cluster %in% unique(C$cluster),]
CX1 <- readRDS("../../iba_data/raw_counts_malaise_se.rds")
CX2 <- readRDS("../../iba_data/raw_counts_litter_se.rds")
CX_clusters=unique(c(CX1$cluster,CX2$cluster))
TX <- TX[TX$cluster %in% CX_clusters,]

T1 <- T
# T1 <- T1[T1$aLWR > 0.8,]
E1 <- data.frame(table(T1$Family))
colnames(E1) <- c("Family","OTUs")
E2 <- data.frame(table(TX$Family))
colnames(E2) <- c("Family","OTUs")

TY <- read.delim("~/dev/figshare-repos/iba/processed_data/v4-prel/cluster_taxonomy_SE.tsv")
TY <- TY[TY$representative==1,]
E3 <- data.frame(table(TY$Family))
colnames(E3) <- c("Family","OTUs")

########### Make plots  #####################

X <- merge(D,E1)
idx <- which(X$Family=="Ectopsocidae")
X$Family[idx] <- ""     # Remove a colliding label
plot_C <- plot_match(X,"none")
plot_B <- plot_match(merge(D,E2),"right")
plot_A <- plot_match(merge(D,E3),"none")

########### Make final figure ###############
ggsave("../figs/Fig_validation.jpg",
       device="jpg",
       width=20.0,
       height=9.0,
       units="in",
#       dpi=300,
       plot = plot_B + plot_C + 
           plot_annotation(tag_levels="A")
        )

