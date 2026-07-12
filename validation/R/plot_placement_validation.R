# Plot match between sintax and epa-ng family assignments
# -------------------------------------------------------

# Dummy data for plot
dummy_data <- data.frame(
    x = c(1, 10),
    y1 = c(1, 10),
    y2 = c(0.5, 5)
)

# Match plot
plot_match <- function(df, lbl_pos) {

    ggplot(df, aes(x=sintax_otus, y=epang_otus)) +
        geom_point(alpha = 0.3, size = 3, colour="blue") +
#       expand_limits(x = 0, y = 0) +
        scale_x_continuous(trans="log10",breaks=c(1,10,100,1000,10000),
                           labels=c("1","10","100","1000","10000")) +  
        scale_y_continuous(trans="log10",breaks=c(1,10,100,1000,10000),
                           labels=c("1","10","100","1000","10000")) +  
        theme_minimal() +
        theme(legend.position = "none") +
        geom_abline(intercept = 0, slope = 1, linetype=1) +
#        geom_abline(intercept = log10(0.5), slope = 1, linetype=3) +
#        geom_line(data = dummy_data, aes(x = x, y = y1, linetype = "100% (y = x)")) +
#        geom_line(data = dummy_data, aes(x = x, y = y2, linetype = "50% (y = 0.5x)")) +
#        scale_linetype_manual(name = "Lines", values = c("100% (y = x)" = 1, "50% (y = 0.5x)" = 3)) +
#        guides(linetype = guide_legend(override.aes = list(colour = c("black", "black")))) +
        coord_cartesian(clip = "off") +
        ylab("Phylogenetic placement (# OTUs) ") +
        xlab("BOLD matching (# OTUs)") +
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
E1 <- data.frame(table(T1$Family))
colnames(E1) <- c("Family","epang_otus")
E2 <- data.frame(table(TX$Family))
colnames(E2) <- c("Family","sintax_otus")

D <- merge(E1,E2)

# Render plot
plot_A <- plot_match(D,"none")

# Save single plot
ggsave("../figs/Fig_sintax_vs_epang.jpg",
       width = 7.0,
       height = 7.0,
       plot = plot_A)

