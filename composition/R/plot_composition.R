# Plot site composition data in one
# joint figure

library(ggplot2)
library(patchwork)

# Read in colours for figures
# ---------------------------
source("../../fig_settings/fig_colours.R")
forest_cols <- forest_cols3

# Read in data
# ------------

otu_site_meta_mg <- read.delim("../data/otu_site_meta_mg.tsv")
otu_site_meta_se <- read.delim("../data/otu_site_meta_se.tsv")
source("filter_se_traps.R")
traps_to_keep <- filter_se_traps(20)    # Only keep se traps with 20 or more samples
otu_site_meta_se <- otu_site_meta_se[otu_site_meta_se$trapID %in% traps_to_keep,]
otu_site_meta_rf <- otu_site_meta_mg[otu_site_meta_mg$trap_habitat %in% c("Rainforest","Montane_Rainforest"),]
otu_site_meta_df <- otu_site_meta_mg[otu_site_meta_mg$trap_habitat == "Dry_Forest",]
otu_site_meta_rf$trap_habitat <- "Rainforest"
otu_site_meta_df$trap_habitat <- "Dry forest"
otu_site_meta_se$trap_habitat <- "Temperate forest"
otu_site_meta_mg$trapID <- paste0(otu_site_meta_mg$trapID,"_mg")
otu_site_meta_mg$trap_habitat <- "All MG forests"
D <- rbind(otu_site_meta_rf,otu_site_meta_df,otu_site_meta_mg,otu_site_meta_se)
D$trap_habitat <- factor(D$trap_habitat,levels=c("Rainforest","Dry forest","All MG forests","Temperate forest"))
idx <- which(D$Habitat=="Temporary habitats")
D$Habitat[idx] <- "Temporary"


# Define plot functions
# ---------------------

box_plot <- function(X, acc_vals, col, plot_title=NULL, y_label=NULL, y_limits=numeric()) {
    X <- X[X$trap_habitat!="All MG forests",]
    acc_vals <- acc_vals[acc_vals$trap_habitat!="All MG forests",]
    idx <- which(colnames(X)==col)
    colnames(X)[idx] <- "count"
    idx <- which(colnames(acc_vals)==col)
    colnames(acc_vals)[idx] <- "count"
    if (length(y_limits)==0)
        y_limits <- c(min(X$count),max(X$count))
    ggplot(X, aes(x=trap_habitat, y=count, fill=trap_habitat)) +
        theme_minimal() +
        geom_boxplot() +
        geom_point(data=acc_vals, aes(x=trap_habitat, y=count), shape=23, size=3, fill=diamond_col) +
        scale_x_discrete(name = NULL, labels = NULL) +  # drop x labels
        scale_fill_manual(values = forest_cols3,
                          labels = c("Madagascan rainforest", "Madagascan dry forest", "Temperate forest"), name = "") +
        ylab(y_label) +
        ylim(y_limits) +
        ggtitle(plot_title) +
        theme(
            legend.position = "bottom",
            legend.title = element_text(size = 8),
            legend.text = element_text(size = 7)
#           plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
#           axis.text.x = element_blank(),
#           axis.ticks.x = element_blank()
    )
}

box_plot_load <- function(D, plot_title = NULL) {
    D <- D[D$trap_habitat!="All MG forests",]
    ggplot(D, aes(x=trap_habitat, y=ph_ratio, fill=trap_habitat)) +
        theme_minimal() +
        geom_boxplot() +
        scale_x_discrete(name = NULL, labels = NULL) +  # drop x labels
        scale_fill_manual(values = forest_cols3,
                          labels = c("Madagascan rainforest","Madagascan dry forest","Temperate forest"), name = "") +
         ylab("Parasitoid load (read ratio)") +
        ggtitle(plot_title) +
        theme(
            legend.position = "bottom",
            legend.title = element_text(size = 8),
            legend.text = element_text(size = 7)
#           plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
#           axis.text.x = element_blank(),
#           axis.ticks.x = element_blank()
    )
}


# Compute and plot niche composition
# ----------------------------------

traps <- unique(D$trapID)
trap_habitat <- D$trap_habitat[match(traps,D$trapID)]
res <- data.frame(list(trapID=traps,trap_habitat=trap_habitat))
hosts <- c("Saprophage","Phytophage","Predator")
parasitoids <- c("Saprophage-parasitoid","Phytophage-parasitoid","Predator-parasitoid")
for (niche in c(hosts,parasitoids)) {
    X <- data.frame(table(D$trapID[D$Niche==niche]))
    colnames(X) <- c("trapID",niche)
    res <- merge(res,X)
}
for (i in 1:length(hosts)) {
    host <- hosts[i]
    parasitoid <- parasitoids[i]
    res[,parasitoid] <- res[,parasitoid] / res[,host]
}
sum_major <- rowSums(res[,hosts])
for (niche in hosts)
    res[,niche] <- res[,niche]/sum_major

# Compute accumulated values
A <- aggregate(D$read_count,by=list(trap_habitat=D$trap_habitat,cluster=D$cluster,Niche=D$Niche),FUN=sum)
acc_res <- data.frame(list(trap_habitat=unique(A$trap_habitat)))
for (niche in c(hosts,parasitoids)) {
    X <- data.frame(table(A$trap_habitat[A$Niche==niche]))
    colnames(X) <- c("trap_habitat",niche)
    acc_res <- merge(acc_res,X)
}
for (i in 1:length(hosts)) {
    host <- hosts[i]
    parasitoid <- parasitoids[i]
    acc_res[,parasitoid] <- acc_res[,parasitoid] / acc_res[,host]
}
sum_major <- rowSums(acc_res[,hosts])
for (niche in hosts)
    acc_res[,niche] <- acc_res[,niche]/sum_major

# Render plots
y_limits1 <- c(0.0,0.8)
y_limits2 <- c(0.0,1.05)
p1 <- box_plot(res, acc_res, "Saprophage", "(A)", "Proportion of major-niche OTUs", y_limits1)
p2 <- box_plot(res, acc_res, "Phytophage", "(B)", "Proportion of major-niche OTUs", y_limits1)
p3 <- box_plot(res, acc_res, "Predator", "(C)", "Proportion of major-niche OTUs", y_limits1)
p4 <- box_plot(res, acc_res, "Saprophage-parasitoid", "(D)", "Parasitoid:host OTU ratio", y_limits2)
p5 <- box_plot(res, acc_res, "Phytophage-parasitoid", "(E)", "Parasitoid:host OTU ratio", y_limits2)
p6 <- box_plot(res, acc_res, "Predator-parasitoid", "(F)", "Parasitoid:host OTU ratio", y_limits2)


# Compute and plot parasitoid load
# --------------------------------

# Aggregate reads for each trap and niche and trap habitat
X <- aggregate(read_count~trapID+Niche+trap_habitat,data=D,FUN=sum)

# Generate box plot data for saprophage communities
B1 <- X[X$Niche=="Saprophage",]
B2 <- X[X$Niche=="Saprophage-parasitoid",]
sapro <- data.frame(trapID=B1$trapID)
sapro$trap_habitat <- B1$trap_habitat[match(sapro$trapID,B1$trapID)]
sapro$trap_habitat <- factor(sapro$trap_habitat, levels=c("Rainforest","Dry forest","All MG forests","Temperate forest"))
sapro$h_reads <- B1$read_count[match(sapro$trapID,B1$trapID)]
sapro$p_reads <- B2$read_count[match(sapro$trapID,B2$trapID)]
sapro$p_reads[is.na(sapro$p_reads)] <- 0
sapro$ph_ratio <- sapro$p_reads/sapro$h_reads

# Generate box plot data for phytophage communities
B1 <- X[X$Niche=="Phytophage",]
B2 <- X[X$Niche=="Phytophage-parasitoid",]
phyto <- data.frame(trapID=B1$trapID)
phyto$trap_habitat <- B1$trap_habitat[match(phyto$trapID,B1$trapID)]
phyto$trap_habitat <- factor(phyto$trap_habitat, levels=c("Rainforest","Dry forest","All MG forests","Temperate forest"))
phyto$h_reads <- B1$read_count[match(phyto$trapID,B1$trapID)]
phyto$p_reads <- B2$read_count[match(phyto$trapID,B2$trapID)]
phyto$p_reads[is.na(phyto$p_reads)] <- 0
phyto$ph_ratio <- phyto$p_reads/phyto$h_reads

# Generate box plot data for predator communities
B1 <- X[X$Niche=="Predator",]
B2 <- X[X$Niche=="Predator-parasitoid",]
pred <- data.frame(trapID=B1$trapID)
pred$trap_habitat <- B1$trap_habitat[match(pred$trapID,B1$trapID)]
pred$trap_habitat <- factor(pred$trap_habitat, levels=c("Rainforest","Dry forest","All MG forests","Temperate forest"))
pred$h_reads <- B1$read_count[match(pred$trapID,B1$trapID)]
pred$p_reads <- B2$read_count[match(pred$trapID,B2$trapID)]
pred$p_reads[is.na(pred$p_reads)] <- 0
pred$ph_ratio <- pred$p_reads/pred$h_reads

# Render plots
p7 <- box_plot_load(sapro, "(G)")
p8 <- box_plot_load(phyto, "(H)")
p9 <- box_plot_load(pred, "(I)")


# Compute and plot habitat composition
# ------------------------------------

res <- data.frame(list(trapID=traps,trap_habitat=trap_habitat))
habitats <- c("Plants","Soil","Water","Wood","Temporary","Fungi")
for (habitat in habitats) {
    X <- data.frame(table(D$trapID[D$Habitat==habitat]))
    colnames(X) <- c("trapID",habitat)
    res <- merge(res,X)
}
sum_habitats <- rowSums(res[,habitats])
for (habitat in habitats)
    res[,habitat] <- res[,habitat]/sum_habitats

# Compute accumulated values
A <- aggregate(D$read_count,by=list(trap_habitat=D$trap_habitat,cluster=D$cluster,Habitat=D$Habitat),FUN=sum)
acc_res <- data.frame(list(trap_habitat=unique(A$trap_habitat)))
for (habitat in habitats) {
    X <- data.frame(table(A$trap_habitat[A$Habitat==habitat]))
    colnames(X) <- c("trap_habitat",habitat)
    acc_res <- merge(acc_res,X)
}
sum_habitats<- rowSums(acc_res[,habitats])
for (habitat in habitats)
    acc_res[,habitat] <- acc_res[,habitat]/sum_habitats

# Render plots
y_limits <- c(0.0,0.65)
p10 <- box_plot(res, acc_res, "Plants", "(J) Plants", "Proportion of OTUs", y_limits)
p11 <- box_plot(res, acc_res, "Soil", "(K) Soil", "Proportion of OTUs", y_limits)
p12 <- box_plot(res, acc_res, "Water", "(L) Water", "Proportion of OTUs", y_limits)
p13 <- box_plot(res, acc_res, "Wood", "(M) Wood", "Proportion of OTUs", y_limits)
p14 <- box_plot(res, acc_res, "Temporary", "(N) Temporary", "Proportion of OTUs", y_limits)
p15 <- box_plot(res, acc_res, "Fungi", "(O) Fungi", "Proportion of OTUs", y_limits)


# Compute and plot taxonomic composition
# --------------------------------------

res <- data.frame(list(trapID=traps,trap_habitat=trap_habitat))
big_five <- c("Diptera","Hymenoptera","Coleoptera","Lepidoptera","Hemiptera")
idx <- which(!(D$Order %in% big_five))
D$Order[idx] <- "Other"
orders <- c(big_five,"Other")
for (order in orders) {
    X <- data.frame(table(D$trapID[D$Order==order]))
    colnames(X) <- c("trapID",order)
    res <- merge(res,X)
}
sum_orders <- rowSums(res[,orders])
for (order in orders)
    res[,order] <- res[,order]/sum_orders

# Compute accumulated values
A <- aggregate(D$read_count,by=list(trap_habitat=D$trap_habitat,cluster=D$cluster,Order=D$Order),FUN=sum)
acc_res <- data.frame(list(trap_habitat=unique(A$trap_habitat)))
for (order in orders) {
    X <- data.frame(table(A$trap_habitat[A$Order==order]))
    colnames(X) <- c("trap_habitat",order)
    acc_res <- merge(acc_res,X)
}
sum_orders<- rowSums(acc_res[,orders])
for (order in orders)
    acc_res[,order] <- acc_res[,order]/sum_orders

# Render plots
y_limits <- c(0.0,0.65)
p16 <- box_plot(res, acc_res, "Diptera", "(P) Diptera", "Proportion of OTUs", y_limits)
p17 <- box_plot(res, acc_res, "Hymenoptera", "(Q) Hymenoptera", "Proportion of OTUs", y_limits)
p18 <- box_plot(res, acc_res, "Coleoptera", "(R) Coleoptera", "Proportion of OTUs", y_limits)
p19 <- box_plot(res, acc_res, "Lepidoptera", "(S) Lepidoptera", "Proportion of OTUs", y_limits)
p20 <- box_plot(res, acc_res, "Hemiptera", "(T) Hemiptera", "Proportion of OTUs", y_limits)
p21 <- box_plot(res, acc_res, "Other", "(U) Other", "Proportion of OTUs", y_limits)


# Put plots together into one composition figure
# ----------------------------------------------

plot1 <- p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9 +
            plot_layout(axis_titles="collect_y", guides="collect", ncol=3) &
            theme(legend.position="none")

plt <- p1 +
        theme(
            legend.position = "bottom",
            legend.title = element_text(size = 22),
            legend.text = element_text(size = 16))

plot2 <- cowplot::get_legend(plt)

p11 <- p11 + theme(axis.title.y=element_blank())
p12 <- p12 + theme(axis.title.y=element_blank())
p14 <- p14 + theme(axis.title.y=element_blank())
p15 <- p15 + theme(axis.title.y=element_blank())

plot3 <- p10 + p11 + p12 + p13 + p14 + p15 +
            plot_layout(guides="collect", ncol=3) &
            theme(legend.position="none")

p17 <- p17 + theme(axis.title.y=element_blank())
p18 <- p18 + theme(axis.title.y=element_blank())
p20 <- p20 + theme(axis.title.y=element_blank())
p21 <- p21 + theme(axis.title.y=element_blank())

plot4 <- p16 + p17 + p18 + p19 + p20 + p21 +
            plot_layout(guides="collect", ncol=3) &
            theme(legend.position="none")

# Save composite figure parts (make final edits in Illustrator)
ggsave(file = "../figs/Fig_composition_upperleft.pdf",
       width = 7.5,
       height = 9.0,
       plot = plot1
      )
ggsave(file = "../figs/Fig_composition_lowerleft.pdf",
       width = 7.5,
       height = 3.0,
       plot = plot2
      )
ggsave(file = "../figs/Fig_composition_upperright.pdf",
       width = 7.5,
       height = 6.0,
       plot = plot3
      )
ggsave(file = "../figs/Fig_composition_lowerright.pdf",
       width = 7.5,
       height = 6.0,
       plot = plot4
      )

