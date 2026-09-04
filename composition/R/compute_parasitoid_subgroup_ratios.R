# Compute species ratios for various subgroups of
# parasitoids for testing nasty host hypothesis

library(ggplot2)
library(patchwork)

# Read in data
# ------------
otu_site_meta_mg <- read.delim("../data/otu_site_meta_mg.tsv")
otu_site_meta_se <- read.delim("../data/otu_site_meta_se.tsv")
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

T1 <- readRDS("../../iba_data/cluster_taxonomy_mg.rds")
T2 <- readRDS("../../iba_data/cluster_taxonomy_se.rds")
D1 <- otu_site_meta_mg
D2 <- otu_site_meta_se
D1$Clade <- T1$Clade[match(D1$cluster,T1$cluster)]
D2$Clade <- T2$Clade[match(D2$cluster,T2$cluster)]


# Define function computing proportions
# -------------------------------------
mg2se_clades <- function(x, D1, D2) { length(unique(D1$cluster[D1$Clade%in%x])) / length(unique(D2$cluster[D2$Clade%in%x])) }
mg2se <- function(D1, D2) { length(unique(D1$cluster)) / length(unique(D2$cluster)) }


# Test various hypotheses
# -----------------------

# Overall ratio
f <- mg2se(D1,D2)
cat("Overall ratio of forest species diversity is:",f,"\n")

# Ratio of phytophage parasitoids in wood
f <- mg2se(D1[D1$Niche=="Phytophage-parasitoid" & D1$Habitat=="Wood",],D2[D2$Niche=="Phytophage-parasitoid" & D2$Habitat=="Wood",])
cat("Ratio of phytophage parasitoids in wood is:",f,"\n")

# Ratio of Scelionidae
f <- mg2se_clades("Scelionidae",D1,D2)
cat("Ratio of Scelionidae is:",f,"\n")

# Ratio of Mymaridae
f <- mg2se_clades("Mymaridae",D1,D2)
cat("Ratio of Mymaridae is:",f,"\n")

# Cryptinae
f <- mg2se_clades("Cryptinae",D1,D2)
cat("Ratio of Cryptinae is:",f,"\n")

# Pimplinae (without Polysphincta group)
f <- mg2se_clades("Pimplinae_s_str",D1,D2)
cat("Ratio of Pimplinae (excluding Polysphincta group) is:",f,"\n")

# Phygadeuontinae including Chirotica group
f <- mg2se_clades(c("Phygadeuontinae_s_str","Chirotica_group"),D1,D2)
cat("Ratio of Phygadeuontinae (including Chirotica group) is:",f,"\n")

