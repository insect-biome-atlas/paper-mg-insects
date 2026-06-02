# Compute reference values for estimated Swedish insect fauna

# Read in data on the estimate
D <- read.delim("../../traits/ronquist_2020_SE_traits.csv",sep=";")

# Remove empty rows
D <- D[!is.na(D$Sorting.Number),]

# Correct coding
D$Main.feeding.niche[D$Main.feeding.niche=="Saprophagous"] <- "Saprophage"
D$Main.feeding.niche[D$Main.feeding.niche=="Phytophagous"] <- "Phytophage"
D$Main.feeding.habitat[D$Main.feeding.habitat=="Temporary Habitats"]<-"Temporary habitats"

# Add and correct data
D$Family <- D$Taxon..Dyntaxa.2017.
D$Main.feeding.niche[D$Family=="Mymaridae"] <- "Phytophage-parasitoid"      # Majority coding
D$Main.feeding.niche[D$Family=="Pteromalidae"] <- "Phytophage-parasitoid"   # Majority coding
D$Order.group <- D$Order
big_five <- c("Hymenoptera","Diptera","Coleoptera","Lepidoptera","Hemiptera")
idx <- which(!(D$Order.group %in% big_five))
D$Order.group[idx] <- "Other"

# Aggregate counts for niche, habitat and taxa
niche_comp <- aggregate(Sweden.estimated.total~Main.feeding.niche, data=D, FUN=sum)
colnames(niche_comp) <- c("niche","species")
habitat_comp <- aggregate(Sweden.estimated.total~Main.feeding.habitat, data=D, FUN=sum)
colnames(habitat_comp) <- c("habitat","species")
taxon_comp <- aggregate(Sweden.estimated.total~Order.group, data=D, FUN=sum)
colnames(taxon_comp) <- c("order","species")

# Write result files
write.table(niche_comp, "../data/niche_comp.tsv", row.names=FALSE, sep="\t")
write.table(habitat_comp, "../data/habitat_comp.tsv", row.names=FALSE, sep="\t")
write.table(taxon_comp, "../data/taxon_comp.tsv", row.names=FALSE, sep="\t")

# Break down niche and habitat by order group
niche_comp_by_order <- aggregate(Sweden.estimated.total~Main.feeding.niche+Order.group, data=D, FUN=sum)
colnames(niche_comp_by_order) <- c("niche","order","species")
habitat_comp_by_order <- aggregate(Sweden.estimated.total~Main.feeding.habitat+Order.group, data=D, FUN=sum)
colnames(habitat_comp_by_order) <- c("habitat","order","species")

# Write result files
write.table(niche_comp_by_order, "../data/niche_comp_by_order.tsv", row.names=FALSE, sep="\t")
write.table(habitat_comp_by_order, "../data/habitat_comp_by_order.tsv", row.names=FALSE, sep="\t")

# Update estimate by assuming that 50% of Diptera and Hymenoptera are sampled in the IBA inventory
T <- readRDS("../../iba_data/cluster_taxonomy_se.rds")
target_hym_OTUs <- sum(T$Order=="Hymenoptera") * 2
target_dip_OTUs <- sum(T$Order=="Diptera") * 2

estimated_hym_OTUs <- sum(D$Sweden.estimated.total[D$Order=="Hymenoptera"])
estimated_dip_OTUs <- sum(D$Sweden.estimated.total[D$Order=="Diptera"])

orig_hym_OTUs <- sum(D$Corrected.values.2017[D$Order=="Hymenoptera"])
orig_dip_OTUs <- sum(D$Corrected.values.2017[D$Order=="Diptera"])

D$prop_unknown <- D$Sweden.estimated.total/D$Corrected.values.2017
idx <- which(D$Corrected.values.2017==0)
D$prop_unknown[idx] <- 1.0

convert_factor <- function(f,newN,oldN,origN) {
    k <- (newN-origN) / (oldN-origN)
    return (1.0 + k*(f-1.0))
}

D$new_factor <- D$prop_unknown
D$new_estimate <- D$Sweden.estimated.total

idx <- which(D$Order=="Hymenoptera")
D$new_factor[idx] <- convert_factor(D$prop_unknown[idx], target_hym_OTUs, estimated_hym_OTUs, orig_hym_OTUs)
D$new_estimate[idx] <- D$new_factor[idx]*D$Corrected.values.2017[idx]

idx <- which(D$Order=="Diptera")
D$new_factor[idx] <- convert_factor(D$prop_unknown[idx], target_dip_OTUs, estimated_dip_OTUs, orig_dip_OTUs)
D$new_estimate[idx] <- D$new_factor[idx]*D$Corrected.values.2017[idx]

# Aggregate counts for niche, habitat and taxa
niche_comp <- aggregate(new_estimate~Main.feeding.niche, data=D, FUN=sum)
colnames(niche_comp) <- c("niche","species")
habitat_comp <- aggregate(new_estimate~Main.feeding.habitat, data=D, FUN=sum)
colnames(habitat_comp) <- c("habitat","species")
taxon_comp <- aggregate(new_estimate~Order.group, data=D, FUN=sum)
colnames(taxon_comp) <- c("order","species")

# Write result files
write.table(niche_comp, "../data/new_niche_comp.tsv", row.names=FALSE, sep="\t")
write.table(habitat_comp, "../data/new_habitat_comp.tsv", row.names=FALSE, sep="\t")
write.table(taxon_comp, "../data/new_taxon_comp.tsv", row.names=FALSE, sep="\t")

