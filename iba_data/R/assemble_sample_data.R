# Assemble sample data (lat, long, time, habitat)
library(lubridate)

# TODO: Set path to the directory for iba metadata on your system
path <- "~/dev/figshare-repos/iba/raw_data/v6/"

meta_mg <- read.delim(paste0(path,"CO1_sequencing_metadata_MG.tsv"))
meta_se <- read.delim(paste0(path,"CO1_sequencing_metadata_SE.tsv"))
malaise_samples_mg <- read.delim(paste0(path,"samples_metadata_malaise_MG.tsv"))
malaise_samples_se <- read.delim(paste0(path,"samples_metadata_malaise_SE.tsv"))
litter_samples_mg <- read.delim(paste0(path,"samples_metadata_litter_MG.tsv"))
litter_samples_se <- read.delim(paste0(path,"samples_metadata_soil_litter_SE.tsv"))
sites_mg <- read.delim(paste0(path,"sites_metadata_MG.tsv"))
sites_se <- read.delim(paste0(path,"sites_metadata_SE.tsv"))

# Read in counts data
counts_malaise_mg <- read.delim("cluster_counts_malaise_mg.tsv")
counts_malaise_se <- read.delim("cluster_counts_malaise_se.tsv")
counts_litter_mg <- read.delim("cluster_counts_litter_mg.tsv")
counts_litter_se <- read.delim("cluster_counts_litter_se.tsv")

set1 <- colnames(counts_malaise_mg)[-1]
set2 <- colnames(counts_malaise_se)[-1]
set3 <- colnames(counts_litter_mg)[-1]
set4 <- colnames(counts_litter_se)[-1]

D1 <- data.frame(list(sampleID_NGI=set1))
D2 <- data.frame(list(sampleID_NGI=set2))
D3 <- data.frame(list(sampleID_NGI=set3))
D4 <- data.frame(list(sampleID_NGI=set4))

add_malaise_sample_meta <- function(D, meta, samples, sites) {
    D$sampleID_FIELD <- meta$sampleID_FIELD[match(D$sampleID_NGI,meta$sampleID_NGI)]
    idx <- match(D$sampleID_FIELD,samples$sampleID_FIELD)
    D$trapID <- samples$trapID[idx]
    D$placing_date <- dmy(samples$placing_date[idx])
    D$collecting_date <- dmy(samples$collecting_date[idx])
    D$mid_date <- D$placing_date + (D$collecting_date - D$placing_date)/2
    idx <- match(D$trapID,sites$trapID)
    D$trap_habitat <- sites$trap_habitat[idx]
    D$latitude <- sites$latitude_WGS84[idx]
    D$longitude <- sites$longitude_WGS84[idx]
    return(D)
}

D1 <- add_malaise_sample_meta(D1, meta_mg, malaise_samples_mg, sites_mg)
D2 <- add_malaise_sample_meta(D2, meta_se, malaise_samples_se, sites_se)
D1 <- D1[complete.cases(D1),]
D2 <- D2[complete.cases(D2),]
write.tsv <- function(D,file) { write.table(D,file,row.names=FALSE,sep="\t") }
write.tsv(D1,"malaise_sample_meta_mg.tsv")
write.tsv(D2,"malaise_sample_meta_se.tsv")

add_litter_sample_meta <- function(D, meta, samples, sites) {
    D$sampleID_FIELD <- meta$sampleID_FIELD[match(D$sampleID_NGI,meta$sampleID_NGI)]
    idx <- match(D$sampleID_FIELD,samples$sampleID_FIELD)
    D$trapID <- samples$trapID[idx]
    D$date <- samples$date[idx]
    idx <- match(D$trapID,sites$trapID)
    D$trap_habitat <- sites$trap_habitat[idx]
    D$latitude <- sites$latitude_WGS84[idx]
    D$longitude <- sites$longitude_WGS84[idx]
    return (D)
}

D3 <- add_litter_sample_meta(D3, meta_mg, litter_samples_mg, sites_mg)
D4 <- add_litter_sample_meta(D4, meta_se, litter_samples_se, sites_se)
D3$date <- mdy(D3$date)     # Note different date format in this dataset
D4$date <- dmy(D4$date)
D3[complete.cases(D3),]
D4[complete.cases(D4),]
write.tsv(D3,"litter_sample_meta_mg.tsv")
write.tsv(D4,"litter_sample_meta_se.tsv")

# Aggregate litter and malaise sample metadata
D1$date <- D1$mid_date
D2$date <- D2$mid_date
idx <- which(colnames(D1) %in% c("placing_date","collecting_date", "mid_date"))
D1 <- D1[,-idx]
idx <- which(colnames(D2) %in% c("placing_date","collecting_date", "mid_date"))
D2 <- D2[,-idx]
D1$sample_type <- "malaise"
D2$sample_type <- "malaise"
D3$sample_type <- "litter"
D4$sample_type <- "litter"
write.tsv(rbind(D1,D3),"malaise_litter_sample_meta_mg.tsv")
write.tsv(rbind(D2,D4),"malaise_litter_sample_meta_se.tsv")

