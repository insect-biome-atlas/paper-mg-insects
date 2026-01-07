# Merge Ichneumonidae taxonomy data and sequences

# 1. Merge taxonomy files
# -----------------------

data_path <- ""
T1 <- read.delim(paste0(data_path,"chesters_ichneumonidae_taxonomy_updated.tsv"))
T2 <- read.delim(paste0(data_path,"missing_ichneumonidae_taxonomy.tsv"))

T2$TipLabel<-sub(" ","_",T2$Species)
T2$Kingdom <- "Animalia"
T2$Phylum <- "Arthropoda"
T2$Class <- "Insecta"
T2$Order <- "Hymenoptera"

cols <- c("TipLabel","Kingdom","Phylum","Class","Order","Family","Subfamily","Tribe","Genus","Species")

T <- rbind(T1[,cols],T2[,cols])

# Add info on Clade, matching life history data file

# Notes from Gavin Broad (in email 2024-03-29):

# Polysphincta group should be separate from Ephialtini as their biology is so different, as koinobiont ectoparasitoids of spiders
# In Chesters et al, the included polypshinctines are Acrodactyla, Acrotaphus, Dreisbachia, Zatypota.

# As Seraina said, the coverage of Cryptinae and Phygadeuontinae by Chesters will overlap very little with the Madagascan sequences.  This is a particular problem for Phygadeuontinae, which seem to be paraphyletic with respect to Ateleutinae, Cryptinae and Ichneumoninae.  There are many phygadeuontines on Madagascar, the majority in the Chirotica group, which don't seem to be represented in Chesters et al.  I am not sure whether these will easily be inferred to be phygadeuontines based on barcode matching.
# Anyway, please include Paraphylax as a key chiroticine.  If available, sequences of Astomaspis and Bodedia would be good.  These seem to be quite common on Madagascar.  Smaller genera include Gabia, Lienella and Mamelia.

T$Clade <- T$Subfamily
for (i in 1:nrow(T)) {
    if (grepl(T$Genus[i],"Acrodactyla|Acrotaphus|Dreisbachia|Zatypota")) {
        T$Clade[i] <- "Polysphincta_group"
    } else if (grepl(T$Genus[i],"Paraphylax|Astomaspis|Lienella")) {
        T$Clade[i] <- "Chirotica_group"
    }
}
T$Clade[T$Clade=="Pimplinae"] <- "Pimplinae_s_str"
T$Clade[T$Clade=="Phygadeuontinae"] <- "Phygadeuontinae_s_str"


cols <- c("TipLabel","Kingdom","Phylum","Class","Order","Family","Subfamily","Tribe","Clade","Genus","Species")
write.table(T[,cols],paste0(data_path,"expanded_ichneumonidae_taxonomy.tsv"),sep="\t",row.names=FALSE)


# 2. Merge sequence files
# -----------------------

# Read in needed functions
source("../code/seq_fxns.R")

# Set data path
data_path <- ""

# Read in sequences and taxonomy info.
chesters_seqs <- read.FASTA(paste0(data_path,"chesters_ichneumonidae.fasta"))
extra_seqs <- read.FASTA(paste0(data_path,"missing_ichneumonidae_CO1.fasta"))
chesters_meta <- T1
extra_meta <- T2

# Construct tip labels
extra_meta$TipLabel <- sub(" ", "_", extra_meta$Species)

# Relabel missing sequences from NCBI strings to tip labels
extra_seqs <- rename_seqs(extra_seqs, extra_meta)

# Extract the coding part
chesters_seqs <- extract_coding(chesters_seqs, chesters_meta)
extra_seqs <- extract_coding(extra_seqs, extra_meta)

# Write the resulting sequences
out_file <- paste0(data_path,"expanded_ichneumonidae.fasta")
write.FASTA(chesters_seqs, out_file)
write.FASTA(extra_seqs, out_file, append=TRUE)

