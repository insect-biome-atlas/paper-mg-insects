# Merge Chalcidoidea taxonomy data and sequences

# We work with two scenarios:
#
# 1. Match sequences to tips in Cruaud et al, drop unmatched tips.
#    Add Mymaromella mira as extra outgroup.
#
# 2. As above, but add all Chesters sequences of species that
#    are not present in the Cruaud-derived dataset.
#

library(ape)

# 1. Merge taxonomy files
# -----------------------

T_cr <- read.delim("cruaud_CO1_taxonomy.tsv")
T_ch <- read.delim("chesters_chalcidoidea_taxonomy_updated.csv",sep=";")

# Update Chesters taxonomy information
T_ch$TipLabel <- T_ch$TipLabel_updated
T_ch$Family <- T_ch$Family_updated
T_ch$Subfamily <- T_ch$Subfamily_updated
T_ch$Tribe <- T_ch$Tribe_updated
T_ch$Clade <- T_ch$Clade_updated
problem_seqs <- c("Aphytis_africanus",         # This sequence now deleted (likely a spider sequence)
                  "Syrphophagus_aphidivorus",  # (Encyrtidae) -- a mymarid sequence
                  "Pseudleptomastix_mexicana", # (Encyrtidae) -- a Telenomus sequence
                  "Brachymeria_compsilurae",   # (Chalcididae2) -- an ant sequence
                  "Pediobomyia_canaliculata",  # (Eulophidae) -- No close match in eulophids
                  "Rhynchentedon_maximus",     # (Eulophidae) -- No close match in eulophids
                  "Acerophagus_papayae"        # (Encyrtidae) -- No close match in encyrtids, odd sequence
                  )
T_ch <- T_ch[!(T_ch$TipLabel %in% problem_seqs),]

# Update taxonomy information at genus and species level for two cases
T_ch$Genus[T_ch$TipLabel=="Valisia_esquirolianae"] <- "Valisia"
T_ch$Genus[T_ch$TipLabel=="Burkseus_vittatus"] <- "Burkseus"
T_ch$Species[T_ch$TipLabel=="Valisia_esquirolianae"] <- "Valisia esquirolianae"
T_ch$Species[T_ch$TipLabel=="Burkseus_vittatus"] <- "Burkseus vittatus"

# Add higher level taxonomic info
T_cr$Kingdom <- "Animalia"
T_cr$Phylum <- "Arthropoda"
T_cr$Class <- "Insecta"
T_cr$Order <- "Hymenoptera"

T_ch$Kingdom <- "Animalia"
T_ch$Phylum <- "Arthropoda"
T_ch$Class <- "Insecta"
T_ch$Order <- "Hymenoptera"

cols <- c("TipLabel","Kingdom","Phylum","Class","Order","Family","Subfamily","Tribe","Genus","Species","Clade")

T_ch1 <- T_ch[T_ch$Family=="Mymarommatidae",]
T_ch2 <- T_ch[!(T_ch$TipLabel %in% T_cr$TipLabel),]

T1 <- rbind(T_cr[,cols], T_ch1[,cols])
T2 <- rbind(T_cr[,cols], T_ch2[,cols])

write.table(T1,"expanded_chalcidoidea1_taxonomy.tsv",sep="\t",row.names=FALSE)
write.table(T2,"expanded_chalcidoidea2_taxonomy.tsv",sep="\t",row.names=FALSE)


# 2. Merge sequence files
# -----------------------

# Set data path
data_path <- ""

# Read in sequences
chesters_seqs <- read.FASTA("chesters_chalcidoidea.fasta")
cruaud_seqs <- read.FASTA("cruaud_CO1.fasta")

# Correct sequence names
names(chesters_seqs)[which(names(chesters_seqs)=="Blastophaga_esquirolianae")] <- "Valisia_esquirolianae"
names(chesters_seqs)[which(names(chesters_seqs)=="Cirrospilus_vittatus")] <- "Burkseus_vittatus"

chesters_seqs1 <- chesters_seqs[names(chesters_seqs) %in% T_ch1$TipLabel]
chesters_seqs2 <- chesters_seqs[names(chesters_seqs) %in% T_ch2$TipLabel]

seqs1 <- c(cruaud_seqs, chesters_seqs1)
seqs2 <- c(cruaud_seqs, chesters_seqs2)

# Write merged sequence files
write.FASTA(seqs1,"expanded_chalcidoidea1.fasta")
write.FASTA(seqs2,"expanded_chalcidoidea2.fasta")

