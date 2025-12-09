# Script correcting the Sundh et al (2024) files ('chesters_new_outgroups...')
# generating updated files ('chesters_new_outgroups_updated...').

# Libraries needed
library(ape)

# Define sets of errors to be removed, identified by TipLabel
# See the README file in this directory and the contents in
# the clean_seq_data folder.

# Apparently correct tabanid sequence placed in the Orthoptera
# in the Chesters (2017) treee and also in Sundh et al (2024)
placement_error <- c("Lilaea_fuliginosa")

# Class hit errors                                      # REMARKS [Chesters placement] -- Real origin (according to blastn hits)
class_hit_errors <- c("Aphytis_africanus",              # [Insecta: Hymenoptera] -- Spider sequence (Arachnida:Araneae:Zodariidae)
                      "Apiophora_quadricinctata",       # [Insecta: Diptera] -- Wolbachia sequence
                      "Arantia_congensis",              # [Insecta: Orthoptera] -- Human sequence
                      "Dichotomius_bos",                # [Insecta: Coleoptera] -- Cichlid sequence (Actinopteri:Cichliformes:Cichlidae)
                      "Dichotomius_geminatus",          # [Insecta: Coleoptera] -- Fish sequence (Actinopteri:Characiformes:Characidae)
                      "Dichotomius_laevicollis",        # [Insecta: Coleoptera] -- Fish sequence (Actinopteri:Siluriformes:Loricariidae)
                      "Dichotomius_nisus",              # [Insecta: Coleoptera] -- Human sequence
                      "Dichotomius_semisquamosus",      # [Insecta: Coleoptera] -- Fish sequence (Actinopteri:Siluriformes:Loricariidae)
                      "Dichotomius_sericeus",           # [Insecta: Coleoptera] -- Fish sequence (Actinopteri:Characiformes:Characidae)
                      "Eusilpha_jakowlewi",             # [Coleoptera: Silphidae] -- ?Uroobovella (Arachnida:Mesostigmata:Urodinychidae) (no match > 83%)
                      "Hoplitis_producta",              # [Insecta: Hymenoptera] -- Xenopus sequence (Amphibia:Anura:Pipidae)
                      "Mycetophila_lunata",             # [Insecta: Diptera] -- Tisiphia endosymbiont (Alphaproteobacteria:Rickettsiales:Rickettsiaceae)
                      "Sophrops_striata",               # [Insecta: Coleoptera] -- Crustacean sequence (Malacostraca:Decapoda:Portunidae)
                      "Synuchus_cycloderus",            # [Inecta: Coleoptera] -- Crustacean sequence (Malacostraca:Decapoda:Porcellanidae)
                      "Timia_nigripes",                 # [Insecta: Diptera] -- Human sequence
                      "Neocondeellum_brachytarsum",     # [Protura: Acerentomata] -- Human sequence
                      "Fujientomon_dicestum"            # [Protura: Sinentomata] -- Fungal sequence (Malasseziomycetes:Malasseziales:Malasseziaceae)
                      )

# Order hit errors (excluding class hit errors)         # REMARKS [Chesters placement] -- Real origin (according to blastn hits)
order_hit_errors <- c("Anoplocnemis_phasianus",         # [Hempitera: Coreidae] -- Sitophilus zeamais equence (Coleoptera: Curculionidae)  
                      "Chaeridiona_metallica",          # [Coleoptera: Chrysomelidae] -- Neohalys (or related) (Hemiptera: Pentatomidae)
                      "Criomorphus_conspicuus",         # [Hemiptera: Delphacidae] -- Trichopalpus punctipes (Diptera: Scatophagidae)
                      "Enicmus_histrio",                # [Coleoptera: Latridiidae] -- Atrichopogon fusculus (Diptera: Ceratopogonidae)
                      "Epicephala_mirivalvata",         # [Lepidoptera: Gracillariidae] -- Eysarcoris annamita (Hemiptera: Pentatomidae)
                      "Fortunia_byrrhoides",            # [Hemiptera: Issidae] -- Cerapanorpa (Mecoptera: Panorpidae)
                      "Gastrimargus_marmoratus",        # [Orthoptera: Acrididae] -- Mantis religiosa (Mantodea: Mantidae)
                      "Helle_longirostris",             # [Diptera: Acroceridae] -- Climacia areolaris (Neuroptera: Sisyriidae) (from another study by the same group...)
                      "Lasia_carbonanicus",             # [Diptera: Acroceridae] -- Climacia areolaris (Neuroptera: Sisyriidae) (from another study by the same group...)
                      "Lasioderma_kiesenwetteri",       # [Coleoptera: Ptinidae] -- Dichopygina ramosa (Diptera: Sciaridae)
                      "Margarinotus_carbonarius",       # [Coleoptera: Histeridae] -- Prosimulium tomosvaryi (Diptera: Simuliidae)
                      "Nesobasis_caerulecaudata",       # [Odonata: Coenagrionidae] -- Rhopalosiphum maidis (Hemiptera: Aphididae)
                      "Opilo_whitei",                   # [Coleoptera: Cleridae] -- Chironomus (Diptera: Chironomidae)
                      "Potamarcha_obscura",             # [Odonata: Libellulidae] -- Pieris canidia (Lepidoptera: Pieridae)
                      "Priacma_serrata",                # [Coleoptera: Cupedidae] -- Probably error: only other close match is Zeugomantispa minuta (see this species...)
                      "Rhagovelia_zela",                # [Hemiptera: Veliidae] -- Chironomidae (Insecta: Diptera) (no match > 89%)
                      "Rhynocoris_fuscipes",            # [Hemiptera: Reduviidae] -- Spodoptera (Lepidoptera: Noctuidae)
                      "Rhynocoris_kumarii",             # [Hemiptera: Reduviidae] -- Spodoptera (Lepidoptera: Noctuidae)
                      "Tholymis_citrina",               # [Odonata: Libellulidae] -- Polyommatus (Lepidoptera: Lycaenidae)
                      "Thyllis_crassa",                 # [Diptera: Acroceridae] -- Climacia areolaris (Neuroptera: Sisyriidae) (from another study by the same group...)
                      "Zeugomantispa_minuta",           # [Neuroptera: Mantispidae] -- Probably error: only other close match is Priacma serrata (see this species...)
                      "Xenopsylla_cheopis"              # [Siphonaptera: Pulicidae] -- Hemipyrellia ligurriens (Diptera: Calliphoridae)
                     )

problem_tips <- c(placement_error, class_hit_errors, order_hit_errors)


# Correct the taxonomic annotations in the modified Chesters tree
# from Sundh et al. (2024)
# ===============================================================

# Read in original data
D <- read.delim("chesters_new_outgroups_taxonomy.tsv")

# Remove erroneous sequences
D <- D[!(D$TipLabel %in% problem_tips),]

# Correct family names (errors discovered when processing trait data)
D$Family[D$Family=="Xylophagaidae"] <- "Xylophagidae"
D$Family[D$Family=="Pemphigidae"] <- "Aphididae"
D$Family[D$Family=="Kerriidae"] <- "Tachardiidae"
D$Family[D$Family=="Synneuridae"] <- "Canthyloscelidae"

# Update family names for Cynipoidea
D$Family[D$Genus=="Diplolepis"] <- "Diplolepididae"
D$Family[D$Genus=="Liebelia"] <- "Diplolepididae"
D$Family[D$Genus=="Cecinothofagus"] <- "Paraulacidae"

# Write updated taxonomy file
write.table(D, "chesters_new_outgroups_taxonomy_updated.tsv", row.names=FALSE, sep="\t")


# Update the sequence data
# ========================

# Read in original data
seqs <- read.FASTA("chesters_new_outgroups.fasta")

# Remove problematic sequences
seqs <- seqs[!(names(seqs) %in% problem_tips)]

# Write updated sequence file
write.FASTA(seqs,"chesters_new_outgroups_updated.fasta")


# Update the tree
# ===============

# Read in original tree
tree <- read.tree("chesters_new_outgroups.nwk")

# Drop tips
tree <- drop.tip(tree, problem_tips)

# Write updated tree
write.tree(tree, "chesters_new_outgroups_updated.nwk")


