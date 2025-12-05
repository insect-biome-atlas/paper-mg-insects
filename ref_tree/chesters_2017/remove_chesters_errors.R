# This script removes tips in the Chesters (2017) tree as
# modified by Sundh et al (2024), which are likely errors

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
                      "Xenopsylla_cheopis",             # [Siphonaptera: Pulicidae] -- Hemipyrellia ligurriens (Diptera: Calliphoridae)
                     )
