# Assemble co1 sequences from Sikora et al (2019) and associated taxonomy information

# List of accession numbers (from their Table 1)
sikora_co1_seqs <- c("KT316846.1","KT316858.1","KT316876.1","KX453757.1","KT316834.1","KT316840.1","MG684785.1","KT316839.1","KT316851.1","KT316867.1","MG684786.1","KT316845.1","KT316853.1","KT316873.1","MG684848.1","MG684804.1","MG684840.1","KX453761.1","MG684799.1","MG684798.1","MG684800.1","MG684834.1","MG684821.1","MG684802.1","MG684803.1","MG684801.1","MG684791.1","MG684795.1","MG684806.1","KT316860.1","EU375697.1","MG684809.1","MG684793.1","MG684794.1","MG684805.1","MG684796.1","MG684792.1","MG684843.1","MG684797.1","MG684807.1","MG684808.1","MG684829.1","KT316850.1","MG684830.1","MG684835.1","MG684790.1","MG684814.1","MG684815.1","MG684813.1","MG684819.1","KT316837.1","MG684811.1","MG684817.1","MG684839.1","MG684822.1","MG684812.1","MG684841.1","MG684842.1","MG684820.1","MG684818.1","MG684816.1","MG684787.1","MG684846.1","MG684788.1","MG684789.1","MG684810.1","KT316859.1","MG684823.1","MG684824.1","MG684836.1","MG684844.1","MG684827.1","MG684826.1","MG684832.1")

# Lists of taxonomy info (from their Table 1)
sikora_species <- c("Bibio marci",
                    "Penthetria funebris",
                    "Ditomyia fasciata",
                    "Symmerus annulatus",
                    "Keroplatus testaceus",
                    "Orfelia nemoralis",
                    "Exechia fusca",
                    "Gnoriste bilineata",
                    "Dolichosciara flavipes",
                    "Zygoneura sciarina",
                    "Heterotricha takkae",
                    "Insulatricha hippai",
                    "Nepaletricha sigma",
                    "Catotricha subobsoleta",
                    "Trichoceromyia oregonensis",
                    "Acodiplosis inulae",
                    "Asphondylia pruniperda",
                    "Asphondylia sarothamni",
                    "Contarinia asclepiadis",
                    "Contarinia loti",
                    "Lestodiplosis polypori",
                    "Loewiola centaureae",
                    "Obolodiplosis robinae",
                    "Odontodiplosis sp.",
                    "Placochela nigripes",
                    "Schizomyia galiorum",
                    "Dasineura kellneri",
                    "Dasineura trifolii",
                    "Lasioptera carophila",
                    "Lasioptera rubi",
                    "Mayetiola destructor",
                    "Mikiola fagi",
                    "Oligotrophus juniperinus",
                    "Oligotrophus panteli",
                    "Ozirhincus longicollis",
                    "Rabdophaga heterobia",
                    "Rhopalomyia tanaceticola",
                    "Rondaniola bursaria",
                    "Spurgia euphorbiae",
                    "Didactylomyia longimana",
                    "Stomatosema obscurum",
                    "Anaretella iola",
                    "Lestremia cinerea",
                    "Lestremia leucophaea",
                    "Wasmanniella sp.",
                    "Acoenonia europaea",
                    "Aprionus brachypterus",
                    "Aprionus spiniger",
                    "Bryomyia apsectra",
                    "Campylomyza flavipes",
                    "Catocha angulata",
                    "Catocha incisa",
                    "Heterogenella bigibbata",
                    "unclassified.Micromyini",
                    "Monardia abnormis",
                    "Neurolyga excavata",
                    "Peromyia fungicola",
                    "Peromyia sp.",
                    "Polyardis silvalis",
                    "Skuhraviana triangulifera",
                    "Tropaprionus sp.",
                    "Camptomyia flavocinerea",
                    "Claspettomyia sp.",
                    "Dicerura dentata",
                    "Dicerura sp.",
                    "Dirhiza lateritia",
                    "Porricondyla nigripennis",
                    "Svenartia spungisi",
                    "Tetraneuromyia hirticornis",
                    "Zaitzeviola rufocinerea",
                    "unclassified.Diallactiini",
                    "Heteropeza pygmea",
                    "Leptosyna nervosa",
                    "Winnertzia nigripennis"
                    )

sikora_genus <- c(  "Bibio",
                    "Penthetria",
                    "Ditomyia",
                    "Symmerus",
                    "Keroplatus",
                    "Orfelia",
                    "Exechia",
                    "Gnoriste",
                    "Dolichosciara",
                    "Zygoneura",
                    "Heterotricha",
                    "Insulatricha",
                    "Nepaletricha",
                    "Catotricha",
                    "Trichoceromyia",
                    "Acodiplosis",
                    "Asphondylia",
                    "Asphondylia",
                    "Contarinia",
                    "Contarinia",
                    "Lestodiplosis",
                    "Loewiola",
                    "Obolodiplosis",
                    "Odontodiplosis",
                    "Placochela",
                    "Schizomyia",
                    "Dasineura",
                    "Dasineura",
                    "Lasioptera",
                    "Lasioptera",
                    "Mayetiola",
                    "Mikiola",
                    "Oligotrophus",
                    "Oligotrophus",
                    "Ozirhincus",
                    "Rabdophaga",
                    "Rhopalomyia",
                    "Rondaniola",
                    "Spurgia",
                    "Didactylomyia",
                    "Stomatosema",
                    "Anaretella",
                    "Lestremia",
                    "Lestremia",
                    "Wasmanniella",
                    "Acoenonia",
                    "Aprionus",
                    "Aprionus",
                    "Bryomyia",
                    "Campylomyza",
                    "Catocha",
                    "Catocha",
                    "Heterogenella",
                    "unclassified.Micromyini",
                    "Monardia",
                    "Neurolyga",
                    "Peromyia",
                    "Peromyia",
                    "Polyardis",
                    "Skuhraviana",
                    "Tropaprionus",
                    "Camptomyia",
                    "Claspettomyia",
                    "Dicerura",
                    "Dicerura",
                    "Dirhiza",
                    "Porricondyla",
                    "Svenartia",
                    "Tetraneuromyia",
                    "Zaitzeviola",
                    "unclassified.Diallactiini",
                    "Heteropeza",
                    "Leptosyna",
                    "Winnertzia"
                    )

sikora_fams <- c(rep("Bibionidae",2),
                 rep("Ditomyiidae",2),
                 rep("Keroplatidae",2),
                 rep("Mycetophilidae",2),
                 rep("Sciaridae",2),
                 rep("",3),
                 rep("Cecidomyiidae",61)
                )

# Note that the text suggests that the correct
# spelling is Micromyinae, and not Micromyiinae,
# as used in the table
sikora_subfams <-c(rep("",13),
                   rep("Catotrichinae",2),
                   rep("Cecidomyiinae",26),
                   rep("Lestremiinae",4),
                   rep("Micromyinae",16),
                   rep("Porricondylinae",9),
                   rep("Winnertziinae",4)
                  )

sikora_supertribes <- c(rep("",15),
                        rep("Cecidomyiidi",11),
                        rep("Lasiopteridi",13),
                        rep("Stomatosematidi",2),
                        rep("",33)
                       )

# Assemble sequences and write to file
library(rentrez)
outfile <- "sikora_CO1.fasta"
records <- entrez_fetch(db = "nuccore", id = sikora_co1_seqs, rettype = "fasta")
cat(records,file=outfile)

# Generate taxonomy file
taxonomy <- data.frame(list(GenBank=sikora_co1_seqs,
                            Family=sikora_fams,
                            Subfamily=sikora_subfams,
                            Supertribe=sikora_supertribes,
                            Genus=sikora_genus,
                            Species=sikora_species
                            )
                      )

# Figure out start and stop of coding sequence
library(ape)
source("../code/seq_fxns.R")

seqs <- read.FASTA(outfile)
start1 <- no_stop_codons(seqs, 1)
start2 <- no_stop_codons(seqs, 2)
start3 <- no_stop_codons(seqs, 3)

taxonomy$Start <- 1
taxonomy$Start[start2] <- 2
taxonomy$Start[start3] <- 3

taxonomy$Stop <- -1
for (i in 1:nrow(taxonomy)) {
    len <- length(seqs[[i]]) - taxonomy$Start[i] + 1

    taxonomy$Stop[i] <- length(seqs[[i]]) - (len %% 3)
}

write.table(taxonomy,"sikora_taxonomy.tsv",row.names=FALSE,sep="\t")

