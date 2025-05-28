# Script testing compliance of Chesters tree with preferred
# backbone (Fig. 3 in Chesters, 2017)

library(ape)

D <- read.delim("chesters_new_outgroups_taxonomy.tsv")
T <- read.tree ("chesters_new_outgroups.nwk")

D <- D[D$TipLabel!="Lilaea_fuliginosa",]
T <- drop.tip(T,"Lilaea_fuliginosa")

# Check monophyly function
check_mon <- function(T,clade,tips) {
    cat(clade,"--", length(tips), "spp --", is.monophyletic(T,tips), "\n")
}

# Check monophyly of all insect orders
orders <- table(D$Order[D$Class=="Insecta"])
x <- list()
for (ord in names(orders))
    x <- c(x,list(D$TipLabel[D$Order==ord]))
names(x) <- names(orders)

cat("\nMonophyly of Insecta orders (after removal of outgroups,\n")
cat ("poor-quality sequences, and Lilaea_fuliginosa)\n")
cat("\nOrder -- spp -- monophyletic?\n")
cat("-----------------------------\n")
for (ord in names(orders)) check_mon(T,ord,x[[ord]])

# Check monohyly of higher clades
cat("\nMonophyly of higher clades\n")
cat("\nClade -- spp -- monophyletic?\n")
cat("-----------------------------\n")

# Get order tips as building blocks
archaeognatha <- D$TipLabel[D$Order=="Archaeognatha"]
zygentoma <- D$TipLabel[D$Order=="Zygentoma"]
ephemeroptera <- D$TipLabel[D$Order=="Ephemeroptera"]
odonata <- D$TipLabel[D$Order=="Odonata"]
zoraptera <- D$TipLabel[D$Order=="Zoraptera"]
dermaptera <- D$TipLabel[D$Order=="Dermaptera"]
plecoptera <- D$TipLabel[D$Order=="Plecoptera"]
orthoptera <- D$TipLabel[D$Order=="Orthoptera"]
blattodea <- D$TipLabel[D$Order=="Blattodea"]
mantodea <- D$TipLabel[D$Order=="Mantodea"]
embioptera <- D$TipLabel[D$Order=="Embioptera"]
phasmatodea <- D$TipLabel[D$Order=="Phasmatodea"]
grylloblattodea <- D$TipLabel[D$Order=="Grylloblattodea"]
mantophasmatodea <- D$TipLabel[D$Order=="Mantophasmatodea"]
psocoptera <- D$TipLabel[D$Order=="Psocoptera"]
phthiraptera <- D$TipLabel[D$Order=="Phthiraptera"]
thysanoptera <- D$TipLabel[D$Order=="Thysanoptera"]
hemiptera <- D$TipLabel[D$Order=="Hemiptera"]
hymenoptera <- D$TipLabel[D$Order=="Hymenoptera"]
lepidoptera <- D$TipLabel[D$Order=="Lepidoptera"]
trichoptera <- D$TipLabel[D$Order=="Trichoptera"]
diptera <- D$TipLabel[D$Order=="Diptera"]
mecoptera <- D$TipLabel[D$Order=="Mecoptera"]
siphonaptera <- D$TipLabel[D$Order=="Siphonaptera"]
strepsiptera <- D$TipLabel[D$Order=="Strepsiptera"]
coleoptera <- D$TipLabel[D$Order=="Coleoptera"]
raphidioptera <- D$TipLabel[D$Order=="Raphidioptera"]
neuroptera <- D$TipLabel[D$Order=="Neuroptera"]
megaloptera <- D$TipLabel[D$Order=="Megaloptera"]

neu.meg <- c(neuroptera, megaloptera)
neuropterida <- c(neuroptera,megaloptera,raphidioptera)
coleopterida <- c(coleoptera,strepsiptera)
neuropteroidea <- c(neuropterida,coleopterida)
amphiesmenoptera <- c(lepidoptera,trichoptera)
mec.sip <- c(mecoptera,siphonaptera)
antliophora <- c(diptera,mec.sip)
panorpida <- c(antliophora,amphiesmenoptera)
aparaglossata <- c(panorpida,neuropteroidea)
holometabola <- c(hymenoptera,aparaglossata)

n <- getMRCA(T,c("Bemisia_tabaci","Aphis_fabae"))
homoptera <- (extract.clade(T,n))$tip.label
n <- getMRCA(T,c("Acanthosoma_giganteum","Velia_caprai"))
heteroptera_ex_ful <- (extract.clade(T,n))$tip.label
n <- getMRCA(T,c("Xenophyes_cascus","Okanagana_lurida"))
fulgoromorpha <- D$TipLabel[D$Order=="Hemiptera" & !(D$TipLabel %in% c(homoptera,heteroptera_ex_ful))]
heteroptera <- c(heteroptera_ex_ful,fulgoromorpha)

condylognatha <- c(thysanoptera, hemiptera)
psocodea <- c(psocoptera,phthiraptera)
paraneoptera <- c(condylognatha, psocodea)
c1 <- c(paraneoptera,holometabola)

notoptera <- c(mantophasmatodea,grylloblattodea)
emb.pha <- c(embioptera,phasmatodea)
emb.pha.not <- c(emb.pha,notoptera)
dictyoptera <- c(mantodea,blattodea)
emb.pha.not.dic <- c(emb.pha.not,dictyoptera)
ort.ple <- c(orthoptera,plecoptera)
ort.ple.emb.pha.not.dic <- c(ort.ple,emb.pha.not.dic)
der.zor <- c(dermaptera,zoraptera)
polyneoptera <- c(der.zor,ort.ple.emb.pha.not.dic)

neoptera <- c(polyneoptera,paraneoptera,holometabola)
paleoptera <- c(odonata,ephemeroptera)
pterygota <- c(neoptera,paleoptera)
dicondylia <- c(pterygota,zygentoma)
insecta <- c(archaeognatha,dicondylia)

check_mon(T,"neu.meg",neu.meg)
check_mon(T,"neuropterida",neuropterida)
check_mon(T,"coleopterida",coleopterida)
check_mon(T,"neuropteroidea (coleopterida+neuropterida)",neuropteroidea)
check_mon(T,"amphiesmenoptera",amphiesmenoptera)
check_mon(T,"mec.sip",mec.sip)
check_mon(T,"antliophora",antliophora)
check_mon(T,"panorpida (antliophora+amphiesmenoptera)",panorpida)
check_mon(T,"aparaglossata (panorpida+neuropteroidea)",aparaglossata)
check_mon(T,"holometabola (hymenoptera+aparaglossata)",holometabola)
check_mon(T,"homoptera",homoptera)
check_mon(T,"heteroptera_ex_ful",heteroptera_ex_ful)
check_mon(T,"fulgoromorpha",fulgoromorpha)
check_mon(T,"heteroptera",heteroptera)
check_mon(T,"condylognatha",condylognatha)
check_mon(T,"psocodea",psocodea)
check_mon(T,"paraneoptera",paraneoptera)
check_mon(T,"paraneoptera+holometabola",c1)
check_mon(T,"notoptera",homoptera)
check_mon(T,"emb.pha",emb.pha)
check_mon(T,"emb.pha.not",emb.pha.not)
check_mon(T,"dictyoptera",dictyoptera)
check_mon(T,"emb.pha.not.dic",emb.pha.not.dic)
check_mon(T,"ort.ple",ort.ple)
check_mon(T,"ort.ple.emb.pha.not.dic",ort.ple.emb.pha.not.dic)
check_mon(T,"der.zor",der.zor)
check_mon(T,"polyneoptera",polyneoptera)
check_mon(T,"neoptera",neoptera)
check_mon(T,"paleoptera",paleoptera)
check_mon(T,"pterygota",pterygota)
check_mon(T,"dicondylia",dicondylia)
check_mon(T,"insecta",insecta)
cat("\n")
