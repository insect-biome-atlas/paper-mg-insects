# Compute sampled diversity
D <- readRDS("../../sampling_bias/data/otu_site_meta_se.rds")
sampled_coleoptera <- length(unique(D$cluster[D$Order=="Coleoptera"]))
sampled_lepidoptera <- length(unique(D$cluster[D$Order=="Lepidoptera"]))

# Compute recorded diversity
T <- read.delim("../../traits/ronquist_2020_SE_traits.csv",sep=";")
recorded_coleoptera <- sum(T$Corrected.values.2017[T$Order=="Coleoptera"])
recorded_lepidoptera <- sum(T$Corrected.values.2017[T$Order=="Lepidoptera"])

bias <- (sampled_coleoptera/sampled_lepidoptera) / (recorded_coleoptera/recorded_lepidoptera)
cat("Estimated sampling bias for Coleoptera is:",bias,"\n")

