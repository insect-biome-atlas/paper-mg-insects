# Notes on mycetophilid data

All observations pertain to 42 unlikely species occurrences in MG data reported by Jostein Kjaerandsen (the e-mail says 41 but the actual number is 42):

> Hei Fredrik,
>
> Jeg følger interessert med på de store datasettene dere har publisert på GBIF fra Sverige og Madagaskar.
>
> Det fra Sverige er jo imponerende med nesten 100 000 soppmygg tilhørende 589 arter. Av disse er 20 nye for Sverige og to er nye for Norden: Docosia morionella og Exechia absurda. Sistnevnte var en overraskelse da den bare er kjent fra Nord-Amerika, men det er jeg som har barkodet og identifisert et par eksemplarer fra Alaska av denne arten i BOLD, så den burde være rett assosiert.
>
> Tilvarende datasett fra norsk insektovervåkning er nå på over 200 000 registreringer av 657 arter over fire år med metabarcoding.
>
> Som jeg vurderer datasettene ser det meste bra ut og til å stole på, og jeg bruker nå det norske datasettet bl.a. som grunnlag for rødlistearbeid.
>
> Når jeg imidlertid sjekker dette leaf-litter datasettet fra Madagaskar så er det noe fundamentalt feil. <https://www.gbif.org/dataset/9dbba635-822e-4270-90ae-908e28ee6cc9> Her listes 41 nordiske arter av Mycetophilidae fra Madagaskar:
>
> Allocotocera pulchella 1
>
> Allodia triangularis 1
>
> Boletina edwardsi 2
>
> Boletina gripha 1
>
> Boletina griphoides 1
>
> Boletina nigricans 3
>
> Boletina onegensis 2
>
> Brevicornu sericoma 2
>
> Brevicornu verralli 10
>
> Coelosia tenella 3
>
> Cordyla flaviceps 2
>
> Cordyla parvipalpis 1
>
> Ectrepesthoneura hirta 2
>
> Exechia confinis 1
>
> Exechia dorsalis 1
>
> Exechia festiva 2
>
> Exechia fusca 6
>
> Exechia macula 1
>
> Exechia pseudocincta 1
>
> Leia cylindrica 1
>
> Leia winthemi 2
>
> Mycetophila ichneumonea 3
>
> Mycetophila signatoides 1
>
> Mycomya affinis 4
>
> Mycomya annulata 3
>
> Mycomya fimbriata 3
>
> Mycomya shermani 6
>
> Mycomya trilineata 3
>
> Opsion nigra 3
>
> Phronia egregia 1
>
> Rymosia domestica 3
>
> Rymosia fasciata 1
>
> Rymosia placida 3
>
> Rymosia signatipes 3
>
> Sciophila hirta 2
>
> Sciophila salassea 1
>
> Synapha vitripennis 1
>
> Syntemna hungarica 1
>
> Tarnania fenestralis 1
>
> Tarnania tarnanii 2
>
> Trichonta atricauda 1
>
> Zygomyia notata 2
>
> Det er utenkelig at alle disse artene finnes på Madagaskar, muligens kan en eller et par av de finnes der men det er svært usannsynlig. Så hva har skjedd her? Det må jo dreie seg om kontaminasjon med det svenske materialet.

-   There is one MG species on this list not present in SE data (based on rep ASV or consensus annotations). It is Brevicornu sericoma. There are two ASVs of it in the MG data and one of them occurs also in the SE data. This ASV is neither consistent with the consensus annotation nor is it the rep ASV.

-   There are three species on Jostein's list that are in neither SE or MG data. They are Mycomya trilineata, Opsion nigra, and Rymosia domestica. Requires research. Likely annotation ambiguities (Genus_X annotations in GBIF, updated annotations resolved presumably).

-   Of 38 unlikely Nordic species present in MG litter data, 12 are also present in MG positive controls.

-   For these 12 species, read ratios and distances are compatible with contamination (\<1/5 of reads in positive controls, at least 2 samples, distance 0) from positive controls only in a few cases: Boletina onegensis, Mycomya fimbriata, and Syntemna hungarica (although the latter in only 1 sample). Of these, only Mycomya fimbriata is present in mg_litter_control samples.

-   However, for the 16 (of 38) species for which there are rep ASVs that are exactly the same, almost half, all occur in \>5000 reads, 14 occur in \>10000 reads and 12 occur in more than 100,000 reads in the Swedish lysate data. But amplified DNA from those samples was not present at Station Linné(?)

-   Only 7 of the 38 species are present in mg_litter_control samples. Three of these are present in mg_litter_pos_con samples. One cluster is present in mg_litter_control samples but not in mg_litter samples nor in mg_litter_pos_con samples.

-   The contamination should come from amplified DNA. Could it come from the Norrland samples? Or is it more likely to come from pre-SIMS samples amplified at Station Linné?

-   Only one of the 38 species is found in mg_lysate samples, indicating that most of the species represent contamination. One could compute the overlap between Malaise trap samples and soil samples in Sweden to double-check this. For the 38 species considered here, all se_litter species (only one actually) were found in se_lysate samples.

-   There is one cluster, Mycetophilidae_cluster196 in MG, which is present in mg_litter_control samples, but only a very distant cluster (judging from rep ASV) is found in Swedish IBA samples. Is this a cluster found in the Norrland samples? Apparently not one of the mycetophilid barcodes from Norrland according to Mårten.
