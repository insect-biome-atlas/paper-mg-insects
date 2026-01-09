# braconidae
Scripts and output data for generating an expanded and more accurate tree for braconids.

## Details
For Braconidae, we retrieved the original subfamily classification from NCBI and updated it using information from the UCE analysis of Jasso-Martínez et al. (2022; Figs. 3-8). The updated classification is in the file `Chesters_Braconidae_taxonomy_updated.xlsx` in the `source_data` directory, and exported in ';'-delimited format in `chesters_braconidae_taxonomy_updated.csv`.

For a few Doryctinae taxa, the placement into the Doryctinae_s_str or Doryctinae_South_America clades was based on preliminary analyses where they were allowed to float in the tree to find their affinities based on the CO1 data. See comments in the xlsx file `Chesters_Braconidae_taxonomy_updated.xlsx`.

We then manually added the available CO1 sequences of the  missing subfamilies. The information about these taxa and the corresponding fasta sequences is in the files:
 - `missing_braconidae_subfamilies_CO1.fasta`
 - `missing_braconidae_subfamilies_taxonomy.tsv`

Finally, we added the available CO1 sequences of the genera included in the UCE analysis of Jasso-Martínez et al. (2022) but not present in the Chesters tree. All UCE genera and their subfamily placement is given in the tsv file below. The genera are partly manually coded from the figures and the text, and partly from the tree given in the Insect Phylogeny Database (Chesters 2025), which is also included here.
 - `uce_braconidae_genera_taxonomy.tsv`: Taxonomy of the genera included in the UCE tree (including genera in the main text figures AND in the IPD tree, see next item).
 - `jassomartinez_2022.nwk`: Tree from Insect Phylogeny Database representing the UCE results (note small discrepancies between this tree and the figures in the main text)

We first added standard 658 bp CO1 barcode sequences of the missing genera using an R script. We then manually added genera with non-standard CO1 sequences. 
 - `add_missing_braconid_genera.R`: Script checking for the availability of standard barcode sequences for missing braconid genera, generating the two files below.
 - `missing_braconidae_genera_taxonomy.tsv`: Taxonomy and sequence accession data for genera with standard barcode sequences.
 - `missing_braconidae_genera_check_needed.tsv`: Taxonomy for genera for which a manual check is needed (CO1 sequences available, but not standard barcode).
 - `missing_braconidae_genera_checked_taxonomy.tsv`: Manually retrieved taxonomy and sequence accession data for genera with nonstandard barcode sequences.

We also manually retrieved full CO1 sequence data from mitochondrial genomes from two ichneumonid representatives suitable as outgroups. The taxonomy and sequence accession info are given in the file
 - `braconidae_outgroups_taxonomy.tsv`

The info in the taxonomy files is used to fetch the actual sequence data from GenBank. The taxonomy and sequence files are then merged.
 - `fetch_braconidae_sequences.R`: Script for fetching Braconidae sequences from GenBank.
 - `missing_braconidae_data.R`: Script merging sequence and taxonomy files, producing the `expanded_braconidae.fasta` and `expanded_braconidae_taxonomy.tsv` files.

The sequences were then aligned as amino acid sequences and the alignment converted back to a nucleotide alignment, and sites with 90% or more gaps were trimmed away (see `make_mb_data.sh` for the scripts used to achieve this). The resulting alignment is in `expanded_braconidae_aligned_trimmed.fasta`.

For the phylogenetic analysis, we used all higher clades (subfamily or above) with more than 95% bootstrap support in the UCE analysis as partial constraints. We added hard constraints for the families. The data and constraint files are generated in the script `generate_braconidae_nexus_files.R`, and the output files are `braconidae.nex` and `braconidae_constraints.nex`. They are placed in the `mb_runs` folder.

We ran the MrBayes analysis for 100 M generations using a strict clock model with a codon-partitioned GTR+Gamma model (see `mb_runs/run2/run.nex`).

## References
Belokoblylskij, Sergey. (2006) ‘Neoheterospilus Gen. n., a New Genus of the Tribe Heterospilini (Hymenoptera: Braconidae, Doryctinae) with Highly Modified Ovipositor and a Worldwide Distribution’. Insect Systematics & Evolution 37 (April 2006): 149–78. https://doi.org/10.1163/187631206788831119.

Belokobylskij, Sergey A., Alejandro Zaldívar-Riverón, and Juana M. Coronado-Blanco. (2014) ‘Phylogenetic Affinities of Monarea Szépligeti, 1904 (Hymenoptera: Braconidae, Doryctinae, with Description of a New Species from Mexico’. Zootaxa 3795, no. 4: 421–30. https://doi.org/10.11646/zootaxa.3795.4.2.

Belokobylskij, Sergey, Fadia Ceccarelli, and Alejandro Zaldívar‐riverón. (2011) ‘Kauriphanes n. Gen., a New Genus of Braconid Parasitoid Wasp (Hymenoptera: Braconidae: Doryctinae) from New Zealand’. Annales- Societe Entomologique de France 47 (September 2011): 394–401. https://doi.org/10.1080/00379271.2011.10697733.

Belokobylskij, Sergey, Ernesto Samaca-Sáenz, and Alejandro Zaldívar‐riverón. (2015) ‘Mexiare Gen. Nov., a New Doryctinae Genus (Hymenoptera: Braconidae) from Mexico with Fused First and Second Metasomal Tergites’. Zootaxa 3914 (January 2015): 122–30. https://doi.org/10.11646/zootaxa.3914.2.2.

Chesters D (2025) Where are the biggest gaps in phylogenetic coverage of insect diversity? Systematic Entomology 50: 221–36. https://doi.org/10.1111/syen.12652.

Iqbal, M., A. D. Austin, and S. A. Belokobylskij. (2006) ‘Systematics of the Australasian Endemic Wasp Genus Syngaster Brullé (Hymenoptera: Braconidae: Doryctinae)’. Journal of Natural History 40, nos 13–14: 819–53. https://doi.org/10.1080/00222930600790653.

Jasso-Martínez JM, Santos BF, Zaldívar-Riverón A, Fernández-Triana JL, Sharanowski BJ, Richter R et al. (2022) Phylogenomics of braconid wasps (Hymenoptera, Braconidae) sheds light on classification and the evolution of parasitoid life history traits. Molecular Phylogenetics and Evolution 173: 107452. https://doi.org/10.1016/j.ympev.2022.107452.

Martínez, Juan José, Fadia Ceccarelli, and Alejandro Zaldívar‐riverón. (2010) ‘The Genus Iare Barbalho and Penteado-Dias (Hymenoptera: Braconidae: Doryctinae) in Mexico, with the Description of Two New Species’. Zootaxa 2685 (November 2010): 30–38. https://doi.org/10.11646/zootaxa.2685.1.2.

van Noort, S. 2025. WaspWeb: Hymenoptera of the World. URL: www.waspweb.org (accessed on 2026-01-06).

Marsh, Paul M. (1988) ‘Revision of the Tribe Odontobraconini in the Western Hemisphere (Hymenoptera: Braconidae: Doryctinae)’. Systematic Entomology 13, no. 4: 443–64. https://doi.org/10.1111/j.1365-3113.1988.tb00255.x.

Nunes, Juliano, Alejandro Zaldivar-Riveron, Clóvis Castro, et al. (2012) ‘Doryctopambolus Nunes & Zaldívar-Riverón (Braconidae), a New Neotropical Doryctine Wasp Genus with Propodeal Spines’. ZooKeys 223 (September 2012): 53–67. https://doi.org/10.3897/zookeys.223.3540.

Wang, Xingeng, and Ellen M. Aparicio. (2020) ‘Reproductive Traits of Ontsira Mellipes (Hymenoptera: Braconidae), a North American Parasitoid, as a Novel Biological Control Agent for Exotic Anoplophora Glabripennis (Coleoptera: Cerambycidae)’. Journal of Economic Entomology 113, no. 5: 2112–19. https://doi.org/10.1093/jee/toaa160.

Zaldívar-Riverón, Alejandro, Mark R Shaw, Alberto G Sáez, et al. ‘Evolution of the Parasitic Wasp Subfamily Rogadinae (Braconidae): Phylogeny and Evolution of Lepidopteran Host Ranges and Mummy Characteristics’. Invertebrate Systematics 8 (2008): 329. https://doi.org/10.1186/1471-2148-8-329.

