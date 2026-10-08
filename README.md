# paper-mg-insects

Repository for data, scripts and figures related to a paper analyzing the size and composition of the MG insect fauna based on data from the Insect Biome Atlas project.

The subdirectories contain the following:
 - `alpha_beta`: Data, scripts and figures for analyses of alpha and beta diversity
 - `clean`: Data and scripts for analyzing contaminations in MG litter data
 - `composition`: Data, scripts and figures for analyses of ecological and taxonomic composition
 - `diversification_tests`: Data, scripts and figures for diversification analyses
 - `family_composition`: Data, scripts and figures for analyses of taxonomic composition at the family and family-clade levels
 - `fig_settings`: Shared settings (colours) used for all figures
 - `iba_data`: Scripts for processing original IBA data so that it is suitable for the analyses in the study
 - `mg_maps`: Data, scripts and figures for mapping spatial features onto the MG map
 - `ms_figures`: Script for assembling and numbering the figures used in the manuscript
 - `parasitoid_load`: Data, scripts and figures for analyzing parasitoid load
 - `placements`: Data, scripts and results of analyzing phylogenetic placement data
 - `ref_tree`: Data, scripts and output from generating the expanded and corrected reference tree used in phylogenetic placement
 - `sampling_bias`: Data, scripts and figures from the analysis of the likely sampling bias with respect to the true insect fauna
 - `source_data`: Source data for life history traits and for the updated classifications used in expanding the original Chesters tree for some taxa
 - `specpool_accumulation`: Data, scripts and figures for species accumulation plots
 - `traits`: Scripts used in processing and assembling life-history trait information
 - `validation`: Contains validation analyses for the pipeline used to process and annotate the metabarcoding data 

Before rerunning or modifying any of the analyses, it is necessary to download the raw data and processed data from the IBA project from the relevant figshare repositoires:
 - [Processed data](https://doi.org/10.17044/scilifelab.27202368) (version 4 was used for the analyses in the paper)
 - [Raw data and sample metadata](https://doi.org/10.17044/scilifelab.25480681) (version 6 was used for the analyses in the paper)

After this, the processing scripts in the `iba_data` directory should be run, as many of the other scripts are dependent on the generated output files, which are not in the repository.

Appropriate dumps from GBIF and from NCBI (the taxonomy) are also needed in some cases. Further information about this and other relevant instructions are found in the README files found within each of the subdirectories.

