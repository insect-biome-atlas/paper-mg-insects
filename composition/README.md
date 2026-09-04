# composition_accumulation

This repo contains scripts used to plot accumulation curves for faunal composition parameters

The R scripts are in the directory `R` and the output figures in the directory `figs`.

The main R scripts are as follows:
- `batch_script.R`: script for running all R scripts in the directory in order
- `prepare_data.R`: prepares data and puts it into the `data` directory.
- `filter_se_traps.R`: help function to filter SE traps based on number of samples
- `compute_site_accumulation_niche.R`: computes site accumulation of niche composition
- `compute_site_accumulation_habitat.R`: computes site accumulation of habitat composition
- `compute_site_accumulation_taxonomic.R`: computes site accumulation of taxonomic composition
- `compute_sample_accumulation_niche.R`: computes sample accumulation of niche composition
- `compute_sample_accumulation_habitat.R`: computes sample accumulation of habitat composition
- `compute_sample_accumulation_taxonomic.R`: computes sample accumulation of taxonomic composition
- `compute_parasitoid_subgroup_ratios.R`: computes Madagascan to Swedish forest species ratios of various parasitoid groups
- `plot_site_accumulation_niche.R`: plot site accumulation of niche composition
- `plot_site_accumulation_habitat.R`: plot site accumulation of habitat composition
- `plot_site_accumulation_taxonomic.R`: plot site accumulation of taxonomic composition
- `plot_sample_accumulation_niche.R`: plot sample accumulation of niche composition
- `plot_sample_accumulation_habitat.R`: plot sample accumulation of habitat composition
- `plot_sample_accumulation_taxonomic.R`: plot sample accumulation of taxonomic composition
- `plot_site_composition_data.R`: plot site composition data per forest type

See the scripts for information on the generated figures. The figure names generally match the name of the generating script.
