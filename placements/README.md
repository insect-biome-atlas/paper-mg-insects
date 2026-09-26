# Placement info

This directory contains information for enerating placement info, placement ages,
useful in analyzing MG radiations

## `source` directory
Contains source data:
 - `gbif_mg_species_list.tsv`: Species list for administrative area Madagascar from GBIF. Citation information:
    GBIF.org (20 April 2026) GBIF Occurrence Download  https://doi.org/10.15468/dl.fa4xsc. The list confirms that
    none of the species included in the reference tree is known from Madagascar, validating its use in the analysis
    of MG radiations.
 - `mg_cluster_placements.tsv`: summary of the epa-ng placements of MG clusters. `cluster_rep` is the representative ASV for the cluster.
 -	`mg_placement_tree_edgenum.nwk`: the reference tree with the indices used by epa-ng + gappa in the placement files as tip and node labels.
    The topology is from the expanded Chesters tree but branch lengths
    have been reoptimized by RaXML for the placement analysis.
    Note that the root has been removed by collapsing the branch leading from the Wolbachia sequence to the rest of the tree.
 -	`mg_placement_tree_tiplabels.nwk`: the same tree with original tip labels.

## `R` directory
Contains the following R scripts:
 - `compute_placement_ages.R`: Uses information in `source` and generates output with placement ages and other info in `data`.
 - `compute_gbif_barcode_overlap.R`: Compute the overlap between IBA MG species and GBIF species recorded from Madagascar / barcoded species
