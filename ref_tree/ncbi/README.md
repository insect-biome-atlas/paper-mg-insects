# ncbi

This directory contains an R script (`extract_arthropod_genus_classification.R`) that extracts information on the family, subfamily, supertribe and tribe placements of all arthropod genera in the NCBI taxonomy. The script was run on a dump of the NCBI taxonomy from 2025-04-17. The NCBI taxonomy dump files are not committed to the repo as they are too large.

The script generates a table for use in obtaining detailed NCBI classification of arthropod genera, `ncbi_arhtropod_genus_classification_detailed.tsv`.

To regenerate the table, first download a local dump of the NCBI taxonomy and unpack it. Then run the R script (uses the dump files `nodes.dmp` and `rankedlineage.dmp`).

The script also generates a more extensive table for use in certain downstream scripts, `ncbi_name_ranked_classification.tsv`. A compressed version of this file is provided in the repo.

