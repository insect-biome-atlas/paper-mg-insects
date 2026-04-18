# iba_data
This repo contains scripts to assemble and process IBA data for downstream analyses of
faunal composition.

To generate all data files, simply run `batch_script.R` in the `R` subdirectory. This will
invoke the relevant R scripts in the correct order and generate all tsv and rds output
files, which will be placed in the parent directory. First make sure that the paths to
the original IBA data files are set correctly (see TODO lines in the R scripts).

Note that the generated taxonomy files are filtered on fairly stringent annotation criteria.
Details on the content of each generated data file can be found in the R scripts generating
them.

To keep the github repo small, the output data files are not committed to the repo.

