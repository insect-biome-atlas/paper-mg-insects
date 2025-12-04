# Code and scripts for cleaning chesters expanded sequence data

This directory contains data and code for filtering out noise from the
expanded Chesters data. All sequences in the file `chesters_expanded_raw.fasta`
were run with blastn against the NCBI 'nucleotide' database; the results are in
`blastn.tsv` for the top 100 hits.

The file `lca_fxns.R` contain functions useful for computing lca (least commoon
ancestor) taxonomy information for sets of blastn hits.

The script `computeblastn_lca.R` applies a percent identity threshold and then
computes the least common ancestor (`lca`) in the ranked classification in the
NCBI taxonomy (using data in the `ncbi` directory and a function in `lca_fxns.R`)

The script `analyze_lca_res.R` finds likely errors by removing the top hit (the
hit to the original record) and then looks for mismatches at the class or 
order level. The files `class_miss_stats.tsv` and `order_miss_stats.tsv`
contain statistics generated in these searches. The detailed taxonomy information
for the hits of the relevant taxa (tip labels in the tree)  are in the folders
`class_miss_taxonomies` and `order_miss_taxonomies`. Note that there is overlap
between these sets, in that all class misses are also necessarily order misses.

Trying to find errors by looking for inconsistencies at the family level was
not deemed meaningful, in that it is quite difficult to separate errors from
real signal at this level.

All potential errors (<95% of the top 100 hits assigned to the wrong order or
class) were checked manually before being definitely identified as errors. The
final set of identified errors is found in the `remove_chesters_errors.R` script
in the `chesters_2017` directory.

Note that there are some order-level misses in the Collembola data, which
were assumed to be due to the relative sparsity of Collembola records in NCBI.

