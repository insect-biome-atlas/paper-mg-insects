# Script for combining cruaud taxonomy data with gb_accn annotations
# with extensice sequence metadata

D <- read.delim("cruaud_taxonomy_coi_min_bc_cov_300.tsv")
E <- read.delim("cruaud_ncbi_coi_survey_steps1-3.tsv")

# Remove a couple of matches that are untenable
# These are matched by NCBI to "Encarsia" and "Ormyrus", respectively,
# and there are several tips in each of these genera already.
D$gb_accn[match("Dirphys",D$Genus)] <- ""
D$gb_accn[match("Ormyrulus",D$Genus)] <- ""

D <- merge(D, E, all.x=TRUE)

D$aln_chalc_bc_start <- D$aln_t_start
D$aln_chalc_bc_end <- D$aln_t_end
D$aln_chalc_bc_cov <- D$aln_t_end

idx <- which(!is.na(D$aln_query))
D$aln_chalc_bc_end[idx] <- 412 - (652 - D$aln_t_end[idx])
for (i in 1:length(idx)) { D$aln_chalc_bc_start[idx[i]] <- max(1, D$aln_t_start[idx[i]] - 240) }
D$aln_chalc_bc_cov[idx] <- D$aln_chalc_bc_end[idx] - D$aln_chalc_bc_start[idx] + 1

write.table(D,"cruaud_taxonomy_coi_min_bc_cov_300_seq_meta.tsv", sep="\t", row.names=FALSE)

