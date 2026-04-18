# Assemble counts data needed for analyses
# NB!! Data assembled here have original taxonomic annotations (sintax+)

# TODO: Set the path to the iba utils repo on your system
utils_path <- "~/dev/ms-repos-iba/utils/"
source(paste0(utils_path,"R/spikes_controls_fxns.R"))
source(paste0(utils_path,"get_iba_co1_data_fxn.R"))

# TODO: Set the path to the iba processed data and metadata repos
# on your system
data_path <- "~/dev/figshare-repos/iba/processed_data/v4-prel/"
metadata_path <- "~/dev/figshare-repos/iba/raw_data/v6/"
malaise_mg <- get_iba_co1_data(data_path=data_path,
                               metadata_path=metadata_path,
                               country="MG",
                               dataset="lysate",
                               calibrate=TRUE)
malaise_se <- get_iba_co1_data(data_path=data_path,
                               metadata_path=metadata_path,
                               country="SE",
                               dataset="lysate",
                               calibrate=TRUE)

litter_mg <- get_iba_co1_data(data_path=data_path,
                              metadata_path=metadata_path,
                              country="MG",
                              dataset="litter",
                              remove_spikes=TRUE)
litter_se <- get_iba_co1_data(data_path=data_path,
                              metadata_path=metadata_path,
                              country="SE",
                              dataset="litter",
                              remove_spikes=TRUE)

malaise_spec_mg <- get_iba_co1_data(data_path=data_path,
                                    metadata_path=metadata_path,
                                    country="MG",
                                    dataset="lysate",
                                    calibrate=TRUE,
                                    calibration_type="specimens")

malaise_spec_se <- get_iba_co1_data(data_path=data_path,
                                    metadata_path=metadata_path,
                                    country="SE",
                                    dataset="lysate",
                                    calibrate=TRUE,
                                    calibration_type="specimens")

malaise_litter_mg <- get_iba_co1_data(data_path=data_path,
                                    metadata_path=metadata_path,
                                    country="MG",
                                    dataset="lysate|litter",
                                    calibrate=FALSE)

malaise_litter_se <- get_iba_co1_data(data_path=data_path,
                                    metadata_path=metadata_path,
                                    country="SE",
                                    dataset="lysate|litter",
                                    calibrate=FALSE)

# Write data tables
write.tsv <- function(D,file) { write.table(D,file,row.names=FALSE,sep="\t") }
write.tsv(malaise_mg,"raw_counts_malaise_mg.tsv")
write.tsv(malaise_se,"raw_counts_malaise_se.tsv")
write.tsv(litter_mg,"raw_counts_litter_mg.tsv")
write.tsv(litter_se,"raw_counts_litter_se.tsv")
write.tsv(malaise_spec_mg,"raw_counts_malaise_spec_mg.tsv")
write.tsv(malaise_spec_se,"raw_counts_malaise_spec_se.tsv")
write.tsv(malaise_litter_mg,"raw_counts_uncal_malaise_litter_mg.tsv")
write.tsv(malaise_litter_se,"raw_counts_uncal_malaise_litter_se.tsv")

