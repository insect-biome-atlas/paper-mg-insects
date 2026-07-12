#!/bin/bash

# TODO: Ensure correct tpplc path on your machine (or compile through python or R package)
tpplc tdbd.tppl --method mcmc --align --kernel --output tdbd_drift_mcmc

echo "run1"
./tdbd_drift_mcmc ../data/div_data.json --iterations 125000 --sampling-period 100 --debug-mcmc --mcmc-global-prob 0.01 > tdbd_lez_125000-100_1.csv 2> tdbd_lez_125000-100_1_logLike.csv
./format_tppl_samples.sh tdbd_lez_125000-100_1.csv
./format_tppl_logLikes.sh tdbd_lez_125000-100_1_logLike.csv

echo "run2"
./tdbd_drift_mcmc ../data/div_data.json --iterations 125000 --sampling-period 100 --debug-mcmc --mcmc-global-prob 0.01 > tdbd_lez_125000-100_2.csv 2> tdbd_lez_125000-100_2_logLike.csv
./format_tppl_samples.sh tdbd_lez_125000-100_2.csv
./format_tppl_logLikes.sh tdbd_lez_125000-100_2_logLike.csv

echo "run3"
./tdbd_drift_mcmc ../data/div_data.json --iterations 125000 --sampling-period 100 --debug-mcmc --mcmc-global-prob 0.01 > tdbd_lez_125000-100_3.csv 2> tdbd_lez_125000-100_3_logLike.csv
./format_tppl_samples.sh tdbd_lez_125000-100_3.csv
./format_tppl_logLikes.sh tdbd_lez_125000-100_3_logLike.csv

echo "run4"
./tdbd_drift_mcmc ../data/div_data.json --iterations 125000 --sampling-period 100 --debug-mcmc --mcmc-global-prob 0.01 > tdbd_lez_125000-100_4.csv 2> tdbd_lez_125000-100_4_logLike.csv
./format_tppl_samples.sh tdbd_lez_125000-100_4.csv
./format_tppl_logLikes.sh tdbd_lez_125000-100_4_logLike.csv
