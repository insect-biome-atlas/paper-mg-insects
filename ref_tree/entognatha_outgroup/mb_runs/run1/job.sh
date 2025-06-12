#!/bin/bash
#
#SBATCH -J entorun1
#SBATCH -t 24:00:00
#SBATCH -n 64
#
mpprun ~/MrBayes/src/mb run.nex > log.txt

