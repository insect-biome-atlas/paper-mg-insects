#!/bin/bash
#
#SBATCH -J zorap1
#SBATCH -t 5:00:00
#SBATCH -n 64
#
mpprun ~/MrBayes/src/mb run.nex > log.txt

