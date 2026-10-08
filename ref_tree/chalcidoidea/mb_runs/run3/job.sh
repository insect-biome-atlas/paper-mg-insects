#!/bin/bash
#
#SBATCH -J charun1
#SBATCH -t 120:00:00
#SBATCH -n 64
#
mpprun ~/MrBayes/src/mb run.nex > log.txt

