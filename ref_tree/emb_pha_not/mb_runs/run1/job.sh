#!/bin/bash
#
#SBATCH -J emphno1
#SBATCH -t 10:00:00
#SBATCH -n 64
#
mpprun ~/MrBayes/src/mb run.nex > log.txt

