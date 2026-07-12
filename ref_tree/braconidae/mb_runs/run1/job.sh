#!/bin/bash
#
#SBATCH -J brac1
#SBATCH -t 168:00:00
#SBATCH -n 64
#
mpprun ~/MrBayes/src/mb run.nex > log.txt

