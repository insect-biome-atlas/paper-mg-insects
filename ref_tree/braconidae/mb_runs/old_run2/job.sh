#!/bin/bash
#
#SBATCH -J bracrun2
#SBATCH -t 168:00:00
#SBATCH -n 64
#
mpprun ~/MrBayes/src/mb run.nex > log.txt

