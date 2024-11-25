#!/bin/bash

#SBATCH --begin=now+12hour
#SBATCH --time=23:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --mem=0
#SBATCH --exclusive

module load gcc
module load r
module load udunits
module load r-rgdal
module load proj

Rscript ./sim_score.R "$1"
