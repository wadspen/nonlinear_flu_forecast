#!/bin/bash

#SBATCH	--begin=now+2hour
#SBATCH --time=6:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=36
#SBATCH --mem=0
#SBATCH --exclusive

module load gcc
module load r
module load udunits
module load r-rgdal
module load proj

Rscript ./sim_eval.R "$1"
