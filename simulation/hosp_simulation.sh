#!/bin/bash

#SBATCH --time=12-23:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --constraint=intel
module load gcc
module load r
module load udunits
module load r-rgdal
module load proj

Rscript ./full_fit_hosp_sim.R "$1" "$2"