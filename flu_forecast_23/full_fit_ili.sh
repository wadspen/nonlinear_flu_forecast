#!/bin/bash

#SBATCH --time=5-23:15:00
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

Rscript ./full_ili_fits.R "$1" "$2"
