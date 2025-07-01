#!/bin/bash

#SBATCH --time=5-23:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --mem=0


module load gcc
module load r
module load udunits
module load proj
module load r-rgdal
#Rscript ./comp_fit_hosp.R "$1"
Rscript ./full_fit_hosp.R "$1" "$2" 
