#!/bin/bash

#SBATCH --time=5-23:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=36
#SBATCH --mem=0
#SBATCH --exclusive

module load gcc
module load r
module load udunits
module load proj
module load r-rgdal

Rscript ./ili_eval.R "$1"
