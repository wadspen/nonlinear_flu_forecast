#!/bin/bash

#SBATCH --time=3:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=0
#SBATCH --exclusive

module load gcc
module load r
module load udunits
module load r-rgdal
module load proj

Rscript ./sim_sample_size_mse.R "$1"
