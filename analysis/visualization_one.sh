#!/bin/bash

#SBATCH -p main
#SBATCH -J run_visualization_one
#SBATCH -t 01:00:00
#SBATCH -c 8
#SBATCH --mem=16G

module load R/3.6.1
mkdir plots_one

Rscript visualization_one.R
