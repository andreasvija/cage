#!/bin/bash

#SBATCH -p main
#SBATCH -J run_visualizations
#SBATCH -t 2-00:00:00
#SBATCH -c 8
#SBATCH --mem=16G

module load R/3.6.1
mkdir plots
Rscript visualization.R
