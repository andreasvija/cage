#!/bin/bash

#SBATCH -p main
#SBATCH -J create_annots
#SBATCH -t 23:00:00
#SBATCH -c 8
#SBATCH --mem=16G

module load R/3.6.1

Rscript annots.R
