#!/bin/bash

#SBATCH -p main
#SBATCH -J make_cage_txrevise_metadata
#SBATCH -t 2-00:00:00
#SBATCH -c 8
#SBATCH --mem=16G

module load R/3.6.1
Rscript makeTxreviseCageMetadata.R
