#!/bin/bash

#SBATCH -p main
#SBATCH -J run_visualizations
#SBATCH -t 3-00:00:00
#SBATCH -c 8
#SBATCH --mem=16G

module load R/3.6.1
mkdir plots

for KIND in cage_over_tx tx_over_cage ann_over_cage
do
  rm -rf plots/$KIND
  mkdir plots/$KIND
  Rscript visualization.R --kind $KIND
done
