#!/bin/bash

#SBATCH -p main
#SBATCH -J bwa
#SBATCH -t 24:00:00
#SBATCH -c 8
#SBATCH --mem=32G

module load python-3.6.0
mkdir SlurmOut
snakemake -s bwa.Snakefile --cluster ./snakemake_submit_UT.py -p out.txt --configfile config.yaml --jobs 8
