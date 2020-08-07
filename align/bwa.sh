#!/bin/bash
# screen, ssh stage1

#SBATCH -p main
#SBATCH -J bwa
#SBATCH -t 24:00:00
#SBATCH -c 2
#SBATCH --mem=4G

module load python-3.6.0
mkdir SlurmOut
snakemake -s bwa.Snakefile --cluster ./snakemake_submit_UT.py -p out.txt --configfile config.yaml --jobs 20
