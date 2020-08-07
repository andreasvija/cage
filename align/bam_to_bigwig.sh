#!/bin/bash
# screen, ssh stage1

#SBATCH -p amd
#SBATCH -J bam_to_bigwig
#SBATCH -t 24:00:00
#SBATCH -c 2
#SBATCH --mem=4G

module load python-3.6.0
mkdir SlurmOut
snakemake -s bam_to_bigwig.Snakefile --cluster ./snakemake_submit_UT.py -p out_bigwig.txt --configfile config.yaml --jobs 20
