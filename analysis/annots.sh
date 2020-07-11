#!/bin/bash

#SBATCH -p main
#SBATCH -J new_annots
#SBATCH -t 23:00:00
#SBATCH -c 8
#SBATCH --mem=16G

module load R/3.6.1
module load python-3.6.0
module load squashfs/4.3-lt2t
module load singularity/3.5.3

#Rscript annots.R

rm -rf txrevise
git clone git@github.com:andreasvija/txrevise.git
cd txrevise
git checkout CAGE

cp ../new_transcripts_25.rds data/
cd scripts
mkdir processed
wget -O processed/Homo_sapiens.GRCh38.96.gtf.gz ftp://ftp.ensembl.org/pub/release-96/gtf/homo_sapiens/Homo_sapiens.GRCh38.96.gtf.gz

mkdir SlurmOut
snakemake -s cage.Snakefile -p processed/Homo_sapiens.GRCh38.96_log.txt --cluster ./snakemake_submit_UT.py --jobs 100 --config fill=TRUE start_end_diff=25 cage=../data/CAGE_promoter_annotations_25.rds --use-singularity
