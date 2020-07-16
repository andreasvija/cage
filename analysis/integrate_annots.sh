#!/bin/bash

#SBATCH -p main
#SBATCH -J integrate_annots
#SBATCH -t 23:00:00
#SBATCH -c 1
#SBATCH --mem=4G

module load R/3.6.1
module load squashfs/4.3-lt2t
module load singularity/3.5.3

Rscript makeTxreviseCageMetadata.R

git clone git@github.com:andreasvija/txrevise.git
cd txrevise
git checkout CAGE
git pull

cd scripts
cp ../../new_transcripts_25.rds processed/
wget -O processed/Homo_sapiens.GRCh38.96.gtf.gz ftp://ftp.ensembl.org/pub/release-96/gtf/homo_sapiens/Homo_sapiens.GRCh38.96.gtf.gz

if [ ! -f txrevise.img ]; then
    singularity build txrevise.img docker://kauralasoo/txrevise:latest
fi

mkdir SlurmOut
source ~/.bashrc
conda activate snakemake
snakemake -ps cage.Snakefile processed/Homo_sapiens.GRCh38.96_log.txt --cluster ./snakemake_submit_UT.py --jobs 100 --config fill=TRUE start_end_diff=25 --use-singularity
conda deactivate
