#!/bin/bash

#SBATCH -p main
#SBATCH -J integrate_annots
#SBATCH -t 23:00:00
#SBATCH -c 4
#SBATCH --mem=8G

ANNOTATION=Homo_sapiens.GRCh38.96
IFS='.' read -r -a array <<< "$ANNOTATION"
RELEASE_N=${array[2]}

module load R/3.6.1
module load squashfs/4.3-lt2t
module load singularity/3.5.3

git clone git@github.com:andreasvija/txrevise.git
cd txrevise
git checkout CAGE
git pull

cd scripts
mkdir processed
mkdir processed/input

if [ ! -f processed/input/promoters.tsv ]; then
    cp ../../../qtlmap_prep/FANTOM5_promoter_annotations.tsv processed/input/promoters.tsv
fi
if [ ! -f processed/input/${ANNOTATION}.gtf.gz ]; then
    wget -O processed/input/${ANNOTATION}.gtf.gz ftp://ftp.ensembl.org/pub/release-${RELEASE_N}/gtf/homo_sapiens/${ANNOTATION}.gtf.gz
fi
if [ ! -f processed/input/gene_metadata.tsv.gz ]; then
    wget -O processed/input/gene_metadata.tsv.gz https://zenodo.org/record/3366011/files/gene_counts_Ensembl_96_phenotype_metadata.tsv.gz
fi
if [ ! -f txrevise.img ]; then
    singularity build txrevise.img docker://andreasvija/txrevise:latest
fi

mkdir SlurmOut
source ~/.bashrc
conda activate snakemake
OUTPUT=processed/${ANNOTATION}_all_completed.txt
rm $OUTPUT
snakemake -ps cage.Snakefile $OUTPUT --cluster ./snakemake_submit_UT.py --jobs 200 --config fill=TRUE Ns="5 10 15 20 25 50 100 200" --use-singularity
conda deactivate
