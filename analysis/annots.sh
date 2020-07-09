#!/bin/bash

#SBATCH -p main
#SBATCH -J new_annots
#SBATCH -t 23:00:00
#SBATCH -c 4
#SBATCH --mem=16G

module load R/3.6.1
module load python-3.6.0

Rscript annots.R

git clone git@github.com:andreasvija/txrevise.git
cd txrevise
git checkout CAGE

cp ../new_transcripts_25.rds data/
cd scripts
mkdir processed
wget -O processed/Homo_sapiens.GRCh38.96.gtf.gz ftp://ftp.ensembl.org/pub/release-96/gtf/homo_sapiens/Homo_sapiens.GRCh38.96.gtf.gz

conda activate
annotation="Homo_sapiens.GRCh38.96"
python extractTranscriptTags.py --gtf processed/"$annotation".gtf.gz > processed/"$annotation".transcript_tags.txt
conda deactivate

Rscript prepareAnnotations.R --gtf processed/"$annotation".gtf.gz --tags processed/"$annotation".transcript_tags.txt --out processed/"$annotation".txrevise_annotations.rds

cd ..
Rscript data/prepare_CAGE_data.R
cd scripts

rm -rf SlurmOut
mkdir SlurmOut
snakemake -p processed/Homo_sapiens.GRCh38.96_log.txt --cluster ./snakemake_submit_UT.py --jobs 100 --config fill=TRUE start_end_diff=25 cage=../data/CAGE_promoter_annotations_25.rds --use-singularity
