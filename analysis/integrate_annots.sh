#!/bin/bash

#SBATCH -p main
#SBATCH -J integrate_annots
#SBATCH -t 23:00:00
#SBATCH -c 8
#SBATCH --mem=24G

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
cp ../../new_transcripts_25.rds processed/input/
wget -O processed/input/${ANNOTATION}.gtf.gz ftp://ftp.ensembl.org/pub/release-${RELEASE_N}/gtf/homo_sapiens/${ANNOTATION}.gtf.gz

if [ ! -f txrevise.img ]; then
    singularity build txrevise.img docker://kauralasoo/txrevise:latest
fi

mkdir SlurmOut
source ~/.bashrc
conda activate snakemake
OUTPUT=processed/${ANNOTATION}_all_completed.txt
rm $OUTPUT
snakemake -pns cage.Snakefile $OUTPUT --cluster ./snakemake_submit_UT.py --jobs 200 --config fill=TRUE Ns="25" --use-singularity
conda deactivate
