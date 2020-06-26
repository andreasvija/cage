#!/bin/bash
#screen, ssh stage1

#SBATCH -p main
#SBATCH -J cage_pipeline
#SBATCH -t 24:00:00
#SBATCH -c 8
#SBATCH --mem=24G

module load java-1.8.0_40
module load singularity
module load nextflow

git clone https://github.com/kerimoff/qtlmap.git

~/cage/nextflow/nextflow \
	run ~/cage/pipeline/qtlmap/main_multi_study.nf \
	-profile tartu_hpc \
	--studyFile ~/cage/qtlmap_prep/sources.tsv \
	--cis_window 200000 \
	--run_permutation true \
	-resume
