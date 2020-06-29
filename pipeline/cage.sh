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

git clone https://github.com/eQTL-Catalogue/qtlmap.git

~/cage/nextflow/nextflow \
	run qtlmap/main.nf \
	-profile tartu_hpc \
	--studyFile sources_cage.tsv \
	--cis_window 200000 \
	--run_permutation true \
	--varid_rsid_map_file /gpfs/hpc/projects/genomic_references/annotations/eQTLCatalogue/v0.1/dbSNP_b151_GRCh38p7_splitted_var_rsid.vcf.gz \
	-resume
