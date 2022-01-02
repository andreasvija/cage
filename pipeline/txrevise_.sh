#!/bin/bash

#SBATCH -p main
#SBATCH -J txrevise_pipeline_
#SBATCH -t 24:00:00
#SBATCH -c 2
#SBATCH --mem=8G

module load java-1.8.0_40
module load singularity/3.5.3
module load squashfs/4.4
module load nextflow

mkdir txrevise
cd txrevise

git clone https://github.com/eQTL-Catalogue/rnaseq
cd rnaseq
find . -type f -exec sed -i -e 's#hpc/projects/#space/projects/#g' {} \;
cd ..
git clone https://github.com/eQTL-Catalogue/qcnorm
cd qcnorm
git checkout d90d2655d5b147e34157caf8544a9445f2748ef2
find . -type f -exec sed -i -e 's#hpc/projects/#space/projects/#g' {} \;
cd ..

set -e

nextflow \
	run rnaseq/main.nf \
	-profile eqtl_catalogue \
	--readPathsFile ../GEUVADIS_readPaths.tsv \
	--unstranded \
	--skip_qc \
	--skip_multiqc \
	--skip_stringtie \
	--run_txrevise \
	--txrevise_gffs '../../analysis/txrevise/scripts/processed/Homo_sapiens.GRCh38.96_regular/txrevise_*.gff3' \
	-resume \
	-executor.queueSize 100

nextflow \
	run qcnorm/normalisation.nf \
	-profile tartu_hpc \
	-resume \
	--study_name txrevise \
	--quant_results_path results \
	--sample_meta_path ../../analysis/GEUVADIS_EUR.tsv \
	--txrev_pheno_meta_path ../../analysis/txrevise/scripts/processed/Homo_sapiens.GRCh38.96_regular/txrevise_regular_phenotype_metadata.tsv.gz \
	--skip_exon_norm \
	--skip_tx_norm \
	--skip_leafcutter_norm \
	--outdir qcnorm_out
