#!/bin/bash
#screen, ssh stage1

#SBATCH -p main
#SBATCH -J txrevise_pipeline
#SBATCH -t 24:00:00
#SBATCH -c 8
#SBATCH --mem=24G

module load java-1.8.0_40
module load singularity
module load nextflow

mkdir txrevise
cd txrevise

git clone https://github.com/eQTL-Catalogue/rnaseq
git clone https://github.com/eQTL-Catalogue/qcnorm
git clone https://github.com/eQTL-Catalogue/qtlmap.git

NXF_VER=18.10.1 nextflow \
	run rnaseq/main.nf \
	-profile eqtl_catalogue \
	--readPathsFile /gpfs/hpc/projects/eQTLCatalogue/SampleArcheology/readPaths/GEUVADIS_readPaths.tsv \
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

../../nextflow/nextflow \
	run qtlmap/main.nf \
	-profile tartu_hpc \
	--studyFile ../sources_txrevise.tsv \
	--cis_window 200000 \
	--run_permutation true \
	--is_imputed FALSE \
	--varid_rsid_map_file /gpfs/hpc/projects/genomic_references/annotations/eQTLCatalogue/v0.1/dbSNP_b151_GRCh38p7_splitted_var_rsid.vcf.gz \
	-resume

gunzip results/sumstats/txrevise.permuted.txt.gz -c > ../../analysis/txrev-25.permuted.txt
cp results/featureCounts/merged_gene_counts.txt ../../analysis/txrev-25_merged_gene_counts.txt
