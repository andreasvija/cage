#!/bin/bash

#SBATCH -p main
#SBATCH -J filter_common_variants
#SBATCH -t 08:00:00
#SBATCH -c 8
#SBATCH --mem=8G

module load bcftools-1.9

bcftools isec -n~11 -w1 cage_variants_common.vcf.gz -w2 txrevise_variants_common.vcf.gz /gpfs/hpc/projects/GENCORD/Garieri_2017/Garieri_filtered.no_DS.vcf.gz /gpfs/hpc/projects/genomic_references/GEUVADIS/genotypes/GEUVADIS_GRCh38_filtered.vcf.gz
