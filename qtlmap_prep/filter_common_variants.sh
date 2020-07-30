#!/bin/bash

#SBATCH -p main
#SBATCH -J filter_common_variants
#SBATCH -t 08:00:00
#SBATCH -c 4
#SBATCH --mem=8G

module load bcftools-1.9

bcftools query -f '%ID' -o variants /gpfs/hpc/projects/GENCORD/Garieri_2017/Garieri_filtered.no_DS.vcf.gz
bcftools query -f '%ID' -o variants_geuvadis /gpfs/hpc/projects/genomic_references/GEUVADIS/genotypes/GEUVADIS_GRCh38_filtered.vcf.gz
