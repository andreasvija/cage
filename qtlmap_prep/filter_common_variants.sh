#!/bin/bash

#SBATCH -p main
#SBATCH -J filter_common_variants
#SBATCH -t 01:00:00
#SBATCH -c 8
#SBATCH --mem=32G

module load bcftools-1.9

cp /gpfs/hpc/projects/GENCORD/Garieri_2017/Garieri_filtered.no_DS.vcf.gz cage.vcf.gz
cp /gpfs/hpc/projects/genomic_references/GEUVADIS/genotypes/GEUVADIS_GRCh38_filtered.vcf.gz txrevise.vcf.gz

bcftools index cage.vcf.gz
bcftools index txrevise.vcf.gz

bcftools isec -n~11 -p temp cage.vcf.gz txrevise.vcf.gz

rm cage.vcf.gz
rm cage.vcf.gz.csi
rm txrevise.vcf.gz
rm txrevise.vcf.gz.csi

cp temp/0000.vcf cage_common.vcf
gzip cage_common.vcf
cp temp/0001.vcf txrevise_common.vcf
gzip txrevise_common.vcf

rm -rf temp
