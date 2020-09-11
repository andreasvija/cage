#!/bin/bash

#SBATCH -p main
#SBATCH -J get_variant_info
#SBATCH -t 24:00:00
#SBATCH -c 8
#SBATCH --mem=32G

gunzip -c ../qtlmap_prep/cage_common.vcf > temp.vcf
gunzip -c ../qtlmap_prep/txrevise_common.vcf > temp_geuvadis.vcf

rm variantinfo.vcf
rm variantinfo_geuvadis.vcf

while read VARIANT
do
  cat temp.vcf | grep $VARIANT >> variantinfo.vcf
  cat temp_geuvadis.vcf | grep $VARIANT >> variantinfo_geuvadis.vcf
done < variants.txt

rm temp.vcf
rm temp_geuvadis.vcf
