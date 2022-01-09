# Improved detection of promoter usage QTLs with augmented transcript annotations

## Order of running code

### Alignment
1. align/bwa.sh
1. aligns/bam_to_bigwig.sh

### qtlmap prep
1. qtlmap_prep/clean_promoters.R
1. qtlmap_prep/promoter_annotations_to_gtf.py
1. qtlmap_prep/cagefc.sh
1. qtlmap_prep/convert.Rmd
1. qtlmap_prep/filter_common_variants.sh

### Pipelines
1. pipelines/cage.sh
1. pipelines/txrevise.sh
1. pipelines/annots_20.sh
1. pipelines/reset.sh

### Analysis
1. analysis/integrate_annots.sh
1. analysis/visualize_cage_annotations.R
1. analysis/get_variant_info.sh
1. analysis/analysis.Rmd
1. analysis/visualization_one.sh
