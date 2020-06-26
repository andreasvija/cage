#!/bin/bash

#SBATCH -p main
#SBATCH -J cage_fc
#SBATCH -t 07:00:00
#SBATCH -c 4
#SBATCH --mem=8G

/gpfs/hpc/projects/genomic_references/software/bin/featureCounts -t promotor -f -O -s 0 -T 4 -a FANTOM5_promoter_annotations.gtf -o bwa.counts ../align/results/ERR2021349.sorted.bam ../align/results/ERR2021350.sorted.bam ../align/results/ERR2021351.sorted.bam ../align/results/ERR2021352.sorted.bam ../align/results/ERR2021353.sorted.bam ../align/results/ERR2021354.sorted.bam ../align/results/ERR2021355.sorted.bam ../align/results/ERR2021356.sorted.bam ../align/results/ERR2021357.sorted.bam ../align/results/ERR2021358.sorted.bam ../align/results/ERR2021359.sorted.bam ../align/results/ERR2021360.sorted.bam ../align/results/ERR2021361.sorted.bam ../align/results/ERR2021362.sorted.bam ../align/results/ERR2021363.sorted.bam ../align/results/ERR2021364.sorted.bam ../align/results/ERR2021365.sorted.bam ../align/results/ERR2021366.sorted.bam ../align/results/ERR2021367.sorted.bam ../align/results/ERR2021368.sorted.bam ../align/results/ERR2021369.sorted.bam ../align/results/ERR2021370.sorted.bam ../align/results/ERR2021371.sorted.bam ../align/results/ERR2021372.sorted.bam ../align/results/ERR2021373.sorted.bam ../align/results/ERR2021374.sorted.bam ../align/results/ERR2021375.sorted.bam ../align/results/ERR2021376.sorted.bam ../align/results/ERR2021377.sorted.bam ../align/results/ERR2021378.sorted.bam ../align/results/ERR2021379.sorted.bam ../align/results/ERR2021380.sorted.bam ../align/results/ERR2021381.sorted.bam ../align/results/ERR2021382.sorted.bam ../align/results/ERR2021383.sorted.bam ../align/results/ERR2021384.sorted.bam ../align/results/ERR2021385.sorted.bam ../align/results/ERR2021386.sorted.bam ../align/results/ERR2021387.sorted.bam ../align/results/ERR2021388.sorted.bam ../align/results/ERR2021389.sorted.bam ../align/results/ERR2021390.sorted.bam ../align/results/ERR2021391.sorted.bam ../align/results/ERR2021392.sorted.bam ../align/results/ERR2021393.sorted.bam ../align/results/ERR2021394.sorted.bam ../align/results/ERR2021395.sorted.bam ../align/results/ERR2021396.sorted.bam ../align/results/ERR2021397.sorted.bam ../align/results/ERR2021398.sorted.bam ../align/results/ERR2021399.sorted.bam ../align/results/ERR2021400.sorted.bam ../align/results/ERR2021401.sorted.bam ../align/results/ERR2021402.sorted.bam ../align/results/ERR2021403.sorted.bam ../align/results/ERR2021404.sorted.bam ../align/results/ERR2021405.sorted.bam ../align/results/ERR2021406.sorted.bam ../align/results/ERR2021407.sorted.bam ../align/results/ERR2021408.sorted.bam ../align/results/ERR2021409.sorted.bam ../align/results/ERR2021410.sorted.bam ../align/results/ERR2021411.sorted.bam ../align/results/ERR2021412.sorted.bam ../align/results/ERR2021413.sorted.bam ../align/results/ERR2021414.sorted.bam ../align/results/ERR2021415.sorted.bam ../align/results/ERR2021416.sorted.bam ../align/results/ERR2021417.sorted.bam ../align/results/ERR2021418.sorted.bam ../align/results/ERR2021419.sorted.bam ../align/results/ERR2021420.sorted.bam ../align/results/ERR2021421.sorted.bam ../align/results/ERR2021422.sorted.bam ../align/results/ERR2021423.sorted.bam ../align/results/ERR2021424.sorted.bam ../align/results/ERR2021425.sorted.bam ../align/results/ERR2021426.sorted.bam ../align/results/ERR2021427.sorted.bam ../align/results/ERR2021428.sorted.bam ../align/results/ERR2021429.sorted.bam ../align/results/ERR2021430.sorted.bam ../align/results/ERR2021431.sorted.bam ../align/results/ERR2021432.sorted.bam ../align/results/ERR2021433.sorted.bam ../align/results/ERR2021434.sorted.bam ../align/results/ERR2021435.sorted.bam ../align/results/ERR2021436.sorted.bam ../align/results/ERR2021437.sorted.bam ../align/results/ERR2021438.sorted.bam ../align/results/ERR2021439.sorted.bam ../align/results/ERR2021440.sorted.bam ../align/results/ERR2021441.sorted.bam ../align/results/ERR2021442.sorted.bam ../align/results/ERR2021443.sorted.bam ../align/results/ERR2021444.sorted.bam ../align/results/ERR2021445.sorted.bam ../align/results/ERR2021446.sorted.bam ../align/results/ERR2021447.sorted.bam ../align/results/ERR2021448.sorted.bam ../align/results/ERR2021449.sorted.bam ../align/results/ERR2021450.sorted.bam ../align/results/ERR2021451.sorted.bam ../align/results/ERR2021452.sorted.bam ../align/results/ERR2021453.sorted.bam ../align/results/ERR2021454.sorted.bam ../align/results/ERR2021455.sorted.bam ../align/results/ERR2021456.sorted.bam ../align/results/ERR2021457.sorted.bam ../align/results/ERR2021458.sorted.bam ../align/results/ERR2021459.sorted.bam ../align/results/ERR2021460.sorted.bam ../align/results/ERR2021461.sorted.bam ../align/results/ERR2021462.sorted.bam ../align/results/ERR2021463.sorted.bam ../align/results/ERR2021464.sorted.bam ../align/results/ERR2021465.sorted.bam ../align/results/ERR2021466.sorted.bam ../align/results/ERR2021467.sorted.bam ../align/results/ERR2021468.sorted.bam ../align/results/ERR2021469.sorted.bam ../align/results/ERR2021470.sorted.bam ../align/results/ERR2021471.sorted.bam ../align/results/ERR2021472.sorted.bam ../align/results/ERR2021473.sorted.bam ../align/results/ERR2021474.sorted.bam ../align/results/ERR2021475.sorted.bam ../align/results/ERR2021476.sorted.bam ../align/results/ERR2021477.sorted.bam ../align/results/ERR2021478.sorted.bam ../align/results/ERR2021479.sorted.bam ../align/results/ERR2021480.sorted.bam ../align/results/ERR2021481.sorted.bam ../align/results/ERR2021482.sorted.bam ../align/results/ERR2021483.sorted.bam ../align/results/ERR2021484.sorted.bam ../align/results/ERR2021485.sorted.bam ../align/results/ERR2021486.sorted.bam ../align/results/ERR2021487.sorted.bam ../align/results/ERR2021488.sorted.bam ../align/results/ERR2021489.sorted.bam ../align/results/ERR2021490.sorted.bam ../align/results/ERR2021491.sorted.bam ../align/results/ERR2021492.sorted.bam ../align/results/ERR2021493.sorted.bam ../align/results/ERR2021494.sorted.bam ../align/results/ERR2021495.sorted.bam ../align/results/ERR2021496.sorted.bam ../align/results/ERR2021497.sorted.bam ../align/results/ERR2021498.sorted.bam ../align/results/ERR2021499.sorted.bam ../align/results/ERR2021500.sorted.bam ../align/results/ERR2021501.sorted.bam ../align/results/ERR2021502.sorted.bam