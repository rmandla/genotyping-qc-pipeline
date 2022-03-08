#!/bin/bash

source /broad/software/scripts/useuse
use .shapeit4-4.2.1
chr=$SGE_TASK_ID
mapfile=/humgen/diabetes/users/josep/refPanels/1000G_phase3/1000GP_Phase3/genetic_map_chr${chr}_combined_b37.txt
if [ $chr == 23 ]
then
mapfile=/humgen/diabetes/users/josep/refPanels/1000G_phase3/genetic_map_chrX_nonPAR_combined_b37.txt
fi
directory=/humgen/florezlab/users/rmandla/genotyping-qc/vcf/
inhead=PHBB_SNPs_35K.clean.snps1.mind.0.02.sexcheck_het_out_no_outliers.clean.snps2_no_GCATs-updated-chr${chr}
infile=$directory${inhead}.vcf.gz
outfile=$directory${inhead}-phased.vcf.gz
shapeit4.2 --B infile \
--map $mapfile \
--region $chrom \
--output $outfile \
--thread 8
fi

