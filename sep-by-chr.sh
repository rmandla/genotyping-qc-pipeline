#!/bin/bash

source /broad/software/scripts/useuse

header=$1
chrom=${SGE_TASK_ID}

/humgen/diabetes/users/josep/apps/plink2/plink --noweb --bfile ${header}.clean.snps1.mind.0.02.sexcheck_het_out.clean.snps2_no_GCATs_no_DUPS.ALL_maf_0.0005 --chr ${chrom} \
--recode vcf --out ${header}.clean.snps1.mind.0.02.sexcheck_het_out.clean.snps2_no_GCATs_no_DUPS.maf_0.0005_chr.${chrom}
