#!/bin/bash

source /broad/software/scripts/useuse

use Tabix
header=$1
MAF=0.0005
chrom=${SGE_TASK_ID}

vcf=${header}.clean.snps1.mind.0.02.sexcheck_het_out_no_GCATs_no_DUPS_maf_${MAF}_chr.${chrom}.vcf

bgzip ${vcf}
tabix -p vcf ${vcf}.gz
