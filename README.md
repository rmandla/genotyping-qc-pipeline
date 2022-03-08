# Pipeline for QC Genotyping Sequencing Data

Ravi M, Philip S, Josep M

This repository contains code and general pipeline of QC-ing, phasing, and imputing sequencing data

## Dependencies

This code was written and tested in Python 3.9 and has the following dependencies:

Python
* Pandas
* Numpy
* seaborn/matplotlib

Other
* plink2
* bcftools
* QCTool
* flashpca
* shapeit4.2

Imputation is performed using TOPMed

## Pipeline

The QC steps follow the following general pipeline:

1. Genotyping vcf files are converted to PLINK bed/bim/fam files
2. Variants from PLINK files are filtered with the PLINK arguments missing, hardy, test-missing, freq, and het
    - Variants within PAR also removed
    - hardy and het are run on ancestry-delimited files. So if a variant fails in one ancestry, then it is removed from every ancestry
3. Sex Check removes additional individuals
4. Run and filter from MDS
5. Remove GCATs
6. Phase
7. Impute

## Scripts

* `general_qc.py` - Python functions useful for removing SNPs and individuals from starting, non-qc'd data
* `mgbb-mega-qc.py` - Python script example of running QC steps 1-5 with the `run_qc` function
* `sep-by-chr.sh` - Shell script for separating plink bed/bim/fam file outputs from `run_qc` by chromosome, and convert to vcf files
* `run-shapeit.sh` - Shell script for running shapeit on `sep-by-chr.sh` outputs
