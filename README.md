# Pipeline for QC Genotyping Sequencing Data

Ravi M, Philip S, Josep M

This repository contains code and general pipeline of QC-ing, phasing, and imputing sequencing data

## Dependencies

This code was written and tested in Python 3.9 and has the following dependencies:

Python
* Pandas
* Numpy

Other
* plink2
* bcftools
* QCTool
* flashpca
* shapeit4.2

Optional
* seaborn/matplotlib

Imputation is performed using TOPMed

## Pipeline

The QC steps follow the following general pipeline:

1. Genotyping vcf files are converted to PLINK bed/bim/fam files
2. Variants from PLINK files are filtered with the PLINK arguments missing, hardy, test-missing, freq, and het
a. Variants within PAR also removed
b. hardy and het are run on ancestry-delimited files. So if a variant fails in one ancestry, then it is removed from every ancestry
3. Sex Check removes additional individuals
4. Run and filter from MDS
5. Remove GCATs
6. Phase
7. Impute
