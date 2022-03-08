import os, subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from general_qc import *

plink2 = '/humgen/diabetes/users/josep/apps/plink2/plink' # plink2 executable
directory = '/humgen/diabetes2/users/josep/PartnersBIOBANK/merge_all_35K_datasets/merge_35K_genotypes/' # directory containing starting plink bed/bim/fam files
header = 'PHBB_SNPs_35K' # header of starting plink bed/bim/fam files
ancestry_file_path = directory + 'biobank_oct2019_with_subject_ID.csv' # tab-delimited file containing
ancestry_matrix = {'EU':['W'],'HIS':['H'],'AFR':['B'],'AS':['A'],'Other':['O','']} # ancestry matrix, to dictate which entries in ancestry_file_path correspond to which ancestries
pheno_files = ['file1.bin.txt','file2.bin.txt','file3.bin.txt','file4.bin.txt','file5.bin.txt','file6.bin.txt','file7.bin.txt','file8.bin.txt'] # phenotype files for test-missing filtering

run_qc(plink2,directory,header,ancestry_file_path,ancestry_matrix,pheno_files)
