# Genotyping QC script

import os, sys, subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from shutil import copyfile
plink2 = # path to plink2 executable
directory = # directory containing genotype files
data = # header of genotype files
data_dir = directory + data
ancestries=['EU','HIS','AFR','AS','Other'] # recorded ancestries of a dataset (used for het and hwe filtering)

def run_IBD_QC(file_header):
    imiss = pd.read_table(file_header + '.imiss',delim_whitespace=True, usecols=['FID','IID','F_MISS'])
    imiss['full_id'] = imiss['FID'] + ' ' + imiss['IID']
    imiss = imiss.set_index('full_id')
    genome = pd.read_table(file_header + '.genome',delim_whitespace=True, usecols=['FID1','IID1','FID2','IID2','PI_HAT'])
    genome = genome[genome['PI_HAT']>0.185]
    genome = genome.drop(columns=['PI_HAT'])

    no_dups = genome.drop_duplicates(subset=['FID1','IID1'])
    to_remove = []
    for fid,iid in enumerate(list(no_dups['FID1']),list(no_dups['IID1'])):
        t_comp = genome[(genome['FID1']==fid) & (genome['IID1']==iid)]
        t_comp_miss = imiss[(imiss['FID']==fid) & (imiss['FID']==fid)]['F_MISS'].to_list()[0]
        t_comp['first_id'] = t_com['FID1'] + ' ' + t_com['IID1']
        t_comp['second_id'] = t_com['FID2'] + ' ' + t_com['IID2']
        t_comp = t_comp.set_index('second_id')
        t_comp['2fmiss'] = imiss['F_MISS']
        f_over_s = t_comp[t_comp['2fmiss']<=t_comp_miss]
        s_over_f = t_comp[t_comp['2fmiss']>t_comp_miss]
        to_remove = list(set(to_remove + list(f_over_s['first_id']) + list(s_over_f.index)))
    with open('fail-IBD-QC.txt','w') as out:
        for i in to_remove:
            out.write(i + '\n')

def create_plink(input,output):
    # take in an input vcf file, convert to an output plink file
    subprocess.run([plink2,'--vcf',input,'--out',output],check=True)

def run_multiple_plinks(inpu_files,outpu_files):
    # run create_plink on many files
    for i,j in enumerate(inpu_files,outpu_files):
        create_plink(i,j)

def reformat_bim(bim,output='same',format='snpid'):
    # take in a bim file, reformat it
    # options for reformat are snpid (default) in chr:pos format or concat in chr_pos_ref_alt
    d = pd.read_table(bim,header=None)
    if format=='snpid':
        d[1] = d[0].astype(str) + ':' + d[3].astype(str)
    elif format=='concat':
        d[1] = d[0].astype(str) + '_' + d[3].astype(str) + '_' + d[4] + '_' + d[5]
    else:
        print('invalid format')
        return
    if output == 'same':
        output = bim
    d.to_csv(output,header=None,index=None,sep='\t')

def create_fams(files):
    # take in a list of files, create individual fam files for each of these files, where the phenotype is whether an individual is in
    # said file
    new_df = ''
    phenos = []

    for i in files:
        i = directory + i
        if new_df == '':
            new_df = pd.read_table(i,header=None,delim_whitespace=True)
            phenos.append(list(t[1]))
            new_df = new_df[[0,1,2,3,4]]
        else:
            t = pd.read_table(i,header=None,delim_whitespace=True)
            new_df = new_df[[0,1,2,3,4]]
            phenos.append(list(t[1]))
            new_df = pd.concat([new_df,t])
    new_df[5] = new_df[1]
    new_df = new_df.set_index(5)
    for i in phenos:
        temp_phen_df = pd.DataFrame({0:i,1:[1]*len(i)}).set_index(0)
        new_df[5] = temp_phen_df[1]
        new_df[5] = new_df[5].replace(np.nan,0).astype(int)
        new_df.to_csv(files[phenos.index(i)],header=None,index=None,sep='\t')

def first_anc_pass(input_bfile=None,output_file_ext=''):
    race_cat_defs = {'EU':['W'],'HIS':['H'],'AFR':['B'],'AS':['A'],'Other':['O','']}
    pheno=pd.read_table(directory + "biobank_oct2019_with_subject_ID.csv")
    FAM_phenos=pd.read_table(data_dir + ".fam",header=None)
    FAM_phenos['final_ID']=FAM_phenos[1].str.split('-',expand=True)[1].astype(int)

    for ANC in ancestries:
    # for a specified ANC, create plink files of just people with that ANC (self-reported)
    # options are EU, HIS, AFR, AS, and Other
        if ANC + '_to_include_PHBB_SNPs_35K.txt' not in os.listdir():
            phenos_to_include_ANC = pheno[pheno['race_cat'].isin(race_cat_defs[ANC])] # update for others
            FAM_phenos_anc=FAM_phenos.merge(phenos_to_include_ANC,left_on="final_ID",right_on="Subject ID")
            phenos_to_include_ANC_fam=FAM_phenos_anc.iloc[:, [0,1]].dropna()
            phenos_to_include_ANC_fam.to_csv(ANC + '_to_include_PHBB_SNPs_35K.txt',sep='\t',header=None,index=None)
        if input_bfile==None:
            input_bfile = data_dir
        outb = data +'.' + ANC + output_file_ext
        subprocess.run([plink2,'--bfile',input_bfile,'--keep',ANC + '_to_include_'+data+'.txt','--make-bed','--out',outb],check=True)

def filter_plinkfiles(input_file,pheno_files,anc_input_file_ext,missing_name=True):
    # run a series of filtering steps with plink
    # pheno_files and ancestry_files must be lists
    if '/' in input_file:
        input_header = input_file.split('/')[-1]
    else:
        input_header = input_file
    for i in ['missing','hardy','test-missing','freq','het']: # missing, test-missing, freq run all together; no hardy/het on others; test-missing to compare array batch
        if i != 'hardy' and i != 'het': # run on missing, test-missing, and freq
            if i == 'test-missing':
                # need phenotype data, so use array batches for initial run
                for j in pheno_files:
                    subprocess.run([plink2,'--bfile',input_file,'--pheno',j,'--'+i,'--out',input_header+'-'+j.split('/')[-1]],check=True)
            if i == 'missing' and missing_name==True: # after hardy, combine fails for all ancestry for exclusion
                subprocess.run([plink2,'--bfile',input_file,'--'+i,'--out',input_header],check=True)
            else:
                subprocess.run([plink2,'--bfile',input_file,'--' + i,'--out','.'],check=True)
        elif i == 'hardy' or i == 'het':
            for ANC in ancestries:
                anc_input = data+'.'+ANC+anc_input_file_ext
                anc_output = input_header + '.' + ANC
                subprocess.run([plink2,'--bfile',anc_input,'--'+i,'--out',anc_output],check=True)

def filter_tables(file='',output_name_partial='to_remove_',snp_outname=''):
    # take in files outputted from filter_plinkfiles, filter based on specific criteria and output a list of SNPs for exclusion
    # possible anc include EU, HIS, AFR, AS
    if file == '':
        header = data
    else:
        header = file
    snps_for_exclusion = []
    for i in os.listdir():
        footer = '.' + i.split('.')[-1]
        output_name = output_name_partial + footer + '.txt'
        if file in i:
            if footer == '.lmiss':
                d = pd.read_table(i,delim_whitespace=True)
                d = d[d['F_MISS']>0.05]
                output='missing'
                snps_for_exclusion += d['SNP'].to_list()
            elif footer == '.missing' and 'txt' in i:
                d = pd.read_table(i,delim_whitespace=True)
                d = d[d['P']<0.00005]
                output='test_missing'
                snps_for_exclusion += d['SNP'].to_list()
            elif footer == '.frq':
                d = pd.read_table(i,delim_whitespace=True)
                d[d['MAF']<0.05].to_csv(header +'.to_removeForPCA_freq.txt',sep='\t') # SNPs to be excluded for PCA
                d = d[d['MAF']<0.0005]
                output='freq'
                snps_for_exclusion += d['SNP'].to_list()
            elif footer == '.hwe':
                d = pd.read_table(i,delim_whitespace=True)
                d = d[(d['P']<1e-10) & (d['TEST']=='UNAFF')]
                output='hardy'
                snps_for_exclusion += d['SNP'].to_list()
            try:
                d.to_csv(output_name,sep='\t',index=None)
            except:
                pass

    # filter out PAR # add diff geno coordiantes
    bim = pd.read_table(data_dir + '.bim',header=None)
    bim_x = bim[bim[0].astype(str)=='23']
    to_remove_par1 = bim_x[(bim_x[3]<2699520) & (bim_x[3]>60001)]
    to_remove_par2 = bim_x[(bim_x[3]<155260560) & (bim_x[3]>154931044)]
    for i in [to_remove_par1,to_remove_par2]:
        snps_for_exclusion += list(i[1])

    # export unique list of SNPs
    snps_for_exclusion = list(set(snps_for_exclusion))
    with open(file + snp_outname + 'all_SNPs_to_remove_u.txt','w') as out:
        for i in snps_for_exclusion:
            out.write(i + '\n')

# after, 865644 are to be removed
def filter_individual_snps():
    # filter out SNPs that need to be removed from previous filtering steps. Remove SNPs for PCA also
    # anc values can be EU, AFR, AS, or HIS
    subprocess.run([plink2,'--bfile',data_dir,'--exclude','all_SNPs_to_remove_u.txt','--make-bed','--out',data+'.clean.snps1'],check=True)
    subprocess.run([plink2,'--bfile',data+'.clean.snps1','--exclude',data+'.to_removeForPCA_freq.txt','--make-bed','--out',data+'.clean.snps1_MAF_05_forMDS'],check=True)

def sex_check():
    # filter out SNPs that will be used for packages
    # anc values can be EU, AFR, AS, or HIS

    subprocess.run([plink2,'--bfile',data+'.clean.snps1','--check-sex','--out',data+'.clean.snps1.sex_check'],check=True)
    subprocess.run([plink2,'--bfile',data+'.clean.snps1','--mind','0.02','--make-bed','--out',data+'.clean.snps1.mind.0.02'],check=True)
    first_anc_pass(data+'.clean.snps1','.clean.snps1')

    sex_check_samples = pd.read_table(data + '.clean.snps1.sex_check.sexcheck',delim_whitespace=True)
    to_remove_sex_check_samples = sex_check_samples[(sex_check_samples['STATUS']=='PROBLEM') & (sex_check_samples['PEDSEX'] != 0)] # when pedsex != genetic sex
    to_remove_sex_check_samples = to_remove_sex_check_samples.iloc[:, [0,1]].drop_duplicates()
    to_remove_sex_check_samples.to_csv(data+'_samples_to_sex_check.txt',sep='\t',index=None)

    subprocess.run([plink2,'--bfile',data+'.clean.snps1.mind.0.02','--remove',data+'_samples_to_sex_check.txt','--make-bed','--out',data+'.clean.snps1.mind.0.02.sexcheck'],check=True)

    for anc in ['EU','HIS','AFR','AS']:
        header = data + '.' + anc
        subprocess.run([plink2,'--bfile',header+'.clean.snps1','--allow-no-sex','--het','--out',header+'.clean.snps1.mind.0.02'],check=True) # by ancestry
        het_samples = pd.read_table(header+'.clean.snps1.mind.0.02.het',delim_whitespace=True)
        mean = np.mean(het_samples['F'])
        std = np.std(het_samples['F'])
        max_val = mean + std * 4
        min_val = mean - std * 4
        het_samples_outliers = het_samples[(het_samples['F']>max_val) | (het_samples['F']<min_val)].iloc[:, [0,1]]
        het_samples_outliers.to_csv(header+'_samples_to_outliers_het.txt',sep='\t',index=None)
        fig, axes = plt.subplots(1, 2)
        sns.histplot(het_samples['F'],ax=axes[0])
        sns.boxplot(het_samples['F'],ax=axes[1])
        plt.savefig(header+'hist_F_het.pdf')
        plt.clf()

    for anc in ['EU','HIS','AFR','AS']:
        # concat all subjects that are outliers from het for removal
        header = data + '.' + anc
        if anc == 'EU':
            subprocess.run('cat ' + header + '_samples_to_outliers_het.txt > ' + data + '_samples_to_outliers_het.txt',shell=True,check=True)
        else:
            subprocess.run('tail -n +2 ' + header + '_samples_to_outliers_het.txt >> ' + data + '_samples_to_outliers_het.txt',shell=True,check=True)
    subprocess.run([plink2,'--bfile',data+'.clean.snps1.mind.0.02.sexcheck','--remove',data+'_samples_to_outliers_het.txt','--make-bed','--out',data+'.clean.snps1.mind.0.02.sexcheck_het_out'],check=True)
    subprocess.run([plink2,'--bfile',data+'.clean.snps1.mind.0.02.sexcheck_het_out','--remove',data+'_samples_to_sex_check.txt','--make-bed','--out',data+'.clean.snps1.mind.0.02.sexcheck'],check=True)
    subprocess.run([plink2,'--bfile',data+'.clean.snps1.mind.0.02.sexcheck_het_out','--exclude',data+'.to_removeForPCA_freq.txt','--make-bed','--out',data+'.clean.snps1.mind.0.02.sexcheck_het_out_MAF_05_forMDS'],check=True)
    subprocess.run([plink2,'--bfile',data+'.clean.snps1.mind.0.02.sexcheck_het_out_MAF_05_forMDS','--exclude','range','/humgen/diabetes/users/josep/PartnersBIOBANK/high-LD-regions_hg19.txt','--indep-pairwise','50','5','0.2','--out',data+'.clean.snps1.mind.0.02.sexcheck_het_out_MAF_05_forMDS_prunnedSNPs'],check=True)
    subprocess.run([plink2,'--bfile',data+'.clean.snps1.mind.0.02.sexcheck_het_out_MAF_05_forMDS','--extract',data+'.clean.snps1.mind.0.02.sexcheck_het_out_MAF_05_forMDS_prunnedSNPs.prune.in','--genome','--out',data+'.clean.snps1.mind.0.02.sexcheck_het_out_MAF_05_forMDS_prunnedSNPs'],check=True)

def check_relatedness():
    subprocess.run("awk '{ if ($10 > 0.185) { print } }' "+data+".clean.snps1.mind.0.02.sexcheck_het_out_MAF_05_forMDS_prunnedSNPs.genome > "+data+"pairs_0.185",shell=True,check=True)
    pairs_UE = pd.read_table(data+"pairs_0.185",delim_whitespace=True)
    intervals = np.linspace(1,0.185,num=10)
    results = []
    for i in range(len(intervals)-1):
        up = intervals[i]
        down = intervals[i+1]
        PI_HAT_interval=(str(round(up,3)) + '-' + str(round(down,3)))
        count = len(pairs_UE[(pairs_UE['PI_HAT']<up) & (pairs_UE['PI_HAT']>=down)])
        results.append([PI_HAT_interval,count])
    pd.DataFrame(results).to_csv(data + '_pi_hat_0.18.txt',index=None,header=None,sep='\t')

    copyfile(data+'.imiss',data+'.clean.snps1.mind.0.02.sexcheck_het_out_MAF_05_forMDS_prunnedSNPs.imiss') #rename for perl script
    subprocess.run(['perl','/humgen/diabetes/users/josep/PartnersBIOBANK/run-IBD-QC.pl',data+'.clean.snps1.mind.0.02.sexcheck_het_out_MAF_05_forMDS_prunnedSNPs'],check=True)
    for i in ['','_MAF_05_forMDS']:
        file = data + '.clean.snps1.mind.0.02.sexcheck_het_out' + i
        outfile = file+'_unrelated'
        subprocess.run([plink2,'--bfile',file,'--remove','fail-IBD-QC.txt','--make-bed','--out',outfile],check=True)

def extract_nonGCAT(return_gcats=False):
    if return_gcats:
        bfile = data + '.clean.snps1.mind.0.02.sexcheck_het_out_no_outliers.clean.snps2.bim'
    else:
        bfile = "/humgen/diabetes/users/josep/PartnersBIOBANK/phase3.pruned.bim" # for 1000G pca
    bim = pd.read_table(bfile,header=None,delim_whitespace=True)
    bim['geno'] = bim[4]+'_'+bim[5]
    if return_gcats:
        bim_gcat=bim[~bim['geno'].isin(['A_C','A_G','C_A','C_T','G_A','G_T'])]
        bim_gcat.drop(columns=['geno']).to_csv(data+'_GCATs_to_remove_RELATED.txt',sep='\t',header=None,index=None)
        return
    bim_no_gcat=bim[bim['geno'].isin(['A_C','A_G','C_A','C_T','G_A','G_T'])]
    bim_no_gcat['SNP'] = bim_no_gcat[0].astype(str) + '_' + bim_no_gcat[3].astype(str) + '_' + bim_no_gcat['geno']
    bim_partners = pd.read_table(data+'.clean.snps1.mind.0.02.sexcheck_het_out_MAF_05_forMDS_unrelated.bim',header=None)
    bim_partners['SNP'] = bim_partners[0].astype(str) + '_' + bim_partners[3].astype(str) + '_' + bim_partners[4]+'_'+bim_partners[5]
    bim_not_gcat_bim_partners = bim_no_gcat.merge(bim_partners,left_on='SNP',right_on='SNP')
    bim_not_gcat_bim_partners[['SNP']].to_csv(data+'.phase3.not_GCAT_SNPs.txt',index=None)

def prep_for_mds():
    # needs a lot of memory (48G)

    subprocess.run([plink2,'--bim',data+'_new_bim_for_merge1000G.bim','--bed',data+'.clean.snps1.mind.0.02.sexcheck_het_out.bed','--fam',data+'.clean.snps1.mind.0.02.sexcheck_het_out.fam','--extract',data+'.phase3.not_GCAT_SNPs.txt','--make-bed','--out',data+'.clean.snps1_MAF_05_forMDS_RELATED_phase3_not_gcat_SNPs'],check=True)
    subprocess.run([plink2,'--allow-no-sex','--bfile',data+'.clean.snps1_MAF_05_forMDS_RELATED_phase3_not_gcat_SNPs','--bmerge','/humgen/diabetes/users/josep/PartnersBIOBANK/phase3.pruned.bed','1000G_chr_pos_alleles.bim','/humgen/diabetes/users/josep/PartnersBIOBANK/phase3.pruned.fam','--extract',data+'.phase3.not_GCAT_SNPs.txt','--make-bed','--out',data+'.clean.snps1_MAF_05_forMDS_RELATED_phase3_not_gcat_SNPs_with_phase3'],check=True)

    subprocess.run(['bash','run_flashpca.sh',data+'.clean.snps1_MAF_05_forMDS_RELATED_phase3_not_gcat_SNPs_with_phase3',data+'_PCs_ALL_RELATED_with_1000G','yes'],check=True)
    # run MDS script

def filter_from_mds():
    subprocess.run([plink2,'--bfile',data+'.clean.snps1.mind.0.02.sexcheck_het_out','--remove',data+'.clean.snps1.mind.0.02.irem','--make-bed','--out',data+'.clean.snps1.mind.0.02.sexcheck_het_out_no_outliers'],check=True)

    subprocess.run([plink2,'--bfile',data+'.clean.snps1_MAF_05_forMDS_RELATED_phase3_not_gcat_SNPs_with_phase3_no_outliers','--logistic','--covar','PHBB_SNPs_35K_PCs_ALL_RELATED_with_1000G.txt','--covar-name','PC1 PC2 PC3 PC4 PC5 PC6 PC7','--out',data+'.logiscic_association_test'])

    first_anc_pass(data+'.clean.snps1.mind.0.02.sexcheck_het_out_no_outliers',output_file_ext='pass2')
    filter_plinkfiles(data+'.clean.snps1.mind.0.02.sexcheck_het_out_no_outliers',['run1/file1.bin.txt','run1/file2.bin.txt','run1/file3.bin.txt','run1/file4.bin.txt','run1/file5.bin.txt','run1/file6.bin.txt','run1/file7.bin.txt','run1/file8.bin.txt'],anc_input_file_ext='pass2')
    filter_tables(file='PHBB_SNPs_35K.clean.snps1.mind.0.02.sexcheck_het_out_no_outliers')

    subprocess.run([plink2,'--bfile',data+'.clean.snps1.mind.0.02.sexcheck_het_out_no_outliers','--exclude',data+'.clean.snps1.mind.0.02.sexcheck_het_out_no_outliers'+'all_SNPs_to_remove_u.txt','--make-bed','--out',data+'.clean.snps1.mind.0.02.sexcheck_het_out_no_outliers.clean.snps2'],check=True)
    subprocess.run([plink2,'--bfile',data+'.clean.snps1.mind.0.02.sexcheck_het_out_no_outliers.clean.snps2','--logistic','--covar','PHBB_SNPs_35K_PCs_ALL_RELATED_with_1000G.txt','--covar-name','PC1 PC2 PC3 PC4 PC5 PC6 PC7','--out',data+'.clean.snps1.mind.0.02.sexcheck_het_out_no_outliers.clean.snps2'])

def identify_gcats():
    extract_nonGCAT(return_gcats=True)
    subprocess.run([plink2,'--bfile',data+'.clean.snps1.mind.0.02.sexcheck_het_out_no_outliers.clean.snps2','--exclude',data+'_GCATs_to_remove_RELATED.txt','--make-bed','--out',data+'.clean.snps1.mind.0.02.sexcheck_het_out_no_outliers.clean.snps2_no_GCATs'],check=True)
    subprocess.run([plink2,'--noweb','--bfile',data+'.clean.snps1.mind.0.02.sexcheck_het_out_no_outliers.clean.snps2_no_GCATs','--chr','23','--make-bed','--out',data+'.clean.snps1.mind.0.02.sexcheck_het_out_unrelated_no_outliers_no_GCATs_no_GCATs_chr_X'],check=True)

def prep_frequency_file(bfile):
    subprocess.run([plink2,'--bfile',data+bfile,'--freq','--out',data+bfile+'.freq'],check=True)
    subprocess.run(['perl','x','-b',data+bfile+'.bim','-f',data+bfile+'.freq.frq','-r','/humgen/diabetes/users/josep/PartnersBIOBANK/HRC_info/HRC.r1-1.GRCh37.wgs.mac5.sites.tab','-h'],check=True)
    subprocess.run(['cp','/humgen/diabetes/users/josep/PartnersBIOBANK/PartnersMerged/Run-plink.sh','Run-plink.sh'],check=True)
    subprocess.run("sed -i 's/x4927_x5352_x4784.clean.snps1.mind.0.02.sexcheck_unrelated_no_outliers_no_GCATs-updated/" + data + ".clean.snps1.mind.0.02.sexcheck_het_out_no_outliers.clean.snps2_no_GCATs/g' Run-plink.sh",shell=True,check=True)
    subprocess.run(['sed','-i',"'s/x4927_x5352_x4784.clean.snps1.mind.0.02.sexcheck_unrelated_no_outliers_no_GCATs/"+data+".clean.snps1.mind.0.02.sexcheck_het_out_no_outliers.clean.snps2_no_GCATs/g'",'Run-plink.sh'],check=True)
    subprocess.run(['sh','Run-plink.sh'],check=True)

    subprocess.run(['bash','run_flashpca.sh',data + '.clean.snps1.mind.0.02.sexcheck_het_out_no_outliers.clean.snps2_no_GCATs',data+'_PCs_RELATED','20'],check=True)
