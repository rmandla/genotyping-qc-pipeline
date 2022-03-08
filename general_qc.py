import os, subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

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

def first_anc_pass(plink2,directory,header,ancestry_file_path,ancestry_matrix,ancestries,ancestry_file_anccolname='race_cat',ancestry_file_idcolname='Subject ID',input_bfile=None):
    pheno=pd.read_table(ancestry_file_path)
    data_dir = directory+header
    FAM_phenos=pd.read_table(data_dir + ".fam",header=None,delim_whitespace=True)
    FAM_phenos['final_ID']=FAM_phenos[1].str.split('-',expand=True)[1].astype(int)
    ancestries = list(ancestry_matrix.keys())

    for ANC in ancestries:
    # for a specified ANC, create plink files of just people with that ANC (self-reported)
        if ANC + '_to_include_' + header + '.txt' not in os.listdir():
            phenos_to_include_ANC = pheno[pheno[ancestry_file_anccolname].isin(ancestry_matrix[ANC])] # update for others
            FAM_phenos_anc=FAM_phenos.merge(phenos_to_include_ANC,left_on="final_ID",right_on=ancestry_file_idcolname)
            phenos_to_include_ANC_fam=FAM_phenos_anc.iloc[:, [0,1]].dropna()
            phenos_to_include_ANC_fam.to_csv(ANC + '_to_include_' + header + '.txt',sep='\t',header=None,index=None)
        if input_bfile==None:
            input_bfile = data_dir
        out = header +'.' + ANC
        subprocess.run([plink2,'--bfile',input_bfile,'--keep',ANC + '_to_include_'+header+'.txt','--make-bed','--out',out],check=True)

def filter_plinkfiles(plink2,directory,header,pheno_files,ancestries,missing=True,test_missing=True,freq=True,hwe=True):
    # run a series of filtering steps with plink
    # pheno_files and ancestry_files must be lists
    data_dir = directory+header

    tests = []
    if missing:
        tests.append('missing')
    if test_missing:
        tests.append('test-missing')
    if freq:
        tests.append('freq')
    if hwe:
        tests.append('hardy')

    for i in tests: # missing, test-missing, freq run all together; no hardy/het on others; test-missing to compare array batch
        if i != 'hardy': # run on missing, test-missing, and freq
            if i == 'test-missing':
                # need phenotype data, so use array batches for initial run
                for j in pheno_files:
                    if '/' in j:
                        out_tail = j.split('/')[-1]
                    else:
                        out_tail = j
                    subprocess.run([plink2,'--bfile',data_dir,'--pheno',j,'--'+i,'--out',header+'-'+out_tail],check=True)
            else:
                subprocess.run([plink2,'--bfile',data_dir,'--' + i,'--out',header],check=True)
        else:
            for ANC in ancestries:
                anc_header = header+'.'+ANC
                subprocess.run([plink2,'--bfile',anc_header,'--'+i,'--out',anc_header],check=True)

def filter_tables(plink2,directory,header,ancestries,miss_cutoff=0.05,test_missing_cutoff=0.00005,maf_cutoff=0.0005,hwe_cutoff=1e-10):
    # take in files outputted from filter_plinkfiles, filter based on specific criteria and output a list of SNPs for exclusion
    snps_for_exclusion = []
    footers = ['.lmiss','.missing','.frq','.hwe']
    for file in os.listdir():
        footer = '.'+file.split('.')[-1]
        if footer in footers:
            out_name = file + '.to_remove.txt'
            if footer == '.lmiss':
                d = pd.read_table(file,delim_whitespace=True)
                d = d[d['F_MISS']>miss_cutoff]
                snps_for_exclusion += d['SNP'].to_list()
            elif footer == '.missing':
                d = pd.read_table(file,delim_whitespace=True)
                d = d[d['P']<test_missing_cutoff]
                snps_for_exclusion += d['SNP'].to_list()
            elif footer == '.frq':
                d = pd.read_table(file,delim_whitespace=True)
                d[d['MAF']<0.05].to_csv(header +'.to_removeForPCA_freq.txt',sep='\t') # SNPs to be excluded for PCA
                d = d[d['MAF']<maf_cutoff]
                snps_for_exclusion += d['SNP'].to_list()
            elif footer == '.hwe':
                for ANC in ancestries:
                    d = pd.read_table(file,delim_whitespace=True)
                    d = d[(d['P']<hwe_cutoff) & (d['TEST']=='UNAFF')]
                    snps_for_exclusion += d['SNP'].to_list()
            try:
                d.to_csv(out_name,sep='\t',index=None)
            except:
                pass

    # filter out PAR # add diff geno coordiantes
    bim = pd.read_table(directory+header + '.bim',header=None)
    bim_x = bim[bim[0].astype(str)=='23']
    to_remove_par1 = bim_x[(bim_x[3]<2699520) & (bim_x[3]>60001)]
    to_remove_par2 = bim_x[(bim_x[3]<155260560) & (bim_x[3]>154931044)]
    for i in [to_remove_par1,to_remove_par2]:
        snps_for_exclusion += list(i[1])

    # export unique list of SNPs
    snps_for_exclusion = list(set(snps_for_exclusion))
    with open(header + '.all_SNPs_to_remove_u.txt','w') as out:
        for i in snps_for_exclusion:
            out.write(i + '\n')

    data_dir = directory + header
    subprocess.run([plink2,'--bfile',data_dir,'--exclude',header + '.all_SNPs_to_remove_u.txt','--make-bed','--out',header+'.clean.snps1'],check=True)
    subprocess.run([plink2,'--bfile',header+'.clean.snps1','--exclude',header+'.to_removeForPCA_freq.txt','--make-bed','--out',header+'.clean.snps1_MAF_05_forMDS'],check=True)

def sex_check(plink2,header,ancestries,ancestry_file_path,ancestry_matrix,ancestry_file_anccolname='race_cat',ancestry_file_idcolname='Subject ID',input_bfile=None,mds=False):
    subprocess.run([plink2,'--bfile',header+'.clean.snps1','--check-sex','--out',header+'.clean.snps1'],check=True)
    subprocess.run([plink2,'--bfile',header+'.clean.snps1','--mind','0.02','--make-bed','--out',header+'.clean.snps1.mind.0.02'],check=True)
    first_anc_pass(plink2,'./',header+'.clean.snps1',ancestry_file_path,ancestry_matrix,ancestries,ancestry_file_anccolname=ancestry_file_anccolname,ancestry_file_idcolname=ancestry_file_idcolname)

    sex_check_samples = pd.read_table(header + '.clean.snps1.sexcheck',delim_whitespace=True)
    to_remove_sex_check_samples = sex_check_samples[(sex_check_samples['STATUS']=='PROBLEM') & (sex_check_samples['PEDSEX'] != 0)] # when pedsex != genetic sex
    to_remove_sex_check_samples = to_remove_sex_check_samples.iloc[:, [0,1]].drop_duplicates()
    to_remove_sex_check_samples.to_csv(header+'_samples_to_sex_check.txt',sep='\t',index=None)

    subprocess.run([plink2,'--bfile',header+'.clean.snps1.mind.0.02','--remove',header+'_samples_to_sex_check.txt','--make-bed','--out',header+'.clean.snps1.mind.0.02.sexcheck'],check=True)

    for ANC in ancestries:
        subprocess.run([plink2,'--bfile',header+'.clean.snps1.'+ANC,'--allow-no-sex','--het','--out',header+'.clean.snps1.mind.0.02.'+ANC],check=True) # by ancestry
        het_samples = pd.read_table(header+'.clean.snps1.mind.0.02.'+ANC+'.het',delim_whitespace=True)
        mean = np.mean(het_samples['F'])
        std = np.std(het_samples['F'])
        max_val = mean + std * 4
        min_val = mean - std * 4
        het_samples_outliers = het_samples[(het_samples['F']>max_val) | (het_samples['F']<min_val)].iloc[:, [0,1]]
        het_samples_outliers.to_csv(ANC+'_samples_to_outliers_het.txt',sep='\t',index=None)
        fig, axes = plt.subplots(1, 2)
        sns.histplot(het_samples['F'],ax=axes[0])
        sns.boxplot(het_samples['F'],ax=axes[1])
        plt.savefig(header+'-hist_F_het.pdf')
        plt.clf()

    for ANC in ancestries:
        # concat all subjects that are outliers from het for removal
        if ANC == 'EU':
            subprocess.run('cat ' + ANC + '_samples_to_outliers_het.txt > ' + header + '_samples_to_outliers_het.txt',shell=True,check=True)
        else:
            subprocess.run('tail -n +2 ' + ANC + '_samples_to_outliers_het.txt >> ' + header + '_samples_to_outliers_het.txt',shell=True,check=True)

    subprocess.run([plink2,'--bfile',header+'.clean.snps1.mind.0.02.sexcheck','--remove',header+'_samples_to_outliers_het.txt','--make-bed','--out',header+'.clean.snps1.mind.0.02.sexcheck_het_out'],check=True)
    if mds:
        subprocess.run([plink2,'--bfile',header+'.clean.snps1.mind.0.02.sexcheck_het_out','--exclude',header+'.to_removeForPCA_freq.txt','--make-bed','--out',header+'.clean.snps1.mind.0.02.sexcheck_het_out_MAF_05_forMDS'],check=True)
        subprocess.run([plink2,'--bfile',header+'.clean.snps1.mind.0.02.sexcheck_het_out_MAF_05_forMDS','--exclude','range','/humgen/diabetes/users/josep/PartnersBIOBANK/high-LD-regions_hg19.txt','--indep-pairwise','50','5','0.2','--out',header+'.clean.snps1.mind.0.02.sexcheck_het_out_MAF_05_forMDS_prunnedSNPs'],check=True)
        subprocess.run([plink2,'--bfile',header+'.clean.snps1.mind.0.02.sexcheck_het_out_MAF_05_forMDS','--extract',header+'.clean.snps1.mind.0.02.sexcheck_het_out_MAF_05_forMDS_prunnedSNPs.prune.in','--genome','--out',header+'.clean.snps1.mind.0.02.sexcheck_het_out_MAF_05_forMDS_prunnedSNPs'],check=True)

def check_relatedness(plink2,header):
    subprocess.run("awk '{ if ($10 > 0.185) { print } }' "+header+".clean.snps1.mind.0.02.sexcheck_het_out_MAF_05_forMDS_prunnedSNPs.genome > "+header+"pairs_0.185",shell=True,check=True)
    pairs_UE = pd.read_table(header+"pairs_0.185",delim_whitespace=True)
    intervals = np.linspace(1,0.185,num=10)
    results = []
    for i in range(len(intervals)-1):
        up = intervals[i]
        down = intervals[i+1]
        PI_HAT_interval=(str(round(up,3)) + '-' + str(round(down,3)))
        count = len(pairs_UE[(pairs_UE['PI_HAT']<up) & (pairs_UE['PI_HAT']>=down)])
        results.append([PI_HAT_interval,count])
    pd.DataFrame(results).to_csv(header + '_pi_hat_0.18.txt',index=None,header=None,sep='\t')

    run_IBD_QC(header+'.clean.snps1.mind.0.02.sexcheck_het_out_MAF_05_forMDS_prunnedSNPs')
    for i in ['','_MAF_05_forMDS']:
        file = header + '.clean.snps1.mind.0.02.sexcheck_het_out' + i
        outfile = file+'_unrelated'
        subprocess.run([plink2,'--bfile',file,'--remove','fail-IBD-QC.txt','--make-bed','--out',outfile],check=True)

def extract_nonGCAT(header,pca_bimfile=None,return_gcats=True):
    if return_gcats:
        bfile = header + '.clean.snps1.mind.0.02.sexcheck_het_out.bim'
    else:
        bfile = pca_bimfile
    bim = pd.read_table(bfile,header=None,delim_whitespace=True)
    bim['geno'] = bim[4]+'_'+bim[5]
    if return_gcats:
        bim_gcat=bim[~bim['geno'].isin(['A_C','A_G','C_A','C_T','G_A','G_T'])]
        bim_gcat[[1]].to_csv(header+'_GCATs_to_remove.txt',sep='\t',header=None,index=None)
        return
    bim_no_gcat=bim[bim['geno'].isin(['A_C','A_G','C_A','C_T','G_A','G_T'])]
    bim_no_gcat['SNP'] = bim_no_gcat[0].astype(str) + '_' + bim_no_gcat[3].astype(str) + '_' + bim_no_gcat['geno']
    bim_geno = pd.read_table(header+'.clean.snps1.mind.0.02.sexcheck_het_out_MAF_05_forMDS_unrelated.bim',header=None)
    bim_geno['SNP'] = bim_geno[0].astype(str) + '_' + bim_geno[3].astype(str) + '_' + bim_geno[4]+'_'+bim_geno[5]
    bim_not_gcat_bim_geno = bim_no_gcat.merge(bim_geno,left_on='SNP',right_on='SNP')
    bim_not_gcat_bim_geno[['SNP']].to_csv(header+'.phase3.not_GCAT_SNPs.txt',index=None)

def prep_for_mds(plink2,header,pca_bedfile):
    subprocess.run([plink2,'--bim',header+'_new_bim_for_merge1000G.bim','--bed',header+'.clean.snps1.mind.0.02.sexcheck_het_out.bed','--fam',header+'.clean.snps1.mind.0.02.sexcheck_het_out.fam','--extract',header+'.phase3.not_GCAT_SNPs.txt','--make-bed','--out',header+'.clean.snps1_MAF_05_forMDS_RELATED_phase3_not_gcat_SNPs'],check=True)
    subprocess.run([plink2,'--allow-no-sex','--bfile',header+'.clean.snps1_MAF_05_forMDS_RELATED_phase3_not_gcat_SNPs','--bmerge','/humgen/diabetes/users/josep/PartnersBIOBANK/phase3.pruned.bed','1000G_chr_pos_alleles.bim','/humgen/diabetes/users/josep/PartnersBIOBANK/phase3.pruned.fam','--extract',data+'.phase3.not_GCAT_SNPs.txt','--make-bed','--out',data+'.clean.snps1_MAF_05_forMDS_RELATED_phase3_not_gcat_SNPs_with_phase3'],check=True)

    subprocess.run(['bash','run_flashpca.sh',header+'.clean.snps1_MAF_05_forMDS_RELATED_phase3_not_gcat_SNPs_with_phase3',header+'_PCs_ALL_RELATED_with_1000G','yes'],check=True)
    # run MDS script

def filter_from_mds(plink2,header):
    filter_plinkfiles(header+'.clean.snps1.mind.0.02.sexcheck_het_out',missing=False,test_missing=False,het=False,hwe=False)
    filter_tables(header+'.clean.snps1.mind.0.02.sexcheck_het_out')

    subprocess.run([plink2,'--bfile',header+'.clean.snps1.mind.0.02.sexcheck_het_out','--exclude',header+'.clean.snps1.mind.0.02.sexcheck_het_out'+'all_SNPs_to_remove_u.txt','--make-bed','--out',header+'.clean.snps1.mind.0.02.sexcheck_het_out.clean.snps2'],check=True)

def remove_gcat_dups(plink2,header):
    extract_nonGCAT(header,return_gcats=True)
    subprocess.run([plink2,'--bfile',header+'.clean.snps1.mind.0.02.sexcheck_het_out','--exclude',header+'_GCATs_to_remove.txt','--make-bed','--out',header+'.clean.snps1.mind.0.02.sexcheck_het_out_no_GCATs'],check=True)
    bim = pd.read_table(header+'.clean.snps1.mind.0.02.sexcheck_het_out_no_GCATs.bim',delim_whitespace=True,header=None)
    dups = bim[bim[3].duplicated()]
    dups[[1]].to_csv(header+'_duplicated_SNPs_to_remove.txt',index=None,header=None)
    subprocess.run([plink2,'--bfile',header+'.clean.snps1.mind.0.02.sexcheck_het_out_no_GCATs','--exclude',header+'_duplicated_SNPs_to_remove.txt','--make-bed','--out',header+'.clean.snps1.mind.0.02.sexcheck_het_out_no_GCATs_no_DUPS'],check=True)

def run_qc(plink2,directory,header,ancestry_file_path,ancestry_matrix,pheno_files,pca_bimfile=None,ancestry_file_anccolname='race_cat',ancestry_file_idcolname='Subject ID',miss_cutoff=0.05,test_missing_cutoff=0.00005,maf_cutoff=0.0005,hwe_cutoff=1e-10,mds=False):
    ancestries = list(ancestry_matrix.keys())
    first_anc_pass(plink2,directory,header,ancestry_file_path,ancestry_matrix,ancestries,ancestry_file_anccolname=ancestry_file_anccolname,ancestry_file_idcolname=ancestry_file_idcolname)
    filter_plinkfiles(plink2,directory,header,pheno_files,ancestries,missing=True,test_missing=True,freq=True,hwe=True)
    filter_tables(plink2,directory,header,ancestries,miss_cutoff=miss_cutoff,test_missing_cutoff=test_missing_cutoff,maf_cutoff=maf_cutoff,hwe_cutoff=hwe_cutoff)
    sex_check(plink2,header,ancestries,ancestry_file_path,ancestry_matrix,ancestry_file_anccolname=ancestry_file_anccolname,ancestry_file_idcolname=ancestry_file_idcolname)
    remove_gcat_dups(plink2,header)
    subprocess.run([plink2,'--bfile',header+'.clean.snps1.mind.0.02.sexcheck_het_out_no_GCATs_no_DUPS','--maf','0.0005','--make-bed','--out',header+'.clean.snps1.mind.0.02.sexcheck_het_out_no_GCATs_no_DUPS_maf_0.0005'])
