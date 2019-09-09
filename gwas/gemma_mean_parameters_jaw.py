# Use this script to calculate mean parameter values across 10 runs of GEMMA's BSLMM for jaw size

import pandas as pd
import numpy as np

trait_name = 'jaw_'
vcf = 'all_pupfish_filtered_snps_passed.Q20.MAF0.05.MAXMISS0.5.plink_sansal_run_'
working_dir = '/pine/scr/j/m/jmcgirr/gemma_emilie/scripts/sansal/jaw/output/'

# start loop with first param file
j=str(1)
param = pd.read_table(working_dir +vcf+j+'.param.txt')
param['beta'+j]= param['beta']
param['gamma'+j]= param['gamma']
param['alpha'+j]= param['alpha']
param = param[['rs','alpha'+j,'beta'+j,'gamma'+j]]

# loop through all param files

iters = np.arange(2,11)

for i in iters: 
    
    j=str(i)
    parami = pd.read_table(working_dir +vcf+j+'.param.txt')
    parami['beta'+j]= parami['beta']
    parami['gamma'+j]= parami['gamma']
    parami['alpha'+j]= parami['alpha']
    parami = parami[['rs','alpha'+j,'beta'+j,'gamma'+j]]
    param = param.merge(parami, on = ['rs'])

alphas = param.filter(regex='alpha')
alphas['rs']=param['rs']
alphas = alphas.assign(alpha_mean=alphas.mean(axis=1))
alphas = alphas[['rs','alpha_mean']]
betas = param.filter(regex='beta')
betas['rs']=param['rs']
betas = betas.assign(beta_mean=betas.mean(axis=1))
betas = betas[['rs','beta_mean']]
gammas = param.filter(regex='gamma')
gammas['rs']=param['rs']
gammas = alphas.assign(gamma_mean=gammas.mean(axis=1))
gammas = gammas[['rs','gamma_mean']]

mean_params = alphas.merge(betas, on = ['rs'])
mean_params = mean_params.merge(gammas, on = ['rs'])

mean_params = mean_params.join(mean_params['rs'].str.split(':', 1, expand=True).rename(columns={0:'CHROM', 1:'POS'}))
mean_params['start'] = pd.to_numeric(mean_params['POS'])
mean_params['stop'] = mean_params['start'] +1
mean_params['ID'] = 'id'
mean_params.head()
mean_params = mean_params[['CHROM', 'start','stop', 'ID','alpha_mean','beta_mean','gamma_mean']]
alphas = mean_params[['CHROM', 'start','stop', 'ID','alpha_mean']]
betas = mean_params[['CHROM', 'start','stop', 'ID','beta_mean']]
gammas = mean_params[['CHROM', 'start','stop', 'ID','gamma_mean']]

# output files with means

mean_params.to_csv(working_dir + 'mean_params_'+trait_name+'.bed',index=False, sep = "\t", header = False)
alphas.to_csv(working_dir + 'mean_alphas_'+trait_name+'.bed',index=False, sep = "\t", header = False)
betas.to_csv(working_dir + 'mean_betas_'+trait_name+'.bed',index=False, sep = "\t", header = False)
gammas.to_csv(working_dir + 'mean_gammas_'+trait_name+'.bed',index=False, sep = "\t", header = False)
