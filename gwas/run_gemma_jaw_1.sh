#!/bin/bash

#SBATCH --job-name=script_run_gemma_low_jaw1
#SBATCH --mem=16G
#SBATCH --ntasks=4
#SBATCH --time=72:00:00 ##
#SBATCH -e script_run_gemma_low_jaw1_%A_%a.err ##error file if you want one (will be per job)
#SBATCH --mail-user=jmcgirr@email.unc.edu ##email you when job starts,ends,etc
#SBATCH --mail-type=ALL

/proj/cmarlab/users/joe/gemma/gemma-0.98.1-linux-static -bfile /pine/scr/j/m/jmcgirr/gemma_emilie/all_pupfish_filtered_snps_passed.Q20.MAF0.05.MAXMISS0.5.plink -w 50000000 -s 100000000 -n 3 -rpace 10000 -wpace 100000 -bslmm 1 -o all_pupfish_filtered_snps_passed.Q20.MAF0.05.MAXMISS0.5.plink_sansal_run_1 

