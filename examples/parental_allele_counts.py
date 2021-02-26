# This script is used after using the following GATK tools:
# 1. ReadBackedPhasing - to phase unbiased reads determined by WASP 
#		-generates phased.vcf
# 2. ASEReadCounter - to get read depth at heterozygous sites
#		-generates counts.csv
# 3. VariantsToTable - to output phased genotypes 
#		-generates snp_table.txt
  
# This script is used to create a dataframe of allele counts 
# at heterozygous sites in F1 hybrids that are alternatively 
# homozygous in their parents. User should make two dictionaries
# relating F1 samples to the name of their mother and father
# as they appear in the header line of the snp_table.txt

import pandas as pd
import numpy as np

#------------------------------------------------------------------------------------------
# set paths to directories / parent dictionaries / F1 sample names (infiles) / scaffold sizes file directory
#------------------------------------------------------------------------------------------

cts_dir = 'C:/Users/jmcgirr/Documents/remote_pups/ase/allele_counts/' # contains sample_counts.csv files
snp_dir ='C:/Users/jmcgirr/Documents/remote_pups/ase/'                # contains snp tables for F1 (rna) and parents (dna)
out_dir = 'C:/Users/jmcgirr/Documents/remote_pups/ase/parental_counts/dad_doublecheck/' 

mom_dict = {"CAE1": "CRPA1000", "CAE2": "CRPA1000", "CAE3": "CRPA1000", "CAE4": "CRPA1000", "CAE5": "CRPA1000", "CAT1": "CRPA1000", "CAT2": "CRPA1000", "CAT3": "CRPA1000", "CME1": "CRPM1000", "CME2": "CRPM1000", "CME5": "CRPM1000", "CMT1": "CRPM1000", "CMT2": "CRPM1000", "CMT3": "CRPM1000", "CPE1": "CRPP1000", "CPE2": "CRPP1000", "CPE3": "CRPP1000", "CPE4": "CRPP1000", "CPE5": "CRPP1000", "CPT1": "CRPP1000", "CPT2": "CRPP1000", "CPT3": "CRPP1000", "CPU1": "CRPA1000", "CPU3": "CRPA1000", "CPU5": "CRPA1000", "CQE1": "CUNP10", "CQE2": "CUNP10", "CQE3": "CUNP10", "CQT1": "CUNP10", "CQT2": "CUNP10", "CUT1": "CRPA1000", "CUT2": "CRPA1000", "CUT3": "CRPA1000", "CVE1": "CRPA1000", "CVE2": "CRPA1000", "CVE5": "CRPA1000", "CVT1": "CRPA1000", "CVT2": "CRPA1000", "CVT3": "CRPA1000", "CWE2": "CRPM1000", "CWE3": "CRPM1000", "CWE4": "CRPM1000", "CWT1": "CRPM1000", "CWT2": "CRPM1000", "CWT3": "CRPM1000", "CXE2": "CRPM1000", "CXE3": "CRPM1000", "CXE4": "CRPM1000", "CXT1": "CRPM1000", "CXT2": "CRPM1000", "CXT3": "CRPM1000", "LFE2": "CRPA1000", "LFE3": "CRPA1000", "LFE4": "CRPA1000", "LFT1": "CRPA1000", "LFT2": "CRPA1000", "LFT3": "CRPA1000", "LGE3": "CRPM1000", "LGE4": "CRPM1000", "LGT1": "CRPM1000", "LGT2": "CRPM1000", "LGT3": "CRPM1000", "LIE2": "OSPA1000", "LIE3": "OSPA1000", "LIE5": "OSPA1000", "LIT1": "OSPA1000", "LIT2": "OSPA1000", "LIT3": "OSPA1000", "LKE1": "OSPP1000", "LKE2": "OSPP1000", "LKE3": "OSPP1000", "LKT1": "OSPP1000", "LKT2": "OSPP1000", "LKT3": "OSPP1000", "NAE1": "NCCA1000", "NAE2": "NCCA1000", "NAE4": "NCCA1000", "NAT1": "NCCA1000", "NAT2": "NCCA1000", "NAT3": "NCCA1000", "NCE1": "NCCA1000", "NCE2": "NCCA1000", "NCE3": "NCCA1000", "NCE4": "NCCA1000", "NCE5": "NCCA1000", "NCT1": "NCCA1000", "NCT2": "NCCA1000", "NCT3": "NCCA1000", "OAE1": "OSPA1000", "OAE2": "OSPA1000", "OAE3": "OSPA1000", "OAE4": "OSPA1000", "OAT1": "OSPA1000", "OAT2": "OSPA1000", "OAT3": "OSPA1000", "OME1": "OSPM1000", "OME2": "OSPM1000", "OME3": "OSPM1000", "OME4": "OSPM1000", "OME5": "OSPM1000", "OMT1": "OSPM1000", "OMT2": "OSPM1000", "OMT3": "OSPM1000", "OPE1": "OSPP1000", "OPE2": "OSPP1000", "OPE3": "OSPP1000", "OPE4": "OSPP1000", "OPE5": "OSPP1000", "OPT1": "OSPP1000", "OPT2": "OSPP1000", "OPT3": "OSPP1000", "OUE1": "OSPA1000", "OUE3": "OSPA1000", "OUE4": "OSPA1000", "OUT1": "OSPA1000", "OUT2": "OSPA1000", "OUT3": "OSPA1000", "OVE1": "OSPA1000", "OVE4": "OSPA1000", "OVE5": "OSPA1000", "OVT1": "OSPA1000", "OVT2": "OSPA1000", "OVT3": "OSPA1000", "OXE2": "OSPM1000", "OXT1": "OSPM1000", "OXT2": "OSPM1000", "OXT3": "OSPM1000", "OYE1": "OSPP1000", "OYE2": "OSPP1000", "OYE3": "OSPP1000", "OYE4": "OSPP1000", "OYE5": "OSPP1000", "OYT1": "OSPP1000", "OYT2": "OSPP1000", "OYT3": "OSPP1000", "OZE2": "OSPP1000", "OZE4": "OSPP1000", "OZE5": "OSPP1000", "OZT1": "OSPP1000", "OZT2": "OSPP1000", "OZT3": "OSPP1000", "PAE1": "CUNP10", "PAE2": "CUNP10", "PAE5": "CUNP10", "PAT1": "CUNP10", "PAT2": "CUNP10", "PAT3": "CUNP10"}
dad_dict = {"CPU1": "CRPM1001","CPU3": "CRPM1001","CPU5": "CRPM1001","CUT1": "CRPM1001","CUT2": "CRPM1001","CUT3": "CRPM1001","CVE1": "CRPP1001","CVE2": "CRPP1001","CVE5": "CRPP1001","CVT1": "CRPP1001","CVT2": "CRPP1001","CVT3": "CRPP1001","CWE2": "CRPP1001","CWE3": "CRPP1001","CWE4": "CRPP1001","CWT1": "CRPP1001","CWT2": "CRPP1001","CWT3": "CRPP1001","CXE2": "CRPA1003","CXE3": "CRPA1003","CXE4": "CRPA1003","CXT1": "CRPA1003","CXT2": "CRPA1003","CXT3": "CRPA1003","LFE2": "OSPA1001","LFE3": "OSPA1001","LFE4": "OSPA1001","LFT1": "OSPA1001","LFT2": "OSPA1001","LFT3": "OSPA1001","LGE3": "OSPM1001","LGE4": "OSPM1001","LGT1": "OSPM1001","LGT2": "OSPM1001","LGT3": "OSPM1001","LIE2": "CRPA1001","LIE3": "CRPA1001","LIE5": "CRPA1001","LIT1": "CRPA1001","LIT2": "CRPA1001","LIT3": "CRPA1001","LKE1": "CRPP1001","LKE2": "CRPP1001","LKE3": "CRPP1001","LKT1": "CRPP1001","LKT2": "CRPP1001","LKT3": "CRPP1001","NAE1": "CRPA1001","NAE2": "CRPA1001","NAE4": "CRPA1001","NAT1": "CRPA1001","NAT2": "CRPA1001","NAT3": "CRPA1001","OUE1": "OSPM1001","OUE3": "OSPM1001","OUE4": "OSPM1001","OUT1": "OSPM1001","OUT2": "OSPM1001","OUT3": "OSPM1001","OVE1": "OSPP1001","OVE4": "OSPP1001","OVE5": "OSPP1001","OVT1": "OSPP1001","OVT2": "OSPP1001","OVT3": "OSPP1001","OXE2": "OSPA1001","OXT1": "OSPA1001","OXT2": "OSPA1001","OXT3": "OSPA1001","OYE1": "OSPA1001","OYE2": "OSPA1001","OYE3": "OSPA1001","OYE4": "OSPA1001","OYE5": "OSPA1001","OYT1": "OSPA1001","OYT2": "OSPA1001","OYT3": "OSPA1001","OZE2": "OSPM1001","OZE4": "OSPM1001","OZE5": "OSPM1001","OZT1": "OSPM1001","OZT2": "OSPM1001","OZT3": "OSPM1001","PAE1": "CRPA1001","PAE2": "CRPA1001","PAE5": "CRPA1001","PAT1": "CRPA1001","PAT2": "CRPA1001","PAT3": "CRPA1001"}

infiles = ["CPU1","CPU3","CPU5","CUT1","CUT2","CUT3","CVE1","CVE2","CVE5","CVT1","CVT2","CVT3","CWE2","CWE3","CWE4","CWT1","CWT2","CWT3","CXE2","CXE3","CXE4","CXT1","CXT2","CXT3","LFE2","LFE3","LFE4","LFT1","LFT2","LFT3","LGE3","LGE4","LGT1","LGT2","LGT3","LIE2","LIE3","LIE5","LIT1","LIT2","LIT3","LKE1","LKE2","LKE3","LKT1","LKT2","LKT3","NAE1","NAE2","NAE4","NAT1","NAT2","NAT3","OUE1","OUE3","OUE4","OUT1","OUT2","OUT3","OVE1","OVE4","OVE5","OVT1","OVT2","OVT3","OXE2","OXT1","OXT2","OXT3","OYE1","OYE2","OYE3","OYE4","OYE5","OYT1","OYT2","OYT3","OZE2","OZE4","OZE5","OZT1","OZT2","OZT3","PAE1","PAE2","PAE5","PAT1","PAT2","PAT3"]

all_rna_snps = pd.read_csv(snp_dir + "rna_snp_table.txt", sep ='\t')
all_dna_snps = pd.read_csv(snp_dir + "rna_parents_filtered_snps_passed.Q20.MAF0.05.MAXMISS0.5.recode_snp_table.txt", sep ='\t')
all_dna_snps['snpIndex'] = all_dna_snps['CHROM'].astype(str) + ':' + all_dna_snps['POS'].astype(str)


# dictionary to match rna samples (PCA outliers removed) with mother dna samples
# mom_dict = {"CAE1": "CAF1", "CAE2": "CAF1", "CAE3": "CAF1", "CAE4": "CAF1", "CAE5": "CAF1", "CAT1": "CAF1", "CAT2": "CAF1", "CAT3": "CAF1", "CME1": "CMF1", "CME2": "CMF1", "CME5": "CMF1", "CMT1": "CMF1", "CMT2": "CMF1", "CMT3": "CMF1", "CPE1": "CPF1", "CPE2": "CPF1", "CPE3": "CPF1", "CPE4": "CPF1", "CPE5": "CPF1", "CPT1": "CPF1", "CPT2": "CPF1", "CPT3": "CPF1", "CPU1": "CAF1", "CPU3": "CAF1", "CPU5": "CAF1", "CQE1": "CUNP", "CQE2": "CUNP", "CQE3": "CUNP", "CQT1": "CUNP", "CQT2": "CUNP", "CUT1": "CAF1", "CUT2": "CAF1", "CUT3": "CAF1", "CVE1": "CAF1", "CVE2": "CAF1", "CVE5": "CAF1", "CVT1": "CAF1", "CVT2": "CAF1", "CVT3": "CAF1", "CWE2": "CMF1", "CWE3": "CMF1", "CWE4": "CMF1", "CWT1": "CMF1", "CWT2": "CMF1", "CWT3": "CMF1", "CXE2": "CMF1", "CXE3": "CMF1", "CXE4": "CMF1", "CXT1": "CMF1", "CXT2": "CMF1", "CXT3": "CMF1", "LFE2": "CAF1", "LFE3": "CAF1", "LFE4": "CAF1", "LFT1": "CAF1", "LFT2": "CAF1", "LFT3": "CAF1", "LGE3": "CMF1", "LGE4": "CMF1", "LGT1": "CMF1", "LGT2": "CMF1", "LGT3": "CMF1", "LIE2": "OAF1", "LIE3": "OAF1", "LIE5": "OAF1", "LIT1": "OAF1", "LIT2": "OAF1", "LIT3": "OAF1", "LKE1": "OPF1", "LKE2": "OPF1", "LKE3": "OPF1", "LKT1": "OPF1", "LKT2": "OPF1", "LKT3": "OPF1", "NAE1": "NAF1", "NAE2": "NAF1", "NAE4": "NAF1", "NAT1": "NAF1", "NAT2": "NAF1", "NAT3": "NAF1", "NCE1": "NAF1", "NCE2": "NAF1", "NCE3": "NAF1", "NCE4": "NAF1", "NCE5": "NAF1", "NCT1": "NAF1", "NCT2": "NAF1", "NCT3": "NAF1", "OAE1": "OAF1", "OAE2": "OAF1", "OAE3": "OAF1", "OAE4": "OAF1", "OAT1": "OAF1", "OAT2": "OAF1", "OAT3": "OAF1", "OME1": "OMF1", "OME2": "OMF1", "OME3": "OMF1", "OME4": "OMF1", "OME5": "OMF1", "OMT1": "OMF1", "OMT2": "OMF1", "OMT3": "OMF1", "OPE1": "OPF1", "OPE2": "OPF1", "OPE3": "OPF1", "OPE4": "OPF1", "OPE5": "OPF1", "OPT1": "OPF1", "OPT2": "OPF1", "OPT3": "OPF1", "OUE1": "OAF1", "OUE3": "OAF1", "OUE4": "OAF1", "OUT1": "OAF1", "OUT2": "OAF1", "OUT3": "OAF1", "OVE1": "OAF1", "OVE4": "OAF1", "OVE5": "OAF1", "OVT1": "OAF1", "OVT2": "OAF1", "OVT3": "OAF1", "OXE2": "OMF1", "OXT1": "OMF1", "OXT2": "OMF1", "OXT3": "OMF1", "OYE1": "OPF1", "OYE2": "OPF1", "OYE3": "OPF1", "OYE4": "OPF1", "OYE5": "OPF1", "OYT1": "OPF1", "OYT2": "OPF1", "OYT3": "OPF1", "OZE2": "OPF1", "OZE4": "OPF1", "OZE5": "OPF1", "OZT1": "OPF1", "OZT2": "OPF1", "OZT3": "OPF1", "PAE1": "CUNP", "PAE2": "CUNP", "PAE5": "CUNP", "PAT1": "CUNP", "PAT2": "CUNP", "PAT3": "CUNP"}
# all samples
# infiles = ["CAE1","CAE2","CAE3","CAE4","CAE5","CAT1","CAT2","CAT3","CME1","CME2","CME5","CMT1","CMT2","CMT3","CPE1","CPE2","CPE3","CPE4","CPE5","CPT1","CPT2","CPT3","CPU1","CPU3","CPU5","CQE1","CQE2","CQE3","CQT1","CQT2","CUT1","CUT2","CUT3","CVE1","CVE2","CVE5","CVT1","CVT2","CVT3","CWE2","CWE3","CWE4","CWT1","CWT2","CWT3","CXE2","CXE3","CXE4","CXT1","CXT2","CXT3","LFE2","LFE3","LFE4","LFT1","LFT2","LFT3","LGE3","LGE4","LGT1","LGT2","LGT3","LIE2","LIE3","LIE5","LIT1","LIT2","LIT3","LKE1","LKE2","LKE3","LKT1","LKT2","LKT3","NAE1","NAE2","NAE4","NAT1","NAT2","NAT3","NCE1","NCE2","NCE3","NCE4","NCE5","NCT1","NCT2","NCT3","OAE1","OAE2","OAE3","OAE4","OME1","OME2","OME3","OME4","OME5","OMT1","OMT2","OMT3","OPE1","OPE2","OPE3","OPE4","OPE5","OPT1","OPT2","OPT3","OUE1","OUE3","OUE4","OUT1","OUT2","OUT3","OVE1","OVE4","OVE5","OVT1","OVT2","OVT3","OXE2","OXT1","OXT2","OXT3","OYE1","OYE2","OYE3","OYE4","OYE5","OYT1","OYT2","OYT3","OZE2","OZE4","OZE5","OZT1","OZT2","OZT3","PAE1","PAE2","PAE5","PAT1","PAT2","PAT3"]
# dad double check to confirm alternatively fixed homozygous (only for hybrids)


for infile in infiles:

    cts = pd.read_csv(cts_dir +infile + "_counts.csv", sep = '\t')
    cts['snpIndex'] = cts['contig'].astype(str) + ':'+ cts['position'].astype(str) 
    cts = cts[['snpIndex','refAllele','altAllele','refCount', 'altCount', 'totalCount']]
    mom = mom_dict[infile] +".GT"
    snps = all_dna_snps[[mom,'snpIndex']]
    
    mom_kid = cts.merge(snps, on='snpIndex')
    mom_kid = mom_kid.join(mom_kid[mom].str.split('/', 1, expand=True).rename(columns={0:'momAllele', 1:'a2'}))
    #only analyze homozygous alleles in mom that are het in offspring
    mom_kid = mom_kid[mom_kid['momAllele'] == mom_kid['a2']]
    mom_kid = mom_kid[(mom_kid['momAllele'] == mom_kid['refAllele']) | (mom_kid['momAllele'] == mom_kid['altAllele'])]
    # set minimum coverage at site
    #mom_kid = mom_kid[(mom_kid['refCount'] >= 10) & (mom_kid['altCount'] >= 10)]
    mom_kid.loc[mom_kid['refAllele'] == mom_kid['momAllele'], 'momCount'] = mom_kid['refCount']
    mom_kid.loc[mom_kid['altAllele'] == mom_kid['momAllele'], 'momCount'] = mom_kid['altCount']
    mom_kid = mom_kid.join(mom_kid['snpIndex'].str.split(':', 1, expand=True).rename(columns={0:'chrom', 1:'position'}))
    mom_kid = mom_kid[['chrom','position','snpIndex','refAllele','altAllele','refCount','altCount','totalCount','momAllele', 'momCount']]
    
    dad = dad_dict[infile] +".GT"
    snps_dad = all_dna_snps[[dad,'snpIndex']]
    
    mom_kid['dadCount'] = mom_kid['totalCount'] - mom_kid['momCount']
    mom_kid['momAllele_is_refAllele'] = np.where(mom_kid['refAllele'] == mom_kid['momAllele'], 'yes', 'no')
    mom_kid['momAllele_is_majorAllele'] = np.where(mom_kid['momCount'] > mom_kid['dadCount'], 'yes', 'no')
    
    # dad allele doublecheck
    dad = dad_dict[infile] +".GT"
    snps_dad = all_dna_snps[[dad,'snpIndex']]
    snps_dad = snps_dad.join(snps_dad[dad].str.split('/', 1, expand=True).rename(columns={0:'dadAllelereal1', 1:'dadAllelereal2'}))
    snps_dad = snps_dad[snps_dad['dadAllelereal1'] == snps_dad['dadAllelereal2']]
    mom_kid = mom_kid.merge(snps_dad, on='snpIndex')
    mom_kid = mom_kid[(mom_kid['momAllele'] != mom_kid['dadAllelereal1']) & (mom_kid['momAllele'] != mom_kid['dadAllelereal2'])]

    
    mom_kid.to_csv(out_dir +infile+'_parental_counts.txt',index=False, sep = "\t")
