#------------------------------------------------------------------------------------
#
# parental_allele_counts.py
#
# This script is used after using the following GATK tools:
# 1. ReadBackedPhasing - to phase unbiased reads determined by WASP 
#		-generates phased.vcf
# 2. ASEReadCounter - to get read depth at heterozygous sites
#		-generates counts.csv
# 3. VariantsToTable - to output phased genotypes 
#		-generates snp_table.txt
# 
# This script is used to create a data frame of allele counts 
# at heterozygous sites in F1 hybrids that are alternatively 
# homozygous in their parents. Make two dictionaries
# relating F1 samples to the name of their mother and father
# as they appear in the header line of the rna_snp_table.txt
#
# example inputs:
# https://github.com/joemcgirr/fishfASE/blob/master/examples/CUT1_counts.csv
# https://github.com/joemcgirr/fishfASE/blob/master/examples/rna_snp_table_example.txt
# https://github.com/joemcgirr/fishfASE/blob/master/examples/rna_parents_snp_table_example.txt
#
# example output:
# https://github.com/joemcgirr/fishfASE/blob/master/examples/CUT1_parental_counts.txt
#
#------------------------------------------------------------------------------------

import pandas as pd
import numpy as np

#------------------------------------------------------------------------------------
# set paths to directories 
#------------------------------------------------------------------------------------

cts_dir = 'C:/Users/jmcgirr/Documents/remote_pups/ase/allele_counts/'  # contains sample_counts.csv files
snp_dir = 'C:/Users/jmcgirr/Documents/remote_pups/ase/'                # contains snp tables for F1 (rna) and F0 parents (dna)
out_dir = 'C:/Users/jmcgirr/Documents/remote_pups/ase/parental_counts/' 

all_rna_snps = pd.read_csv(snp_dir + "rna_snp_table.txt", sep ='\t')
all_dna_snps = pd.read_csv(snp_dir + "rna_parents_snp_table.txt", sep ='\t')
all_dna_snps['snpIndex'] = all_dna_snps['CHROM'].astype(str) + ':' + all_dna_snps['POS'].astype(str)


# dictionary to match F1 rna samples with F0 maternal dna sample
mom_dict = {"CAE1": "CRPA1000", "CAE2": "CRPA1000", "CAE3": "CRPA1000", "CAE4": "CRPA1000", "CAE5": "CRPA1000", "CAT1": "CRPA1000", "CAT2": "CRPA1000", "CAT3": "CRPA1000", "CME1": "CRPM1000", "CME2": "CRPM1000", "CME5": "CRPM1000", "CMT1": "CRPM1000", "CMT2": "CRPM1000", "CMT3": "CRPM1000", "CPE1": "CRPP1000", "CPE2": "CRPP1000", "CPE3": "CRPP1000", "CPE4": "CRPP1000", "CPE5": "CRPP1000", "CPT1": "CRPP1000", "CPT2": "CRPP1000", "CPT3": "CRPP1000", "CPU1": "CRPA1000", "CPU3": "CRPA1000", "CPU5": "CRPA1000", "CQE1": "CUNP10", "CQE2": "CUNP10", "CQE3": "CUNP10", "CQT1": "CUNP10", "CQT2": "CUNP10", "CUT1": "CRPA1000", "CUT2": "CRPA1000", "CUT3": "CRPA1000", "CVE1": "CRPA1000", "CVE2": "CRPA1000", "CVE5": "CRPA1000", "CVT1": "CRPA1000", "CVT2": "CRPA1000", "CVT3": "CRPA1000", "CWE2": "CRPM1000", "CWE3": "CRPM1000", "CWE4": "CRPM1000", "CWT1": "CRPM1000", "CWT2": "CRPM1000", "CWT3": "CRPM1000", "CXE2": "CRPM1000", "CXE3": "CRPM1000", "CXE4": "CRPM1000", "CXT1": "CRPM1000", "CXT2": "CRPM1000", "CXT3": "CRPM1000", "LFE2": "CRPA1000", "LFE3": "CRPA1000", "LFE4": "CRPA1000", "LFT1": "CRPA1000", "LFT2": "CRPA1000", "LFT3": "CRPA1000", "LGE3": "CRPM1000", "LGE4": "CRPM1000", "LGT1": "CRPM1000", "LGT2": "CRPM1000", "LGT3": "CRPM1000", "LIE2": "OSPA1000", "LIE3": "OSPA1000", "LIE5": "OSPA1000", "LIT1": "OSPA1000", "LIT2": "OSPA1000", "LIT3": "OSPA1000", "LKE1": "OSPP1000", "LKE2": "OSPP1000", "LKE3": "OSPP1000", "LKT1": "OSPP1000", "LKT2": "OSPP1000", "LKT3": "OSPP1000", "NAE1": "NCCA1000", "NAE2": "NCCA1000", "NAE4": "NCCA1000", "NAT1": "NCCA1000", "NAT2": "NCCA1000", "NAT3": "NCCA1000", "NCE1": "NCCA1000", "NCE2": "NCCA1000", "NCE3": "NCCA1000", "NCE4": "NCCA1000", "NCE5": "NCCA1000", "NCT1": "NCCA1000", "NCT2": "NCCA1000", "NCT3": "NCCA1000", "OAE1": "OSPA1000", "OAE2": "OSPA1000", "OAE3": "OSPA1000", "OAE4": "OSPA1000", "OAT1": "OSPA1000", "OAT2": "OSPA1000", "OAT3": "OSPA1000", "OME1": "OSPM1000", "OME2": "OSPM1000", "OME3": "OSPM1000", "OME4": "OSPM1000", "OME5": "OSPM1000", "OMT1": "OSPM1000", "OMT2": "OSPM1000", "OMT3": "OSPM1000", "OPE1": "OSPP1000", "OPE2": "OSPP1000", "OPE3": "OSPP1000", "OPE4": "OSPP1000", "OPE5": "OSPP1000", "OPT1": "OSPP1000", "OPT2": "OSPP1000", "OPT3": "OSPP1000", "OUE1": "OSPA1000", "OUE3": "OSPA1000", "OUE4": "OSPA1000", "OUT1": "OSPA1000", "OUT2": "OSPA1000", "OUT3": "OSPA1000", "OVE1": "OSPA1000", "OVE4": "OSPA1000", "OVE5": "OSPA1000", "OVT1": "OSPA1000", "OVT2": "OSPA1000", "OVT3": "OSPA1000", "OXE2": "OSPM1000", "OXT1": "OSPM1000", "OXT2": "OSPM1000", "OXT3": "OSPM1000", "OYE1": "OSPP1000", "OYE2": "OSPP1000", "OYE3": "OSPP1000", "OYE4": "OSPP1000", "OYE5": "OSPP1000", "OYT1": "OSPP1000", "OYT2": "OSPP1000", "OYT3": "OSPP1000", "OZE2": "OSPP1000", "OZE4": "OSPP1000", "OZE5": "OSPP1000", "OZT1": "OSPP1000", "OZT2": "OSPP1000", "OZT3": "OSPP1000", "PAE1": "CUNP10", "PAE2": "CUNP10", "PAE5": "CUNP10", "PAT1": "CUNP10", "PAT2": "CUNP10", "PAT3": "CUNP10"}

# dictionary to match F1 rna samples with F0 maternal dna sample
dad_dict = {"CPU1": "CRPM1001","CPU3": "CRPM1001","CPU5": "CRPM1001","CUT1": "CRPM1001","CUT2": "CRPM1001","CUT3": "CRPM1001","CVE1": "CRPP1001","CVE2": "CRPP1001","CVE5": "CRPP1001","CVT1": "CRPP1001","CVT2": "CRPP1001","CVT3": "CRPP1001","CWE2": "CRPP1001","CWE3": "CRPP1001","CWE4": "CRPP1001","CWT1": "CRPP1001","CWT2": "CRPP1001","CWT3": "CRPP1001","CXE2": "CRPA1003","CXE3": "CRPA1003","CXE4": "CRPA1003","CXT1": "CRPA1003","CXT2": "CRPA1003","CXT3": "CRPA1003","LFE2": "OSPA1001","LFE3": "OSPA1001","LFE4": "OSPA1001","LFT1": "OSPA1001","LFT2": "OSPA1001","LFT3": "OSPA1001","LGE3": "OSPM1001","LGE4": "OSPM1001","LGT1": "OSPM1001","LGT2": "OSPM1001","LGT3": "OSPM1001","LIE2": "CRPA1001","LIE3": "CRPA1001","LIE5": "CRPA1001","LIT1": "CRPA1001","LIT2": "CRPA1001","LIT3": "CRPA1001","LKE1": "CRPP1001","LKE2": "CRPP1001","LKE3": "CRPP1001","LKT1": "CRPP1001","LKT2": "CRPP1001","LKT3": "CRPP1001","NAE1": "CRPA1001","NAE2": "CRPA1001","NAE4": "CRPA1001","NAT1": "CRPA1001","NAT2": "CRPA1001","NAT3": "CRPA1001","OUE1": "OSPM1001","OUE3": "OSPM1001","OUE4": "OSPM1001","OUT1": "OSPM1001","OUT2": "OSPM1001","OUT3": "OSPM1001","OVE1": "OSPP1001","OVE4": "OSPP1001","OVE5": "OSPP1001","OVT1": "OSPP1001","OVT2": "OSPP1001","OVT3": "OSPP1001","OXE2": "OSPA1001","OXT1": "OSPA1001","OXT2": "OSPA1001","OXT3": "OSPA1001","OYE1": "OSPA1001","OYE2": "OSPA1001","OYE3": "OSPA1001","OYE4": "OSPA1001","OYE5": "OSPA1001","OYT1": "OSPA1001","OYT2": "OSPA1001","OYT3": "OSPA1001","OZE2": "OSPM1001","OZE4": "OSPM1001","OZE5": "OSPM1001","OZT1": "OSPM1001","OZT2": "OSPM1001","OZT3": "OSPM1001","PAE1": "CRPA1001","PAE2": "CRPA1001","PAE5": "CRPA1001","PAT1": "CRPA1001","PAT2": "CRPA1001","PAT3": "CRPA1001"}

# F1 rna samples
infiles = ["CPU1","CPU3","CPU5","CUT1","CUT2","CUT3","CVE1","CVE2","CVE5","CVT1","CVT2","CVT3","CWE2","CWE3","CWE4","CWT1","CWT2","CWT3","CXE2","CXE3","CXE4","CXT1","CXT2","CXT3","LFE2","LFE3","LFE4","LFT1","LFT2","LFT3","LGE3","LGE4","LGT1","LGT2","LGT3","LIE2","LIE3","LIE5","LIT1","LIT2","LIT3","LKE1","LKE2","LKE3","LKT1","LKT2","LKT3","NAE1","NAE2","NAE4","NAT1","NAT2","NAT3","OUE1","OUE3","OUE4","OUT1","OUT2","OUT3","OVE1","OVE4","OVE5","OVT1","OVT2","OVT3","OXE2","OXT1","OXT2","OXT3","OYE1","OYE2","OYE3","OYE4","OYE5","OYT1","OYT2","OYT3","OZE2","OZE4","OZE5","OZT1","OZT2","OZT3","PAE1","PAE2","PAE5","PAT1","PAT2","PAT3"]

for infile in infiles:

    cts = pd.read_csv(cts_dir +infile + "_counts.csv", sep = '\t')
    cts['snpIndex'] = cts['contig'].astype(str) + ':'+ cts['position'].astype(str) 
    cts = cts[['snpIndex','refAllele','altAllele','refCount', 'altCount', 'totalCount']]
    mom = mom_dict[infile] +".GT"
    snps = all_dna_snps[[mom,'snpIndex']]
    
    mom_kid = cts.merge(snps, on='snpIndex')
    mom_kid = mom_kid.join(mom_kid[mom].str.split('/', 1, expand=True).rename(columns={0:'momAllele', 1:'a2'}))
    
    # only analyze homozygous alleles in maternal samples that are heterozygous in offspring
    mom_kid = mom_kid[mom_kid['momAllele'] == mom_kid['a2']]
    mom_kid = mom_kid[(mom_kid['momAllele'] == mom_kid['refAllele']) | (mom_kid['momAllele'] == mom_kid['altAllele'])]
    
    # set minimum coverage at site (>= 10 counts at each site)
    mom_kid = mom_kid[(mom_kid['refCount'] >= 10) & (mom_kid['altCount'] >= 10)]
    mom_kid.loc[mom_kid['refAllele'] == mom_kid['momAllele'], 'momCount'] = mom_kid['refCount']
    mom_kid.loc[mom_kid['altAllele'] == mom_kid['momAllele'], 'momCount'] = mom_kid['altCount']
    mom_kid = mom_kid.join(mom_kid['snpIndex'].str.split(':', 1, expand=True).rename(columns={0:'chrom', 1:'position'}))
    mom_kid = mom_kid[['chrom','position','snpIndex','refAllele','altAllele','refCount','altCount','totalCount','momAllele', 'momCount']]
    
    mom_kid['dadCount'] = mom_kid['totalCount'] - mom_kid['momCount']
    mom_kid['momAllele_is_refAllele'] = np.where(mom_kid['refAllele'] == mom_kid['momAllele'], 'yes', 'no')
    mom_kid['momAllele_is_majorAllele'] = np.where(mom_kid['momCount'] > mom_kid['dadCount'], 'yes', 'no')
    
    # confirm that the paternal allele is different from maternal allele
    dad = dad_dict[infile] +".GT"
    snps_dad = all_dna_snps[[dad,'snpIndex']]
    snps_dad = snps_dad.join(snps_dad[dad].str.split('/', 1, expand=True).rename(columns={0:'dadAllelereal1', 1:'dadAllelereal2'}))
    snps_dad = snps_dad[snps_dad['dadAllelereal1'] == snps_dad['dadAllelereal2']]
    mom_kid = mom_kid.merge(snps_dad, on='snpIndex')
    mom_kid = mom_kid[(mom_kid['momAllele'] != mom_kid['dadAllelereal1']) & (mom_kid['momAllele'] != mom_kid['dadAllelereal2'])]

    
    mom_kid.to_csv(out_dir +infile+'_parental_counts.txt',index=False, sep = "\t")