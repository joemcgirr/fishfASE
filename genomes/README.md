# Genomic Analyses
## genomes.ipynb jupyter notebook generates scripts to:
1. trim and align reads
2. deduplicate .bam files
3. call snps with gatk 3.8
4. calculate Fst, Tajima's D, and pi with vcftools
5. calculate Dxy with simon martin scripts
6. find hard sweeps with SweeD
7. create RaxML plylogeny
8. GWAS with GEMMA

# Commands
## 1. unzip, trim, and align reads
```
gzip -d sample.fq.gz
trim_galore -q 20 --paired --illumina sample_R1.fq sample_R2.fq
bwa mem -aM -t 4 -R "@RG\\tID:group1\\tSM:'+infile+'\\tPL:illumina\\tLB:lib1" reference.fasta sample_trim_R1.fq sample_trim_R2.fq > sample.sam
samtools view -Shu sample.sam > sample.bam
samtools index sample.bam
samtools sort sample.bam -o sample.sort.bam
samtools index sample.sort.bam
rm sample.sam
rm sample.bam
```

## 2. deduplicate .bam files with picard.jar
```
java -Xmx10g -jar picard.jar MarkDuplicates INPUT=sample.sort.bam OUTPUT=sample.sort.dedup.bam METRICS_FILE=sample.metrics.txt MAX_FILE_HANDLES=1000
samtools index sample.sort.dedup.bam
```
## 3. call snps with gatk 3.8
```
gatk -T HaplotypeCaller -ERC GVCF -drf DuplicateRead -R reference.fasta -I sample.sort.dedup.bam -dontUseSoftClippedBases -stand_call_conf 20.0 -nct 4 -o sample_raw_variants.g.vcf
gatk -T GenotypeGVCFs -R reference.fasta --variant sample1_raw_variants.g.vcf --variant sample2_raw_variants.g.vcf -o merged_raw_variants.vcf
vcftools --vcf merged_raw_variants.vcf --maf 0.05 --max-missing 0.9 --remove-indels --recode --recode-INFO-all --out filtered_snps.vcf
```
## 4. calculate Fst, Tajima's D, and pi with vcftools
```
vcftools --vcf filtered_snps.vcf --keep populations_1_and_2.txt --out population_1_vs_2.weir.fst --weir-fst-pop population_1.txt --weir-fst-pop population_2.txt
vcftools --vcf filtered_snps.vcf --TajimaD 20000 --out taj_d_20kb_windows.txt 
vcftools --vcf filtered_snps.vcf --window-pi 20000 --out pi_20kb_windows.txt
```
> bash commands to remove negative values from Fst output and calculate genome-wide mean Fst
```
sed \'s/-[0-9].*/0/g\' population_1_vs_2.weir.fst | sed \'s/-nan/0/g\' > population_1_vs_2.weir.fst
awk -F\'\\t\' \'{ sum += $3 } END { print sum / NR }\' population_1_vs_2.weir.fst > genome_wide_avg.txt
```
## 5. calculate Dxy with simon martin scripts
> see [simonhmartin/genomics_general](https://github.com/simonhmartin/genomics_general/tree/master/VCF_processing) for `parseVCF.py` and `popgenWindows.py`
```
python parseVCF.py -i input.vcf.gz --skipIndels --minQual 30 --gtf flag=DP min=5 | bgzip > output.geno.gz
python popgenWindows.py -w 10000 -m 10 -g output.geno.gz -o popgen_stats.csv -f phased -T 4 -p population_1 sample_1 sample_2 -p population_2.txt sample_3 sample_4
```
## 6. find hard sweeps with SweeD
> see [SweeD documentation](https://cme.h-its.org/exelixis/resource/download/software/sweed3.0_manual.pdf) and [Pavlidis et al 2013](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3748355/) for details
>
> -eN flag used to specify demography specific to San Salvador pupfish system (population decrease 10kya)
```
SweeD -input filtered_snps.vcf -s 100 -name sweeps.txt -eN .005 .01 -strictPolymorphic -folded -grid 1000 -osfs sweeps_sfs.txt
```
## 7. create RaxML plylogeny
> see [RaxML documentation](https://cme.h-its.org/exelixis/resource/download/NewManual.pdf) for details
>
> first create biallelic `.vcf` and convert to `.phy` with `vcf2phylip.py` from [edgardomortiz](https://github.com/edgardomortiz/vcf2phylip)
```
vcftools --vcf filtered_snps.vcf --min-alleles 2 --man-alleles 2 --recode --out biallelic_snps.vcf

python vcf2phylip.py --input biallelic_snps.vcf

standard-RAxML-master/raxmlHPC -x 76345 -f a -m GTRGAMMA -o outgroup_sample_name -p 28393 -s biallelic_snps.phy -# 10 -n output.name
```
## 8. GWAS with GEMMA
> see [GEMMA documentation](https://www.xzlab.org/software/GEMMAmanual.pdf) for details
>
> first convert biallelic `.vcf` to [plink](http://zzz.bwh.harvard.edu/plink/) format
```
vcftools --vcf biallelic_snps.vcf --plink --out biallelic_snps.plink
plink.exe --file biallelic_snps.plink --make-bed --out biallelic_snps.plink
gwas/gemma-0.98.1-linux-static -bfile biallelic_snps.plink -w 50000000 -s 100000000 -n 5 -rpace 10000 -wpace 100000 -bslmm 1 -o run_1.output
```




## Species matched to sample name
San Salvador Island Generalists
["RHPA1","SPPA1","CRPA3","MRKA1","PIGA3","WDPA1","OSPA4","ME2A2","OSPA7","OYSA1","PIGA1","OSPA6","OSPA5","CRPA1","MERA2","GNYA1","LILA1","ME2A1","GREA1","OSPA1","OSPA9","OYSA2","CLRA1","OSPA11","NLLA1","OSPA13","GREA2","OSPA12","OSPA8","PAIA1","MERA3","OSPA10"]
### San Salvador Island Molluscivores
["OYSM8","CRPM8","CRPM7","OSPM6","OSPM7","LILM5","OYSM6","OSPM5","CRPM3","CRPM6","MRKM5","OYSM1","CRPM9","CRPM2","LILM4","CRPM11","OSPM11","CRPM10","OSPM3","OSPM9","OYSM3","MRKM3","OYSM2","OSPM2","OYSM5","WDPM2","LILMQ","OSPM4","OSPM10","OSPM8","LILM3","OYSM4","OSPM1","OYSM7","MRKM1","MRKM4","MRKM2","CRPM1","CRPM5"]
### San Salvador Island Scale-eaters

["ME2P1","CRPP3","CRPP7","OSPP9","LILP4","CRPP8","OSPP1","LILP3","OYSP6","OSPP8","LILPQ","OYSP7","OYSP1","CRPP5","OSPP3","OSPP11","CRPP9","CRPP2","OSPP5","CRPP4","CRPPQ","LILP5","OSPP6","OSPP7","OSPP10","OSPP2","OSPP4","OYSP3","OYSP4"]
### San Salvador Island Small-jawed Scale-eaters

["OSPS8","OSPS2","GRES4","OYSS3","GREP1","OSPS11","OSPS5","LILS1","GRES3","OYSS5","GREP2","OSPS3","LILS3","LILS4","MERP1","LILS2","OSPS7","OSPS10","OSPS1","OSPS9","OYSS4","OYSS6","MERP2","OYSS9","OSPS6","OSPS4"]
### Caribbean Generalists and Outgroups
["MEGQ1","GEO2A8","CATA1","CUNP5","BAVA8","ARTA2","GEO2A5","VENA1","CURA21","VENA12","GEO2A9","NCCA11","NCCA4","GEOA2","NCCA2","NCCA5","NCCA12","CUNP7","NCCA1","NCCA15","GEOA10","BAVA11","CUNA6","NCCA9","BAVA2","FCTA1","BAVA4","VENA5","CUNA1","GEOA6","VENA10","MAY1","BAVA6","CUNP4","VENA2","BAVA5","CUNA2","BAVA10","CUNA7","VENA3","GEO2A1","CUNA10","EXUA2","NBIA1","BAVA14","CUNA4","FLSA1","CUNP11","MAFA1","VENA9","CUNA3","VENA7","PWLA1","ETA1","GEOA7","GEO2A6","BAVA42","BAVA7","VENA8","BAVA9","GEO2A7","NCCA8","SIM1","GEOA1","ARTA1","EPLA1","NCCA3","ACKA1","GEOA5","NCCA10","VENA13","GEOA4","CUNA9","KILA1","GEOA11","BAVA13","VENA4","CAIA1","BAVA41","CUNP3","LGIA1","EXUA1","CUNA8","DEAA1","BAVA12","CURA1","CUNP6","CURA2","BONA1","CUNA5","NCCA7","GEO2A3","GEO2A10","GEO2A4","SCLA1","GEO2A2","SALA1","NCCA6"]
fdir = "/pine/scr/j/m/jmcgirr/pupfish_genomes/Caribbean_pups/"
### San Salvador Island Breeders Used to Generate F1s for RNAseq
["CRPM1001","CRPP1000","CRPM1000","NCCA1000","CRPP1001","CRPA1000","OSPM1001","CRPA1001","CUNP10.2","OSPA1001","OSPP1000","OSPA1000","CRPA1003","OSPM1000","OSPP1001"]
