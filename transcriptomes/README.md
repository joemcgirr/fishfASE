# Transcriptomic Analyses
## transcriptomes.ipynb jupyter notebook generates scripts to:
1. trim and align reads with STAR
2. call snps with gatk 3.8
3. remove biased reads with WASP
4. phase reads
5. calculate quality metrics with RseQC
6. create gene feature file with stringtie
7. count reads with featureCounts
8. run mbased for allele specific expression


# Commands
## 1. unzip, trim, and align reads
```
gzip -d sample.fq.gz
trim_galore -q 20 --paired --illumina sample_R1.fq sample_R2.fq
```
> create genome directory with [STAR](https://github.com/alexdobin/STAR)
```
STAR --runThreadN 4 \
--runMode genomeGenerate 
--genomeDir /path/to/genome_dir \
--genomeFastaFiles /path/to/reference.fasta \
--sjdbOverhang 149 \
```
> align with STAR
```
STAR --runThreadN 4 \
--genomeDir /path/to/genome_dir \
--readFilesIn sample_trimmed_R1.fq sample_trimmed_R2.fq \
--outFileNamePrefix sample \
--outSAMtype BAM SortedByCoordinate
```
> process `.bam`
```
samtools view -Sb sample_Aligned.out.sam > sample.bam 
samtools sort sample.bam  -o sample.sort.bam 
samtools index sample.sort.bam 
```
## 2. call snps with gatk 3.8
```
java -jar picard.jar AddOrReplaceReadGroups I=sample.sort.bam O=sample.sort.RG.bam RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=sample
samtools index sample.sort.RG.bam
gatk -T SplitNCigarReads -R reference.fasta -I sample.sort.RG.bam -o sample.sort.RG.split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS
samtools index sample.sort.RG.split.bam
gatk -T HaplotypeCaller -ERC GVCF -drf DuplicateRead -R reference.fasta -I sample.sort.RG.split.bam -dontUseSoftClippedBases -stand_call_conf 20.0 -nct 4 -o sample_raw_variants.g.vcf
bgzip sample_raw_variants.g.vcf
tabix -p vcf sample_raw_variants.g.vcf.gz
vcf-merge sample1_raw_variants.g.vcf.gz sample2_raw_variants.g.vcf.gz | bgzip -c > merged.vcf.gz
```

## 3. remove biased reads with WASP
> see WASP documentation at [bmvdgeijn/WASP](https://github.com/bmvdgeijn/WASP) for details
>
> STAR also added a flag implementing WASP filters
```
snp2h5 --chrom reference.scaffold.sizes --format vcf --haplotype haplotypes.h5 --snp_index snp_index.h5 --snp_tab snp_tab.h5 /pine/scr/j/m/jmcgirr/pupfish_transcriptomes/vcf/all_rna_filtered_snps_pre_wasp.reheader.vcf \n')

```




