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
> create genome directory with STAR
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
--outFileNamePrefix sample_name \
--outSAMtype BAM SortedByCoordinate
```
