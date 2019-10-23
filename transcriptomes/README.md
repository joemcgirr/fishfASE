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
vcf-merge sample1_raw_variants.g.vcf.gz sample2_raw_variants.g.vcf.gz | bgzip -c > merged_raw.vcf.gz
vcftools --vcf merged_raw_variants.vcf --maf 0.05 --max-missing 0.9 --remove-indels --recode --out filtered_snps.vcf
```

## 3. remove biased reads with WASP
> see WASP documentation at [bmvdgeijn/WASP](https://github.com/bmvdgeijn/WASP) for details
>
> STAR also added a flag implementing WASP filters
```
snp2h5 --chrom reference.scaffold.sizes --format vcf --haplotype haplotypes.h5 --snp_index snp_index.h5 --snp_tab snp_tab.h5 merged.vcf
python3 find_intersecting_snps.py --is_paired_end --is_sorted --output_dir out_dir --snp_tab snp_tab.h5 --snp_index snp_index.h5 --haplotype haplotypes.h5 --samples sample.txt sample.sort.RG.split.bam \n')
gzip -d sample.sort.remap.fq1.gz
gzip -d sample.sort.remap.fq2.gz
```
> remap `.fq1` and `.fq2` with STAR (see step 1)
```
python3 /nas/longleaf/apps/wasp/2018-07/WASP/mapping/filter_remapped_reads.py sample.sort.to.remap.bam sample.remapped.sort.bam sample.keep.to.merge.bam
samtools merge sample.keep.merged.bam sample.keep.to.merge.bam sample.sort.keep.bam
samtools sort -o sample.filtered.merged.sort.bam sample.keep.merged.bam
samtools index sample.filtered.merged.sort.bam
```
> repeat step 2 with unbiased `.bam` files (call snps again) 

## 4. phase reads
```
gatk -T ReadBackedPhasing -R reference.fasta -I sample.filtered.merged.sort.bam --variant wasp_unbiased_filtered_snps.vcf -o /pine/scr/j/m/jmcgirr/pupfish_transcriptomes/wasp/vcf/'+infile+'_wasp_unbiased_phased.vcf --phaseQualityThresh 20.0 \n')
    outfile.write('gatk -T ASEReadCounter -R /proj/cmarlab/users/joe/Cyprinodon/bronto/asm.racon.fasta -U ALLOW_SEQ_DICT_INCOMPATIBILITY -o '+counts_dir+infile+'_counts.csv -I '+rg_bams_dir+infile+'_filtered.merged.sort.RG.split.bam -sites /pine/scr/j/m/jmcgirr/pupfish_transcriptomes/wasp/vcf/'+infile+'_wasp_unbiased_phased.vcf \n')           
    outfile.write('#gatk -T VariantsToTable -R /proj/cmarlab/users/joe/Cyprinodon/bronto/asm.racon.fasta -V /pine/scr/j/m/jmcgirr/pupfish_transcriptomes/wasp/vcf/'+infile+'_wasp_unbiased_phased.vcf -F CHROM -F POS -GF GT -GF HP -o '+counts_dir+infile+'_snp_table.txt')

```


