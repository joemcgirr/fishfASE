
### SKELLY  ###

### Skelly Drafts ###
###############
library(reshape2)
  #### RNA ####

setwd("c:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/axm/ASE/skelly/")
gatk_files <- c("axmE1",
                "axmE2",
                "axmJ1",
                "axmJ4",
                "axmJ5",
                "axmJ6",
                "CPAE1",
                "CPAE2",
                "CPAE3",
                "CPAJ1",
                "CPAJ2",
                "CPAJ3",
                "CPME1",
                "CPME2",
                "CPME3",
                "CPMJ1",
                "CPMJ2",
                "CPMJ3",
                "CPPE1",
                "CPPE2",
                "CPPE3",
                "CPPJ1",
                "CPPJ2",
                "LLAE1",
                "LLAE2",
                "LLAE3",
                "LLAJ1",
                "LLAJ2",
                "LLAJ3",
                "LLME1",
                "LLME2",
                "LLME3",
                "LLMJ1",
                "LLMJ2",
                "LLMJ3")
# make haplotype counts file
for(file in gatk_files) 
{
  
  
  infile1 <- paste(file,"_snp_table_haplotypes.file",sep="")
  infile2 <- paste(file,"_allele_counts_phased_gatk.table",sep="")
  outfile <- paste(file,"_haplotype_counts.txt",sep="")
  
  #infile1 <- "axmE1_snp_table_haplotypes.file"
  #infile2 <- "axmE1_allele_counts_phased_gatk.table"
  
  haps <- read.table(infile1, header = TRUE, stringsAsFactors = FALSE)
  head(haps)
  nrow(haps)
  haps$snpIndex <- paste(haps$CHROM, haps$POS, sep = ":")
  counts <- read.table(infile2, header = TRUE, stringsAsFactors = FALSE)
  head(counts)
  nrow(counts)
  counts$snpIndex <- paste(counts$contig, counts$position, sep = ":")
  haps_counts <- merge(counts, haps, by = c("snpIndex"))
  head(haps_counts)
  nrow(haps_counts)
  keeps <- c("snpIndex", "refAllele", "altAllele", "refCount", "altCount", "totalCount", "sample1.SNPs.haplotype_A", "sample1.SNPs.haplotype_B")
  haps_counts <- haps_counts[keeps]
  
  
  refs <- haps_counts$refAllele
  alts <- haps_counts$altAllele
  ref_counts <- haps_counts$refCount
  alt_counts <- haps_counts$altCount
  alleleOnes <- haps_counts$sample1.SNPs.haplotype_A
  alleleTwos <- haps_counts$sample1.SNPs.haplotype_B
  alleleOneCounts <- c()
  alleleTwoCounts <- c()
  for (i in c(1:nrow(haps_counts)))
  {
    
    if (refs[i] == alleleOnes[i] & alts[i] == alleleTwos[i])
    {
      alleleOneCount <- ref_counts[i]
      alleleTwoCount <- alt_counts[i]
    }
    else if (refs[i] == alleleTwos[i] & alts[i] == alleleOnes[i])
    {
      alleleOneCount <- alt_counts[i]
      alleleTwoCount <- ref_counts[i]
    }

    else {
      alleleOneCount <- "NA"
      alleleTwoCount <- "NA"
    }
    
    alleleOneCounts <- c(alleleOneCounts, alleleOneCount)
    alleleTwoCounts <- c(alleleTwoCounts, alleleTwoCount) 
  }
  
  haps_counts$alleleOneCount <- alleleOneCounts
  haps_counts$alleleTwoCount <- alleleTwoCounts
  
  keeps <- c("snpIndex", "alleleOneCount", "alleleTwoCount")
  haps_counts <- haps_counts[keeps]
  
  #write.table(haps_counts,quote=FALSE, row.names = FALSE, col.names = TRUE, outfile)
}

## merge genes with haplotype counts file (use python snp_genes_overlap_haplotype_counts)

  #### DNA ####

# make haplotype counts file
setwd("c:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/axm/ASE/skelly/")
gatk_files <- c("CRPA1_DNA",
                "CRPM1_DNA",
                "LILA1_DNA",
                "LILM3_DNA")
gatk_files <- c("LLA-QTL",
                "LLM-QTL")
  

for(file in gatk_files) 
{
  
  infile <- paste(file,"_snp_table_nw.txt",sep="")
  outfile <- paste(file,"_snp_table_haplotypes.file",sep="")


  dat <- read.table(infile, header = TRUE, stringsAsFactors = FALSE)
  head(dat)
  dat <- na.omit(dat)
  dat <- cbind(dat, colsplit(dat[,4], ",", c("HP1", "HP2")))
  dat <- cbind(dat, colsplit(dat[,5], "-", c("junk1", "HP1.1")))
  dat <- cbind(dat, colsplit(dat[,6], "-", c("junk2", "HP2.1")))
  dat <- cbind(dat, colsplit(dat[,3], "/", c("GT_1", "GT_2")))
  dat1 <- dat[which(dat$HP1.1 == 1 ),]
  dat2 <- dat[which(dat$HP1.1 == 2 ),]
  dat1["sample1.SNPs.haplotype_A"] <- dat1$GT_1
  dat2["sample1.SNPs.haplotype_A"] <- dat2$GT_2
  dat1["sample1.SNPs.haplotype_B"] <- dat1$GT_2
  dat2["sample1.SNPs.haplotype_B"] <- dat2$GT_1
  final <- rbind(dat1, dat2)
  tail(final,20)
  keeps <- c("CHROM", "POS", "sample1.SNPs.haplotype_A", "sample1.SNPs.haplotype_B")
  final <- final[keeps]
  write.table(final,quote=FALSE,sep="\t", row.names = FALSE, col.names = TRUE, outfile)
}


names <- read.table("C:/Users/Joseph McGirr Lab/Desktop/Cyprinodon/scaffold_names", header = TRUE, stringsAsFactors = FALSE)
head(names)
head(fixed)
keeps <- c("RefSeq_Accession.version", "GenBank_Accession.version")
scaffolds <- names[keeps]
head(scaffolds)
scaffolds$genomic_accession <- scaffolds$RefSeq_Accession.version
scaffolds$CHROM <- scaffolds$GenBank_Accession.version
keeps <- c("CHROM", "genomic_accession")
scaffolds_key <- scaffolds[keeps]
scaffolds_key$Chr <- scaffolds_key$genomic_accession
head(scaffolds_key)


for(file in gatk_files) 
{
  
  
  infile1 <- paste(file,"_snp_table_haplotypes.file",sep="")
  infile2 <- paste(file,"_allele_counts_gatk_nw.csv",sep="")
  outfile <- paste(file,"_haplotype_counts.txt",sep="")
  
  #infile1 <- "LLM-QTL_snp_table_haplotypes.file"
  #infile2 <- "LLM-QTL_allele_counts_gatk_nw.csv"
  
  haps <- read.table(infile1, header = TRUE, stringsAsFactors = FALSE)
  head(haps)
  nrow(haps)
  #haps <- merge(haps, scaffolds_key, by = c("CHROM"))
  #haps$CHROM <- haps$Chr
  
  haps$snpIndex <- paste(haps$CHROM, haps$POS, sep = ":")
  counts <- read.table(infile2, header = TRUE, stringsAsFactors = FALSE)
  head(counts)
  nrow(counts)
  counts$CHROM <- counts$contig
  #counts <- merge(counts, scaffolds_key, by = c("CHROM"))
  #counts$CHROM <- counts$Chr
  counts$snpIndex <- paste(counts$CHROM, counts$position, sep = ":")
  haps_counts <- merge(counts, haps, by = c("snpIndex"))
  head(haps_counts)
  nrow(haps_counts)
  keeps <- c("snpIndex", "refAllele", "altAllele", "refCount", "altCount", "totalCount", "sample1.SNPs.haplotype_A", "sample1.SNPs.haplotype_B")
  haps_counts <- haps_counts[keeps]
  
  
  refs <- haps_counts$refAllele
  alts <- haps_counts$altAllele
  ref_counts <- haps_counts$refCount
  alt_counts <- haps_counts$altCount
  alleleOnes <- haps_counts$sample1.SNPs.haplotype_A
  alleleTwos <- haps_counts$sample1.SNPs.haplotype_B
  alleleOneCounts <- c()
  alleleTwoCounts <- c()
  for (i in c(1:nrow(haps_counts)))
  {
    
    if (refs[i] == alleleOnes[i] & alts[i] == alleleTwos[i])
    {
      alleleOneCount <- ref_counts[i]
      alleleTwoCount <- alt_counts[i]
    }
    else if (refs[i] == alleleTwos[i] & alts[i] == alleleOnes[i])
    {
      alleleOneCount <- alt_counts[i]
      alleleTwoCount <- ref_counts[i]
    }
    
    else {
      alleleOneCount <- "NA"
      alleleTwoCount <- "NA"
    }
    
    alleleOneCounts <- c(alleleOneCounts, alleleOneCount)
    alleleTwoCounts <- c(alleleTwoCounts, alleleTwoCount) 
  }
  
  haps_counts$alleleOneCount <- alleleOneCounts
  haps_counts$alleleTwoCount <- alleleTwoCounts
  
  keeps <- c("snpIndex", "alleleOneCount", "alleleTwoCount")
  haps_counts <- haps_counts[keeps]
  
  write.table(haps_counts,quote=FALSE, row.names = FALSE, col.names = TRUE, outfile)
}

# python snp_genes_overlap_haplotype_counts


### set minimum coverage
for(file in gatk_files) 
{
  
  
  infile1 <- paste(file,"_gene_haplotype_counts.txt",sep="")
  outfile <- paste(file,"_min_cov_20.txt",sep="")
  
  #infile1 <- "CRPA1_DNA_gene_haplotype_counts_no_zeros.txt"
  #outfile <- "CRPA1_DNA_gene_haplotype_counts_no_zeros.txt"
  
  haps <- read.table(infile1, header = TRUE, stringsAsFactors = FALSE)
  #head(haps)
  #nrow(haps)
  haps <- haps[which(haps$alleleOneCount >= 1 & haps$alleleTwoCount >= 1),]
  
  #write.table(haps, outfile, quote=FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  
}


#### ase stats results ####
#### binomial comparison ###

parental_species_jaw <- c("CPAJ1",
                          "CPAJ2",
                          "CPAJ3",
                          "CPMJ1",
                      "CPMJ2",
                      "CPMJ3",
                      "LLAJ1",
                      "LLAJ2",
                      "LLAJ3",
                      "LLMJ1",
                      "LLMJ2",
                      "LLMJ3")
parental_species_a_jaw <- c("LLAJ1",
                          "LLAJ2",
                          "LLAJ3",
                          "CPAJ1",
                          "CPAJ2",
                          "CPAJ3")
parental_species_m_jaw <- c("LLMJ1",
                            "LLMJ2",
                            "LLMJ3",
                            "CPMJ1",
                            "CPMJ2",
                            "CPMJ3")
hybrids <- c("axmJ1",
             "axmJ4",
             "axmJ5",
             "axmJ6")

sum_stats <- c()
#parental_ase_genes <- c()
for(file in hybrids) 
{
  
  
  infile1 <- paste(file,"_gene_haplotype_counts.txt",sep="")
  outfile <- paste(file,"_binomial_all.txt",sep="")
  
  #infile1 <- "axmJ1_min_cov_20.txt"
  #infile1 <- "CPAJ3_gene_haplotype_counts.txt"
  
  
  haps <- read.table(infile1, header = TRUE, stringsAsFactors = FALSE)
  head(haps)
  nrow(haps)
  haps$totalCount <- haps$alleleOneCount + haps$alleleTwoCount
  haps <- haps[which(haps$totalCount != 0),]
  
  binom_tests <- c()
  #log2_a_divided_by_ms <- c()
  
  for (i in c(1:nrow(haps)))
    
  {
    
    binom_test <- binom.test(haps$alleleOneCount[i], haps$totalCount[i], p = 0.5)$p.value
    binom_tests <- c(binom_tests, binom_test)
    #log2_a_divided_by_m <- log2(counts$LLA_allele_count[i] / counts$LLM_allele_count[i])
    #log2_a_divided_by_ms <- c(log2_a_divided_by_ms, log2_a_divided_by_m)
    
  }
  
  
  haps$binom_test_p_value <- binom_tests
  geneIndex <- haps$geneIndex
  counter <- 0
  ase_gene_names <- c()
  
  #for (j in geneIndex)
  #{
  #  
  #  geneIndex_table <- haps[which(haps$geneIndex == j),]
  #  p_vals <- geneIndex_table$binom_test_p_value
  #  
  #  if (all(p_vals <= 0.05))
  #  {
  #    counter <- counter + 1
  #    ase_gene_names <- c(ase_gene_names, j)
  #  }
  #  
  #}
  
  #ase_genes <- counter
  #total_genes <- length(unique(haps$geneIndex))
  
  #sum_stat <- c(file, (ase_genes/total_genes))
  #sum_stats <- c(sum_stats, sum_stat)
  #parental_ase_genes <- c(parental_ase_genes, ase_gene_names)
  
  #print(file)
  #print(ase_genes)
  #print(total_genes)
  #print(ase_genes/total_genes)
  
  #haps <- haps[haps$geneIndex %in% ase_gene_names, ]

  #write.table(haps, outfile, quote=FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  
}  
length(setdiff(hybrid_ase_genes, parental_ase_genes))
setdiff(hybrid_ase_genes, parental_ase_genes)

setwd("c:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/axm/ASE/skelly/min_cov_20/")
setwd("c:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/axm/ASE/skelly")

ase_skelly <- read.table("RNAresults_axmJ1_500000_20_200_hat10.gz_signASE.csv", header = FALSE, stringsAsFactors = FALSE) 
ase_skelly <- ase_skelly[order(ase_skelly$V2),]
head(ase_skelly)
nrow(ase_skelly)
length(intersect(ase_gene_names, ase_skelly$V1))
am$tag <- paste(am$symbol, am$Geneid, sep = ":")
length(intersect(ase_gene_names, am$tag))
length(intersect(ase_skelly$V1, am$tag))


# binomial sugggests ~15% of genes show ase for every snp across gene



# identify genes that are ase in hybrids only

setwd("c:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/axm/ASE/skelly")

axmJ1 <- read.csv("axmJ6_gene_haplotype_counts.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
axmJ4 <- read.csv("axmJ6_gene_haplotype_counts.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
axmJ5 <- read.csv("axmJ6_gene_haplotype_counts.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
axmJ6 <- read.csv("axmJ6_gene_haplotype_counts.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

axmJ1 <- read.csv("axmJ1_binomial.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
axmJ4 <- read.csv("axmJ4_binomial.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
axmJ5 <- read.csv("axmJ5_binomial.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
axmJ6 <- read.csv("axmJ6_binomial.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
axmJ1$snp_tag <- paste(axmJ1$geneIndex, axmJ1$snpIndex, sep = ":")
axmJ4$snp_tag <- paste(axmJ4$geneIndex, axmJ4$snpIndex, sep = ":")
axmJ5$snp_tag <- paste(axmJ5$geneIndex, axmJ5$snpIndex, sep = ":")
axmJ6$snp_tag <- paste(axmJ6$geneIndex, axmJ6$snpIndex, sep = ":")
axmJ1$totalCount <- axmJ1$alleleOneCount + axmJ1$alleleTwoCount
axmJ4$totalCount <- axmJ4$alleleOneCount + axmJ4$alleleTwoCount
axmJ5$totalCount <- axmJ5$alleleOneCount + axmJ5$alleleTwoCount
axmJ6$totalCount <- axmJ6$alleleOneCount + axmJ6$alleleTwoCount
axmJ1 <- axmJ1[which(axmJ1$totalCount >= 20 & axmJ1$alleleOneCount > 0 & axmJ1$alleleTwoCount > 0),]
axmJ4 <- axmJ4[which(axmJ4$totalCount >= 20 & axmJ4$alleleOneCount > 0 & axmJ4$alleleTwoCount > 0),]
axmJ5 <- axmJ5[which(axmJ5$totalCount >= 20 & axmJ5$alleleOneCount > 0 & axmJ5$alleleTwoCount > 0),]
axmJ6 <- axmJ6[which(axmJ6$totalCount >= 20 & axmJ6$alleleOneCount > 0 & axmJ6$alleleTwoCount > 0),]

range(axmJ6$alleleTwoCount)
head(axmJ6)
a <- merge(axmJ1, axmJ4, by = c("snp_tag"))
a <- merge(a, axmJ5, by = c("snp_tag"))
a <- merge(a, axmJ6, by = c("snp_tag"))
head(a)
nrow(a)
a <- merge(axmJ1, axmJ4, by = c("geneIndex"))
a <- merge(a, axmJ5, by = c("geneIndex"))
a <- merge(a, axmJ6, by = c("geneIndex"))
head(a)
nrow(a)
range(a$binom_test_p_value.x)

count_genes_with_het_snps <-  count(a, 'geneIndex.x')
head(count_genes_with_het_snps)
hist(count_genes_with_het_snps$freq, ylim = c(0,250))
nrow(count_genes_with_het_snps[which(count_genes_with_het_snps$freq == 1),])
nrow(count_genes_with_het_snps)

a <- cbind(a, colsplit(a[,1], ":", c("symbol", "Geneid")))

hybrid_ase_genes <- unique(a$geneIndex)
length(unique(a$geneIndex))
length(setdiff(hybrid_ase_genes, parental_ase_genes))



ham <- read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/axm/conditions/stringent_filtering/DE_ham_jaw_genes_sig.csv", header = TRUE, stringsAsFactors = FALSE)
ha <- read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/axm/conditions/stringent_filtering/DE_ha_jaw_genes_sig.csv", header = TRUE, stringsAsFactors = FALSE)
hm <- read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/axm/conditions/stringent_filtering/DE_hm_jaw_genes_sig.csv", header = TRUE, stringsAsFactors = FALSE)
amp <- read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/blast/all_DE_master.csv", header = TRUE, stringsAsFactors = FALSE)
head(ham)
am <- amp[which(amp$comp == "am_jaw"),]
head(axmJ1)
head(am)
nrow(am)

all_cis <- merge(a, ham, by = c("symbol"))
head(all_cis)
nrow(all_cis)
unique(all_cis$symbol)
range(amp$padj)

axmJ1_cis <- merge(axmJ1, am, by = c("Geneid"))
axmJ4_cis <- merge(axmJ4, am, by = c("Geneid"))
axmJ5_cis <- merge(axmJ5, am, by = c("Geneid"))
axmJ6_cis <- merge(axmJ6, am, by = c("Geneid"))

#### make plot like McManus 2010 4B

ham <- read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/axm/conditions/stringent_filtering/DE_ham_jaw_genes.csv", header = TRUE, stringsAsFactors = FALSE)
ha <- read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/axm/conditions/stringent_filtering/DE_ha_jaw_genes.csv", header = TRUE, stringsAsFactors = FALSE)
hm <- read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/axm/conditions/stringent_filtering/DE_hm_jaw_genes.csv", header = TRUE, stringsAsFactors = FALSE)
amp <- read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/mapping_correction/DESeq_genes_all.csv", header = TRUE, stringsAsFactors = FALSE)
head(amp)
am <- amp[which(amp$comp == "am_jaw"),]
am <- am[which(am$feature == "gene"),]
nrow(am)
head(am)
am$tag <- paste(am$symbol, am$Geneid, sep = ":")
length(unique(am$tag))

am$log2FoldChange_am <- am$log2FoldChange
am$padj_am <- am$padj
ham$log2FoldChange_ham <- (ham$log2FoldChange) * -1
ham$tag <- paste(ham$symbol, ham$Geneid, sep = ":")
ha$log2FoldChange_ha <- ha$log2FoldChange
ha$tag <- paste(ha$symbol, ha$Geneid, sep = ":")
hm$log2FoldChange_hm <- (hm$log2FoldChange) * -1
hm$tag <- paste(hm$symbol, hm$Geneid, sep = ":")
ham$padj_ham <- ham$padj
ha$padj_ha <- ha$padj
hm$padj_hm <- hm$padj
head(ham)
keeps <- c("tag", "log2FoldChange_ham", "log2FoldChange_ha", "padj_ham", "padj_ha")
all_genes <- merge(ham, ha, by = ("tag"))
all_genes <- all_genes[keeps]
all_genes <- merge(all_genes, hm, by = ("tag"))
keeps <- c("tag", "log2FoldChange_ham", "log2FoldChange_ha", "padj_ham", "padj_ha", "padj_hm", "log2FoldChange_hm")
all_genes <- all_genes[keeps]
all_genes <- merge(all_genes, am, by = ("tag"))
keeps <- c("tag", "log2FoldChange_ham", "log2FoldChange_ha", "padj_ham", "padj_ha", "padj_hm", "log2FoldChange_hm", "padj_am", "log2FoldChange_am", "Start", "End")
all_genes <- all_genes[keeps]
head(all_genes)
nrow(all_genes)
length(unique(all_genes$tag))

total_genes <- nrow(all_genes)
con <- all_genes[which(all_genes$padj_ham > 0.05 & all_genes$padj_ha > 0.05 & all_genes$padj_hm > 0.05 & all_genes$padj_am > 0.05),]
con$inheritance <- "conserved"
con$inheritance <- 1
nrow(con)
add <- all_genes[which(all_genes$padj_ham > 0.05 & all_genes$padj_am < 0.05),]
add1 <- add[which(add$log2FoldChange_ha < 0 & add$log2FoldChange_hm > 0),]
add2 <- add[which(add$log2FoldChange_ha > 0 & add$log2FoldChange_hm < 0),]
add <- rbind(add1, add2)
add$inheritance <- "additive"
add$inheritance <- 2
nrow(add)
a_dom <- all_genes[which(all_genes$padj_ha > 0.05 & all_genes$padj_am < 0.05 & all_genes$padj_hm < 0.05),]
a_dom$inheritance <- "a_dominant"
a_dom$inheritance <- 3
nrow(a_dom)
m_dom <- all_genes[which(all_genes$padj_ha < 0.05 & all_genes$padj_am < 0.05 & all_genes$padj_hm > 0.05),]
m_dom$inheritance <- "m_dominant"
m_dom$inheritance <- 4
nrow(m_dom)
over_dom <- all_genes[which(all_genes$padj_ham < 0.05 & all_genes$log2FoldChange_ha > 0 & all_genes$log2FoldChange_hm > 0),]
over_dom$inheritance <- "over_dominant"
over_dom$inheritance <- 5
nrow(over_dom)
under_dom <- all_genes[which(all_genes$padj_ham < 0.05 & all_genes$log2FoldChange_ha < 0 & all_genes$log2FoldChange_hm < 0),]
under_dom$inheritance <- "under_dominant"
under_dom$inheritance <- 6
nrow(under_dom)

inheritance <- rbind(con, add, a_dom, m_dom, over_dom, under_dom)
head(inheritance)
nrow(inheritance)
length(unique(inheritance$tag))
inheritance_and_ase <- inheritance[inheritance$tag %in% hybrid_ase_genes,]
length(unique(inheritance_and_ase$tag))

head(inheritance_and_ase)

inheritance$gene_length <- inheritance$End - inheritance$Start
gene_lengths <- merge(inheritance, inheritance_and_ase, suffixes = c("",".y"), by = ("tag"))
nrow(gene_lengths)
head(gene_lengths)
inheritance_and_ase <- gene_lengths
unique(inheritance_and_ase$ase)
trans_gene_lengths <- ase_type[which(ase_type$ase_type == 4),]
am$tag <- paste(am$symbol,am$Geneid, sep = ":" )
am$gene_length <- am$End - am$Start
am <- am[which(am$padj <= 0.05),]
head(am)
other_gene_lengths <- am[! am$tag %in% trans_gene_lengths$tag,] 
boxplot(unique(trans_gene_lengths$gene_length),unique(other_gene_lengths$gene_length), names=c("DE am hets", "DE am no hets"), 
        col="grey")
t.test(unique(trans_gene_lengths$gene_length),unique(other_gene_lengths$gene_length))
hets_test <- read.table("c:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/axm/ASE/fisher_test_hets.txt", header = FALSE, stringsAsFactors = FALSE)
fisher.test(hets_test)

#write.table(inheritance_and_ase, "c:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/axm/ASE/inheritance_and_ase.txt", quote=FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
inheritance <- read.table("c:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/axm/ASE/inheritance.txt", header = TRUE, stringsAsFactors = FALSE)
inheritance_and_ase <- read.table("c:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/axm/ASE/inheritance_and_ase.txt", header = TRUE, stringsAsFactors = FALSE)


cols <- 	c("#2E5894", "#C46210", "#D98695",
           "#000000", "#9C2542", "#319177")

cols <- 	c("#0067a7", "#e97600", "#964f8e",
           "#000000", "#bd1e24", "#f6c700")
cols <- 	c("#f6c700", "#76D7EA", "#FF00CC",
           "#000000", "#bd1e24", "#0067a7")

cols <- 	c("#f6c700", "#E77200", "#FF00CC",
           "#000000", "#bd1e24", "#0067a7")

#(11120)	w3-safety-red	#bd1e24
#(12300)	w3-safety-orange	#e97600
#(13591)	w3-safety-yellow	#f6c700
#(14120)	w3-safety-green	#007256
#(15092)	w3-safety-blue	#0067a7
#(17155)	w3-safety-purple	#964f8e

cols_in <- cols[inheritance$inheritance]

#tiff("c:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/axm/ASE/skelly/plots/inheritance_patterns.tif", width = 7, height = 7, units = 'in', res = 1000)
plot(inheritance$log2FoldChange_ha, inheritance$log2FoldChange_hm, pch =16, col = cols_in, 
     ylab = "log2 fold change Hybrids vs. Snail-eaters",
     xlab = "log2 fold change Hybrids vs. Generalists",cex.axis=1.5)
abline(v=0, col = "black", lty = 3, lwd = 1.8)
abline(h=0, col = "black", lty = 3, lwd = 1.8)
#dev.off()

cols_bar <- 	c("#f6c700", "#0067a7", "#bd1e24",
               "#000000", "#FF00CC", "#E77200")
inheritance_cats <- c(nrow(con), nrow(under_dom), nrow(over_dom), nrow(m_dom), nrow(a_dom), nrow(add))
#tiff("c:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/axm/ASE/skelly/plots/inheritance_barplot.tif", width = 7, height = 5, units = 'in', res = 1000)
barplot(inheritance_cats,horiz=TRUE, xlim = c(0,7000), col = cols_bar)
#dev.off()


trans <- inheritance_and_ase[which(inheritance_and_ase$padj_am <= 0.05 & inheritance_and_ase$ase_all_h == 1),]
nrow(trans)
length(unique(trans$tag))
trans$ase_type <- 4
# no compensatory
# comp <- inheritance_and_ase[which(inheritance_and_ase$padj_am <= 0.05 & inheritance_and_ase$ase != 4),]
#number of DE axm genes included in f1 ase analyses
#de_am <- am[which(am$padj <= 0.05),]
#de_am$tag <- paste(de_am$symbol, de_am$Geneid, sep = ":")
all <- c(axmJ1$geneIndex, axmJ4$geneIndex,axmJ5$geneIndex,axmJ6$geneIndex)
length(intersect(unique(a$geneIndex), de_am$tag))
#nrow(inheritance_and_ase[inheritance_and_ase$tag %in% de_am$tag,])
mis_over <- inheritance_and_ase[which(inheritance_and_ase$padj_ham <= 0.05 & inheritance_and_ase$log2FoldChange_ham > 0 & inheritance_and_ase$ase_all_h == 1),]
mis_over$ase_type <- 2
mis_under <- inheritance_and_ase[which(inheritance_and_ase$padj_ham <= 0.05 & inheritance_and_ase$log2FoldChange_ham < 0 & inheritance_and_ase$ase_all_h == 1),]
mis_under$ase_type <- 3
mis_ase_over <- inheritance_and_ase[which(inheritance_and_ase$padj_ham <= 0.05 & inheritance_and_ase$log2FoldChange_ham > 0 & inheritance_and_ase$ase_all_h == 2),]
nrow(mis_ase_over)
mis_ase_over$ase_type <- 5
mis_ase_under <- inheritance_and_ase[which(inheritance_and_ase$padj_ham <= 0.05 & inheritance_and_ase$log2FoldChange_ham < 0 & inheritance_and_ase$ase_all_h == 2),]
nrow(mis_ase_under)
mis_ase_under$ase_type <- 6
con <- inheritance_and_ase[which(inheritance_and_ase$padj_ham >= 0.05 & inheritance_and_ase$padj_am >= 0.05 & inheritance_and_ase$ase_all_h == 1),]
con$ase_type <- 1
#cis <- inheritance_and_ase[which(inheritance_and_ase$padj_am <= 0.05 & inheritance_and_ase$ase_all_h == 2),]
#cis$ase_type <- 5
# no cis
nrow(inheritance_and_ase[which(inheritance_and_ase$padj_am <= 0.05),])

ase_type <- rbind(con, mis_over, mis_under, trans, mis_ase_over, mis_ase_under)
head(ase_type)
nrow(ase_type)

#cols <- 	c("#f6c700", "#E77200", "#FF00CC","#bd1e24", "#0067a7")
cols <- 	c("#f6c700", "#bd1e24", "#0067a7","#E77200", "#000000", "#FF00CC")
cols_ase <- cols[ase_type$ase_type]

#tiff("c:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/axm/ASE/skelly/plots/cis_trans_patterns.tif", width = 7, height = 7, units = 'in', res = 1000)
plot(ase_type$log2FoldChange_am, ase_type$log2FoldChange_hm, col = cols_ase, 
     pch = c(16, 16, 16, 17,17,17)[ase_type$ase_type],
     cex = c(1, 1, 1, 1.6,1.6,1.6)[ase_type$ase_type],
     ylab = "log2 fold change Hybrids vs. Parental Species",
     xlab = "log2 fold change Generalists vs. Snail-eaters",cex.axis=1.5)
abline(v=0, col = "black", lty = 3, lwd = 1.8)
abline(h=0, col = "black", lty = 3, lwd = 1.8)
#dev.off()


# PCA #
pca_data <- read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/axm/axm_pca.csv", header=T, stringsAsFactors = F)
head(pca_data,100)

cols <- c("#bd1e24", "#0067a7", "#964f8e")
cols_pca <- cols[pca_data$species]
#setwd("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/axm/")
tiff("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/Manuscripts/misexpression/ASE_plots/pca_axm_jaw.tif", width = 6, height = 5.5, units = 'in', res = 600)
plot(pca_data$PC1, pca_data$PC2, col= cols_pca, 
     pch = c(16, 17)[as.numeric(pca_data$lake_num)], ylab = 'PC 2', xlab = 'PC 1',cex = 2.1, cex.axis = 1.3)
dev.off()

head(inheritance)
down <- inheritance[which(inheritance$log2FoldChange_ham < 0 & inheritance$padj_ham <= 0.05),]
down$type <- "underdominant"
#down$type <- 2
down$log2FoldChange_ham <- abs(down$log2FoldChange_ham)
up <- inheritance[which(inheritance$log2FoldChange_ham > 0 & inheritance$padj_ham <= 0.05),]
up$type <- "overdominant"
#up$type <- 1
mags <- rbind(up, down)
head(mags)
nrow(mags)
boxplot(mags$log2FoldChange_ham~mags$type, main = "jaws")

library(ggplot2)
# Basic violin plot

#tiff("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/axm/over_under_magnitude.tif", width = 5.5, height = 6, units = 'in', res = 600)
p <- ggplot(mags, aes(x=type, y=log2FoldChange_ham, fill = type)) + 
  geom_violin(trim = TRUE)
p+ scale_fill_manual(values=c("#bd1e24", "#0067a7")) + theme(axis.text.y = element_text(size=18, color = '#000000'),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),axis.line = element_line(colour = "black"),
          axis.text.x = element_text(size=18, color = '#000000'), legend.position="none",
          axis.title.y=element_blank(),axis.title.x=element_blank()) + stat_summary(fun.data=data_summary, col = "black")
#dev.off()
wilcox.test(mags$log2FoldChange_ham~mags$type)
install.packages('pwr')
library(pwr)
wilcoxsign_test(log2FoldChange_ham~type, data = mags, distribution="exact")
pwr.f2.test(u = 2, f2 = 0.3/(1 - 0.3), sig.level = 0.05, power = 0.8)

# variance and bartlett test
axm_var <- read.table("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/axm/axm_variance_table.txt", header=T, stringsAsFactors = F)
axm_var$baseVar_axm <- axm_var$baseVar
axm_var$baseMean_axm <- axm_var$baseMean
axm_var$species <- 2
a_var   <- read.table("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/axm/a_variance_table.txt", header=T, stringsAsFactors = F)
a_var$baseVar_a <- a_var$baseVar
a_var$baseMean_a <- a_var$baseMean
a_var$species <- 1
m_var   <- read.table("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/axm/m_variance_table.txt", header=T, stringsAsFactors = F)
m_var$baseVar_m <- m_var$baseVar
m_var$baseMean_m <- m_var$baseMean
m_var$species <- 3
head(m_var)
head(axm_var)

am <- merge(a_var, m_var, by = c("Geneid"))
all_var <- merge(am, axm_var, by = c("Geneid"))
keeps <- c("Geneid", "species", "baseVar_axm", "baseVar_a", "baseVar_m", "baseMean_axm", "baseMean_a", "baseMean_m")
all_var <- all_var[keeps]
head(all_var)

#tiff("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/axm/over_under_magnitude.tif", width = 5.5, height = 6, units = 'in', res = 600)
library(vioplot)
#install.packages("vioplot")
mean(all_var$baseVar_axm)
mean(all_var$baseVar_a)
mean(all_var$baseVar_m)
median(all_var$baseVar_axm)
median(all_var$baseVar_a)
median(all_var$baseVar_m)
range(all_var$baseVar_axm)
range(all_var$baseVar_a)
range(all_var$baseVar_m)

bartlett.test(list(all_var$baseVar_axm, all_var$baseVar_a, all_var$baseVar_m))
bartlett.test(list(all_var$baseVar_axm, all_var$baseVar_a))
bartlett.test(list(all_var$baseVar_axm, all_var$baseVar_m))
bartlett.test(list(all_var$baseVar_a, all_var$baseVar_m))

bartlett.test(list(all_var$baseMean_axm, all_var$baseMean_a, all_var$baseVar_m))
bartlett.test(list(all_var$baseMean_axm, all_var$baseMean_a))
bartlett.test(list(all_var$baseMean_axm, all_var$baseMean_m))
bartlett.test(list(all_var$baseMean_a, all_var$baseMean_m))


boxplot(all_var$baseMean_axm,all_var$baseMean_a,all_var$baseMean_m, names=c("4 cyl", "6 cyl", "8 cyl"), ylim = c(0,200), 
        col="gold")

#dev.off()


hist(inheritance$log2FoldChange_am, breaks = 50)
nrow(inheritance[which(inheritance$log2FoldChange_am > 0),]) / nrow(inheritance)
nrow(inheritance[which(inheritance$log2FoldChange_am < 0),]) / nrow(inheritance)

####### pleiotropy #####
library(Rcpp)
zfin <- read.table("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/pleiotropy/gene_association.jam.zfin", fill=TRUE,header = FALSE,quote = "", stringsAsFactors = FALSE, sep = "\t")
head(zfin)
zfin <- zfin[which(zfin$V9 == "P"),]
zfin <- zfin[which(zfin$V7 == "EXP" | zfin$V7 == "IDA"| zfin$V7 == "IPI" | zfin$V7 == "IMP" | zfin$V7 == "IGI" | zfin$V7 == "IEP"),]
head(zfin)

over <- read.table("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/axm/ASE/over_cranial_genes.txt",header = FALSE,quote = "", stringsAsFactors = FALSE, sep = "\t")
under <- read.table("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/axm/ASE/under_cranial_genes.txt",header = FALSE,quote = "", stringsAsFactors = FALSE, sep = "\t")
con <- read.table("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/axm/ASE/conserved_cranial_genes.txt",header = FALSE,quote = "", stringsAsFactors = FALSE, sep = "\t")
over <- read.table("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/axm/ASE/over_genes.txt",header = FALSE,quote = "", stringsAsFactors = FALSE, sep = "\t")
under <- read.table("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/axm/ASE/under_genes.txt",header = FALSE,quote = "", stringsAsFactors = FALSE, sep = "\t")
con <- read.table("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/axm/ASE/conserved_genes.txt",header = FALSE,quote = "", stringsAsFactors = FALSE, sep = "\t")

head(over)
nrow(under)
over_zeb <- zfin[zfin$V3 %in% over$V1,]
over_zeb$go_index <- paste(over_zeb$V3, over_zeb$V5, sep = ";")
over_zeb <- subset(over_zeb, !duplicated(over_zeb[,18]))
over_zeb_count <- count(over_zeb, "V3")
over_zeb_count$type <- "over"
head(over_zeb_count)
nrow(over_zeb_count)
under_zeb <- zfin[zfin$V3 %in% under$V1,]
under_zeb$go_index <- paste(under_zeb$V3, under_zeb$V5, sep = ";")
under_zeb <- subset(under_zeb, !duplicated(under_zeb[,18]))
under_zeb_count <- count(under_zeb, "V3")
under_zeb_count$type <- "under"
head(under_zeb_count)
nrow(under_zeb_count)
con_zeb <- zfin[zfin$V3 %in% con$V1,]
con_zeb$go_index <- paste(con_zeb$V3, con_zeb$V5, sep = ";")
con_zeb <- subset(con_zeb, !duplicated(con_zeb[,18]))
con_zeb_count <- count(con_zeb, "V3")
con_zeb_count$type <- "con"
con_zeb_count$type2 <- "con"
under_zeb_count$type2 <- "mis"
over_zeb_count$type2 <- "mis"
head(con_zeb_count)
nrow(con_zeb_count)
all_cranial_go <- rbind(con_zeb_count, over_zeb_count, under_zeb_count)
#library(ggplot2)
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  #xzy <- ymin * -1
  #ymin <- ymin + xzy +1
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

#data_summary <- function(x) {
#  m <- mean(x)
#  sem<-sd(x)/sqrt(length(x))
#  ymin <- m-sem
#  ymax <- m+sem
#  return(c(y=m,ymin=ymin,ymax=ymax))
#}

median(con_zeb_count$freq)
median(over_zeb_count$freq)
median(under_zeb_count$freq)
mean(con_zeb_count$freq)
mean(over_zeb_count$freq)
mean(under_zeb_count$freq)
library(ggplot2)
#tiff("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/axm/pleio_all_go.tif", width = 5.5, height = 6, units = 'in', res = 600)
p <- ggplot(all_cranial_go, aes(x=type, y=freq, fill = type)) +  #ylim(0,20) +
  #geom_boxplot()
geom_violin(trim = TRUE, bw = 2)#bw = 4
p+ scale_fill_manual(values=c("#f6c700", "#bd1e24", "#0067a7")) + theme(axis.text.y = element_text(size=18, color = '#000000'),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                             panel.background = element_blank(),axis.line = element_line(colour = "black"),
                                                             axis.text.x = element_text(size=18, color = '#000000'), legend.position="none",
                                                             axis.title.y=element_blank(),axis.title.x=element_blank()) + stat_summary(fun.data=data_summary, col = "black")
#dev.off()

ppi_lm <- glm(all_cranial_go$freq~all_cranial_go$type2, family = poisson)
ppi_lm <- glm.nb(all_cranial_go$freq~all_cranial_go$type)
summary(ppi_lm)
glm_res <- resid(ppi_lm)
plot(all_cranial_go$freq, glm_res)

ppi <- read.table("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/pleiotropy/danio_ppi_experimental", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
head(ppi)
range(ppi$experimental)
range(ppi$combined_score)
quantile(ppi$experimental,.95)
quantile(ppi$combined_score,.95)
ppi_sig <- ppi[which(ppi$combined_score >= 481 & ppi$experimental >= 265),]
#ppi_sig <- ppi
ppi_counts <-  count(ppi_sig, 'protein1')
head(ppi_counts)

#ensemble <- read.table("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/axm/ASE/ppi_ensembl_ids.txt",header = TRUE,quote = "", stringsAsFactors = FALSE, sep = "\t")
ensemble <- read.table("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/axm/ASE/all_con_mis_ensemble_ids.txt",header = TRUE,quote = "", stringsAsFactors = FALSE, sep = "\t")

head(ensemble)
ensemble$tag <- paste(ensemble$Gene, ensemble$Symbol, sep= ":")
ensemble$protein1 <- ensemble$Symbol
ppi_genes <- merge(ensemble, ppi_counts, by = c("protein1"))
head(ppi_genes)
nrow(ppi_genes)
length(unique(ppi_genes$protein1))
over$Gene <- over$V1
under$Gene <- under$V1
con$Gene <- con$V1
over_ppi <- merge(over, ppi_genes, by = ("Gene"))
over_ppi$type <- "over"
under_ppi <- merge(under, ppi_genes, by = ("Gene"))
under_ppi$type <- "under"
con_ppi <- merge(con, ppi_genes, by = ("Gene"))
con_ppi$type <- "con"
con_ppi$type2 <- "con"
under_ppi$type2 <- "mis"
over_ppi$type2 <- "mis"
nrow(con_ppi)
nrow(over_ppi)
nrow(under_ppi)
#under_ppi <- under_ppi[which(under_ppi$freq < 1000),]
all_ppi <- rbind(con_ppi, under_ppi, over_ppi)


#tiff("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/axm/pleio_all_ppi.tif", width = 5.5, height = 6, units = 'in', res = 600)
p <- ggplot(all_ppi, aes(x=type, y=freq, fill = type)) + 
  #geom_boxplot()
  geom_violin(trim = TRUE, bw = 10)#+geom_dotplot(binaxis='y', stackdir='center',
                                   #position=position_dodge(1),dotsize = .3)
p+ scale_fill_manual(values=c("#f6c700", "#bd1e24", "#0067a7")) + theme(axis.text.y = element_text(size=18, color = '#000000'),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                       panel.background = element_blank(),axis.line = element_line(colour = "black"),
                                                                       axis.text.x = element_text(size=18, color = '#000000'), legend.position="none",
                                                                       axis.title.y=element_blank(),axis.title.x=element_blank()) + stat_summary(fun.data=data_summary, col = "black")
#dev.off()
median(con_ppi$freq)
median(under_ppi$freq)
median(over_ppi$freq)
mean(con_ppi$freq)
mean(under_ppi$freq)
mean(over_ppi$freq)

ppi_lm <- glm(all_ppi$freq~all_ppi$type2, family = poisson)
ppi_lm <- glm.nb(all_ppi$freq~all_ppi$type)
summary(ppi_lm)

# tissue expression
tis <- read.table("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/pleiotropy/Danio_rerio_expr_simple_development.tsv", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
head(tis)
nrow(tis)
tis2 <- cbind(tis, colsplit(tis$Developmental.stage.ID, ":", c("stage", "junk")))
head(tis2)
unique(tis2$stage)
tis2 <- tis2[which(tis2$Call.quality == "gold quality" &  tis2$stage == "ZFS"),]
nrow(tis2)
unique(tis2$Anatomical.entity.name)
tis2$tag <- paste(tis2$Gene.name, tis2$Developmental.stage.ID, sep = ":")
length(tis2)
head(tis2)
tis2 <- subset(tis2    , !duplicated(tis2[,12]))
unique(tis2$Developmental.stage.name)
tis2 <- tis2[grep("Day", tis2$Developmental.stage.name), ]
#unique(tis2$Developmental.stage.name)
tis2$Gene <- tis2$Gene.name
tis <- tis2
head(tis)
tis_counts <- count(tis, 'Gene')
head(tis_counts)
range(tis_counts$freq)

over_tis <- merge(over, tis_counts, by = ("Gene"))
over_tis$type <- "over"
under_tis <- merge(under, tis_counts, by = ("Gene"))
under_tis$type <- "under"
con_tis <- merge(con, tis_counts, by = ("Gene"))
con_tis$type <- "con"
con_tis$type2 <- "con"
under_tis$type2 <- "mis"
over_tis$type2 <- "mis"
all_tis <- rbind(con_tis, under_tis, over_tis)

#tiff("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/axm/over_under_magnitude.tif", width = 5.5, height = 6, units = 'in', res = 600)
p <- ggplot(all_tis, aes(x=type, y=freq, fill = type)) + 
  geom_boxplot()
p+ scale_fill_manual(values=c("#f6c700", "#bd1e24", "#0067a7")) + theme(axis.text.y = element_text(size=18, color = '#000000'),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                       panel.background = element_blank(),axis.line = element_line(colour = "black"),
                                                                       axis.text.x = element_text(size=18, color = '#000000'), legend.position="none",
                                                                       axis.title.y=element_blank(),axis.title.x=element_blank()) #+ stat_summary(fun.data=data_summary, col = "black")
#dev.off()



ppi_lm <- glm(all_tis$freq~all_tis$type2, family = poisson)
ppi_lm <- glm.nb(all_tis$freq~all_tis$type2)
summary(ppi_lm)








###### Drafts #####

####  CLUSTER Skelly  ####

# RNA #

gatk_files <- c("axmE1",
                "axmE2",
                "axmJ1",
                "axmJ4",
                "axmJ5",
                "axmJ6",
                "CPAE1",
                "CPAE2",
                "CPAE3",
                "CPAJ1",
                "CPAJ2",
                "CPAJ3",
                "CPME1",
                "CPME2",
                "CPME3",
                "CPMJ1",
                "CPMJ2",
                "CPMJ3",
                "CPPE1",
                "CPPE2",
                "CPPE3",
                "CPPJ1",
                "CPPJ2",
                "LLAE1",
                "LLAE2",
                "LLAE3",
                "LLAJ1",
                "LLAJ2",
                "LLAJ3",
                "LLME1",
                "LLME2",
                "LLME3",
                "LLMJ1",
                "LLMJ2",
                "LLMJ3")
#"CRPA1_DNA",
#"CRPM1_DNA",
#"LILA1_DNA",
#"LILM3_DNA")
for(file in gatk_files) 
{
  
  
  infile1 <- paste(file,"_allele_counts_phased_gatk.table",sep="")
  infile2 <- paste(file,"_snp_table",sep="")
  outfile <- paste(file,"_allele_counts.csv",sep="")
  
  #setwd("c:/Users/jmcgirr/Desktop/")
  #counts <- read.table("test_allele_counts_phased_gatk.table", header = TRUE)
  #dat <- read.table("test_snp_table.txt", header = TRUE)
  #counts <- read.table("CPAE1_allele_counts_phased_gatk.table", header = TRUE)
  #dat <- read.table("CPAE1_snp_table.txt", header = TRUE)
  
  
  dat <- read.table(infile, header = TRUE, stringsAsFactors = FALSE)
  dat <- na.omit(dat)
  dat <- cbind(dat, colsplit(dat[,4], ",", c("HP1", "HP2")))
  dat <- cbind(dat, colsplit(dat[,5], "-", c("junk1", "HP1.1")))
  dat <- cbind(dat, colsplit(dat[,6], "-", c("junk2", "HP2.1")))
  dat <- cbind(dat, colsplit(dat[,3], "/", c("GT_1", "GT_2")))
  dat1 <- dat[which(dat$HP1.1 == 1 ),]
  dat2 <- dat[which(dat$HP1.1 == 2 ),]
  dat1["sample1.SNPs.haplotype_A"] <- dat1$GT_1
  dat2["sample1.SNPs.haplotype_A"] <- dat2$GT_2
  dat1["sample1.SNPs.haplotype_B"] <- dat1$GT_2
  dat2["sample1.SNPs.haplotype_B"] <- dat2$GT_1
  final <- rbind(dat1, dat2)
  final$contig <- final$CHROM
  final$position <- final$POS
  keeps <- c("contig", "position", "sample1.SNPs.haplotype_A", "sample1.SNPs.haplotype_B")
  final_snps <- final[keeps]
  head(final_snps)
  
  counts <- read.table(infile2, header = TRUE, stringsAsFactors = FALSE)
  keeps <- c("contig", "position", "refAllele", "altAllele", "refCount", "altCount")
  final_counts <- counts[keeps]
  head(final_counts)
  
  merge <- merge(final_counts, final_snps, by = c("contig","position"))
  head(merge)
  ref_allele <- merge$refAllele
  alt_allele <- merge$altAllele
  ref_count <- merge$refCount
  alt_count <- merge$altCount
  hap_a <- merge$sample1.SNPs.haplotype_A
  hap_b <- merge$sample1.SNPs.haplotype_B
  hap_a_counts <- c()
  hap_b_counts <- c()
  for (i in c(1:nrow(merge)))
  {
    
    if (ref_allele[i] == hap_a[i])
    {
      hap_a_count <- ref_count[i]
      hap_b_count <- alt_count[i]
    }
    else if (ref_allele[i] != hap_a[i])
    {
      hap_a_count <- alt_count[i]
      hap_b_count <- ref_count[i]
    }
    
    hap_a_counts <- c(hap_a_counts, hap_a_count)
    hap_b_counts <- c(hap_b_counts, hap_b_count) 
  }
  merge$hap_a_counts <- hap_a_counts
  merge$hap_b_counts <- hap_b_counts
  
  keeps <- c("contig", "position","hap_a_counts","hap_b_counts")
  final <- merge[keeps]
  
  #write.table(final,quote=FALSE, sep = "\t", row.names = FALSE, col.names = TRUE, outfile)
}




install.packages("runjags")
library(coda)


dnaD <- na.omit(read.table(file, header = T, na.strings = 'NA'))
dnaD$total <- apply(dnaD[, c(3,4)], 1, sum)



source('readGzippedMcmcOutput.R')
result <- read.mcmc('test_unbiased.txt.gz')
n.iter <- result$n.iter/result$thin
burnin <- 0.1*n.iter
a.hat <- median(exp(result$mcmc$logA[burnin:n.iter]))
d.hat <- median(exp(result$mcmc$logD[burnin:n.iter]))














setwd("c:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/axm/ASE/skelly/")
gatk_files <- c("axmE1",
                "axmE2",
                "axmJ1",
                "axmJ4",
                "axmJ5",
                "axmJ6",
                "CPAE1",
                "CPAE2",
                "CPAE3",
                "CPAJ1",
                "CPAJ2",
                "CPAJ3",
                "CPME1",
                "CPME2",
                "CPME3",
                "CPMJ1",
                "CPMJ2",
                "CPMJ3",
                "CPPE1",
                "CPPE2",
                "CPPE3",
                "CPPJ1",
                "CPPJ2",
                "LLAE1",
                "LLAE2",
                "LLAE3",
                "LLAJ1",
                "LLAJ2",
                "LLAJ3",
                "LLME1",
                "LLME2",
                "LLME3",
                "LLMJ1",
                "LLMJ2",
                "LLMJ3")
#file <- "LLMJ3_allele_counts_gatk"
#infile <- "test_snp_table.txt"
#outfile <- paste(file,"_haplotypes.file",sep="")

for(file in gatk_files) 
{
  
  infile_haps <- paste(file,"_snp_table_haplotypes.file",sep="")
  infile_counts <- paste(file,"_allele_counts_gatk.table",sep="")
  outfile <- paste(file,"_allele_counts_phased_gatk.table",sep="")
  
  infile_haps <- "axmE1_snp_table_haplotypes.file"
  infile_counts <- "axmE1_allele_counts_gatk.table" 
  
  haps <- read.table(infile_haps, header = TRUE, stringsAsFactors = FALSE)
  counts <- read.table(infile_counts, header = TRUE, stringsAsFactors = FALSE)
  haps["contig"] <- haps$CHROM
  haps["position"] <- haps$POS
  final <- merge(haps, counts, by = c("contig", "position"))
  keeps <- c("contig", "position", "variantID", "refAllele", "altAllele", "refCount", "altCount", "totalCount", "lowMAPQDepth", "lowBaseQDepth", "rawDepth", "otherBases", "improperPairs")
  final <- final[keeps]
  write.table(final,quote=FALSE,sep="\t", row.names = FALSE, col.names = TRUE, outfile)
}





dat <- read.table("axmE1_allele_counts_gatk.table", header = TRUE, stringsAsFactors = FALSE)
dat <- na.omit(dat)
dat <- cbind(dat, colsplit(dat[,4], ",", c("HP1", "HP2")))
dat <- cbind(dat, colsplit(dat[,5], "-", c("junk1", "HP1.1")))
dat <- cbind(dat, colsplit(dat[,6], "-", c("junk2", "HP2.1")))
dat <- cbind(dat, colsplit(dat[,3], "/", c("GT_1", "GT_2")))
dat1 <- dat[which(dat$HP1.1 == 1 ),]
dat2 <- dat[which(dat$HP1.1 == 2 ),]
dat1["sample1.SNPs.haplotype_A"] <- dat1$GT_1
dat2["sample1.SNPs.haplotype_A"] <- dat2$GT_2
dat1["sample1.SNPs.haplotype_B"] <- dat1$GT_2
dat2["sample1.SNPs.haplotype_B"] <- dat2$GT_1
final <- rbind(dat1, dat2)
keeps <- c("CHROM", "POS", "sample1.SNPs.haplotype_A", "sample1.SNPs.haplotype_B")
final <- final[keeps]

final <- final[order(final$POS),]
dat <- dat[order(dat$POS),]
head(dat)
head(final)

#write.table(final, "c:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/axm/ASE/skelly/axmE1_haplotypes.file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)















#### MAMBA ####
###############
#### convert gatk ASEReadCounter datatables into mamba format ####

#gatk_files <- c("axmE1",
#                "axmE2",
#                "axmJ1",
#                "axmJ4",
#                "axmJ5",
#                "axmJ6",
#                "CPAE1",
#                "CPAE2",
#                "CPAE3",
#                "CPAJ1",
#                "CPAJ2",
#                "CPAJ3",
#                "CPME1",
#                "CPME2",
#                "CPME3",
#                "CPMJ1",
#                "CPMJ2",
#                "CPMJ3",
#                "CPPE1",
#                "CPPE2",
#                "CPPE3",
#                "CPPJ1",
#                "CPPJ2",
#                "LLAE1",
#                "LLAE2",
#                "LLAE3",
#                "LLAJ1",
#                "LLAJ2",
#                "LLAJ3",
#                "LLME1",
#                "LLME2",
#                "LLME3",
#                "LLMJ1",
#                "LLMJ2",
#                "LLMJ3")
#for(file in gatk_files) 
#{


infile <- paste(file,"_allele_counts_gatk.table",sep="")

dat <- read.table("mamba_r_input_test.txt", header = TRUE, stringsAsFactors = FALSE)

#choose a test gene region
#dat <-  dat[which(dat$contig == "NW_015150787.1" & dat$position > 100000 & dat$position < 109000),]
dat <-  dat[which(dat$refCount > 0 & dat$altCount > 0),]
#head(dat)
#nrow(dat)
#tail(dat)

dat["RSID"] <- paste(dat$contig, dat$position, sep = ":")
dat["INDIVIDUAL"] <- "ind1"
dat["REF_COUNT"] <- dat$refCount
dat["NONREF_COUNT"] <- dat$altCount
dat["TOTAL_COUNT"] <- dat$totalCount
dat["EXON_INFO"] <- "SEVERE_IMPACT=STOP_GAINED"
dat["TISSUE"] <- "TIS1"
keeps <- c("INDIVIDUAL", "RSID", "REF_COUNT", "NONREF_COUNT", "TOTAL_COUNT", "EXON_INFO", "TISSUE")

#create dataframe where tis1 and tis2 are identical for read counts. Only care about ASE in individual and comparing to ASE in non-hybrid
final <- dat[keeps]
#final <- head(final,((nrow(final))/3))
#final <- head(final,1000)
#head(final)
final2 <- final
final2["TISSUE"] <- "TIS2"
#head(final2)
final3 <- rbind(final,final2)
#tail(final3)
#head(final3)

write.table(final3, "mamba_input_test.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)






gatk_files <- c("axmE1",
                "axmE2",
                "axmJ1",
                "axmJ4",
                "axmJ5",
                "axmJ6",
                "CPAE1",
                "CPAE2",
                "CPAE3",
                "CPAJ1",
                "CPAJ2",
                "CPAJ3",
                "CPME1",
                "CPME2",
                "CPME3",
                "CPMJ1",
                "CPMJ2",
                "CPMJ3",
                "CPPE1",
                "CPPE2",
                "CPPE3",
                "CPPJ1",
                "CPPJ2",
                "LLAE1",
                "LLAE2",
                "LLAE3",
                "LLAJ1",
                "LLAJ2",
                "LLAJ3",
                "LLME1",
                "LLME2",
                "LLME3",
                "LLMJ1",
                "LLMJ2",
                "LLMJ3")
for(file in gatk_files) 
{
  
  
  infile <- paste(file,"_allele_counts_gatk.table",sep="")
  
  dat <- read.table(infile, header = TRUE, stringsAsFactors = FALSE)
  
  #choose a test gene region
  #dat <-  dat[which(dat$contig == "NW_015150787.1" & dat$position > 100000 & dat$position < 109000),]
  dat <-  dat[which(dat$refCount > 0 & dat$altCount > 0),]
  #head(dat)
  #nrow(dat)
  #tail(dat)
  
  dat["RSID"] <- paste(dat$contig, dat$position, sep = ":")
  dat["INDIVIDUAL"] <- "ind1"
  dat["REF_COUNT"] <- dat$refCount
  dat["NONREF_COUNT"] <- dat$altCount
  dat["TOTAL_COUNT"] <- dat$totalCount
  dat["EXON_INFO"] <- "SEVERE_IMPACT=STOP_GAINED"
  dat["TISSUE"] <- "TIS1"
  keeps <- c("INDIVIDUAL", "RSID", "REF_COUNT", "NONREF_COUNT", "TOTAL_COUNT", "EXON_INFO", "TISSUE")
  
  #create dataframe where tis1 and tis2 are identical for read counts. Only care about ASE in individual and comparing to ASE in non-hybrid
  final <- dat[keeps]
  #final <- head(final,((nrow(final))/3))
  #final <- head(final,1000)
  #head(final)
  final2 <- final
  final2["TISSUE"] <- "TIS2"
  #head(final2)
  final3 <- rbind(final,final2)
  #tail(final3)
  #head(final3)
  indiv_name <- paste("mamba_input_", file, sep = "")
  indiv_name <- paste(indiv_name, ".txt", sep = "")
  #write.table(final3, indiv_name, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  
  
  
  #### parse posterioir files ####
  
  setwd("c:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/axm/ASE/MAMBA/posteriors/")
  
  CRPA1_DNA <- read.table("CRPA1_DNA_posteriors.txt", header = TRUE, stringsAsFactors = FALSE)
  CRPM1_DNA <- read.table("CRPM1_DNA_posteriors.txt", header = TRUE, stringsAsFactors = FALSE)
  LILA1_DNA <- read.table("LILA1_DNA_posteriors.txt", header = TRUE, stringsAsFactors = FALSE)
  LILM3_DNA <- read.table("LILM3_DNA_posteriors.txt", header = TRUE, stringsAsFactors = FALSE)
  
  all_posterior_files <- c("axmE1","axmE2","axmJ1","axmJ4","axmJ5","axmJ6","CPAE1","CPAE2","CPAE3","CPAJ1","CPAJ2",
                           "CPAJ3","CPME1","CPME2","CPME3","CPMJ1","CPMJ2","CPMJ3","LLAE1","LLAE2","LLAE3","LLAJ1",
                           "LLAJ2","LLAJ3","LLME1","LLME2","LLME3","LLMJ1","LLMJ2","LLMJ3","CRPA1_DNA","CRPM1_DNA",
                           "LILA1_DNA","LILM3_DNA")
  all_posterior_files_rna <- c("axmE1","axmE2","axmJ1","axmJ4","axmJ5","axmJ6","CPAE1","CPAE2","CPAE3","CPAJ1","CPAJ2",
                               "CPAJ3","CPME1","CPME2","CPME3","CPMJ1","CPMJ2","CPMJ3","LLAE1","LLAE2","LLAE3","LLAJ1",
                               "LLAJ2","LLAJ3","LLME1","LLME2","LLME3","LLMJ1","LLMJ2","LLMJ3")
  hybrid_posterior_files <- c("axmE1","axmE2","axmJ1","axmJ4","axmJ5","axmJ6")
  parent_posterior_files <- c("CPAE1","CPAE2","CPAE3","CPAJ1","CPAJ2",
                              "CPAJ3","CPME1","CPME2","CPME3","CPMJ1","CPMJ2","CPMJ3","LLAE1","LLAE2","LLAE3","LLAJ1",
                              "LLAJ2","LLAJ3","LLME1","LLME2","LLME3","LLMJ1","LLMJ2","LLMJ3")
  
  hybrid_jaw_posterior_files <- c("axmJ1","axmJ4","axmJ5","axmJ6")
  parent_jaw_posterior_files <- c("CPAJ1","CPAJ2","CPAJ3","CPMJ1","CPMJ2","CPMJ3","LLAJ1","LLAJ2","LLAJ3","LLMJ1","LLMJ2","LLMJ3")
  hybrid_wbl_posterior_files <- c("axmE1","axmE2")
  parent_wbl_posterior_files <- c("CPAE1","CPAE2","CPAE3","CPME1","CPME3","LLAE1","LLAE2","LLAE3","LLME1","LLME2","LLME3")
  dna_control_poterior_files <- c("CRPA1_DNA","CRPM1_DNA","LILA1_DNA","LILM3_DNA")
  
  high_ase <- c()
  mod_ase <- c()
  any_ase <- c()
  indiv <- c()
  proportion_high_ase <- c()
  proportion_mod_ase <- c()
  proportion_any_ase <- c()
  for (file in all_posterior_files)
  {
    
    infile <- paste(file,"_posteriors.txt",sep="")  
    dat <- read.table(infile, header = TRUE, stringsAsFactors = FALSE)
    #set allele count threshold
    #dat["total_count"] <- (dat$refCount + dat$altCount)
    #dat <- dat[which(dat$total_count >=10),]
    h_ase <- (nrow(dat[which(dat$SNGASE>= 0.9),]))
    m_ase <- (nrow(dat[which(dat$MODASE>= 0.9),]))
    a_ase <- h_ase + m_ase
    proportion_high_ase <- c(proportion_high_ase, (h_ase/(nrow(dat))))
    proportion_mod_ase <- c(proportion_mod_ase, (m_ase/(nrow(dat))))
    proportion_any_ase <- c(proportion_any_ase, (a_ase/(nrow(dat))))
    high_ase <- c(high_ase, h_ase)
    mod_ase <- c(mod_ase, m_ase)
    any_ase <- c(any_ase, a_ase)
    indiv <- c(indiv, file)
    
  }
  ase_table <- do.call(rbind, Map(data.frame, indiv=indiv, alleles_showing_high_ase=high_ase, alleles_showing_moderate_ase=mod_ase,
                                  alleles_showing_moderate_to_high_ase=any_ase,
                                  proportion_alleles_showing_high_ase=proportion_high_ase,
                                  proportion_alleles_showing_moderate_ase=proportion_mod_ase,
                                  proportion_alleles_showing_moderate_to_high_ase=proportion_any_ase))
  
  barplot(ase_table$alleles_showing_high_ase, names = ase_table$indiv, cex.axis = .9,cex.names = 1, las=2)
  barplot(ase_table$proportion_alleles_showing_high_ase, names = ase_table$indiv, cex.axis = .9,cex.names = 1, las=2)
  barplot(ase_table$alleles_showing_moderate_ase, names = ase_table$indiv, cex.axis = .9,cex.names = 1, las=2)
  barplot(ase_table$proportion_alleles_showing_moderate_ase, names = ase_table$indiv, cex.axis = .9,cex.names = 1, las=2)        
  barplot(ase_table$alleles_showing_moderate_to_high_ase, names = ase_table$indiv, cex.axis = .9,cex.names = 1, las=2)
  barplot(ase_table$proportion_alleles_showing_moderate_to_high_ase, names = ase_table$indiv, cex.axis = .9,cex.names = 1, las=2)        
  
  head(high_ase_table)
  
  ###### remove biased alleles that show ase in parent DNA counts #####
  
  high_tags <- c()
  for (file in dna_control_poterior_files)
  {
    
    infile <- paste(file,"_posteriors.txt",sep="")  
    dat <- read.table(infile, header = TRUE, stringsAsFactors = FALSE)
    dat["total_count"] <- (dat$refCount + dat$altCount)
    dat <- dat[which(dat$total_count >=10),]
    #dat <- (dat[which(dat$SNGASE>= 0.9),])
    dat <- (dat[which(dat$MODASE>= 0.9 | dat$SNGASE>= 0.9),])
    high_tag <- as.vector((paste(dat$contig, dat$position, sep= ":")))
    high_tags <- c(high_tags, high_tag)
    
  }
  length(high_tags)
  high_ase_dna <- as.data.frame(table(high_tags))
  hist(high_ase_dna$Freq)
  snps_to_remove_dna <- high_ase_dna[which(high_ase_dna$Freq >= 2),]
  snps_to_remove_dna <- as.vector(snps_to_remove_dna$high_tags)
  #number of snps filtered that show high ase in dna_counts
  length(snps_to_remove_dna)
  test <- dat[order(dat$SNGASE, decreasing = TRUE),]
  head(test1,100)
  test1 <- test[which(test$refCount >=10 & test$altCount >= 10),]
  test$diff <- abs(test$refCount - test$altCount)
  hist(test1$diff, xlim = c(0,500), breaks = 1000)
  nrow(test1)
  nrow(test1[which(test1$diff < 2),])
  
  ###### remove biased alleles that show ase in parent RNA counts #####
  
  high_tags <- c()
  for (file in parent_jaw_posterior_files)
  {
    
    infile <- paste(file,"_posteriors.txt",sep="")  
    dat <- read.table(infile, header = TRUE, stringsAsFactors = FALSE)
    dat["total_count"] <- (dat$refCount + dat$altCount)
    dat <- dat[which(dat$total_count >=10),]
    #dat <- (dat[which(dat$SNGASE>= 0.9),])
    dat <- (dat[which(dat$MODASE>= 0.9 | dat$SNGASE>= 0.9),])
    tag <- as.vector((paste(dat$contig, dat$position, sep= ":")))
    high_tag <- as.vector((paste(dat$contig, dat$position, sep= ":")))
    high_tags <- c(high_tags, high_tag)
    
  }
  length(high_tags)
  high_ase_rna <- as.data.frame(table(high_tags))
  hist(high_ase_rna$Freq)
  snps_to_remove_rna <- high_ase_rna[which(high_ase_rna$Freq >= 2),]
  snps_to_remove_rna <- as.vector(snps_to_remove_rna$high_tags)
  #number of snps filtered that show high ase in rna_counts
  length(snps_to_remove_rna)
  
  
  high_tags <- c()
  for (file in hybrid_jaw_posterior_files)
  {
    
    infile <- paste(file,"_posteriors.txt",sep="")  
    dat <- read.table(infile, header = TRUE, stringsAsFactors = FALSE)
    dat["total_count"] <- (dat$refCount + dat$altCount)
    dat <- dat[which(dat$total_count >=10),]
    #dat <- (dat[which(dat$SNGASE>= 0.9),])
    dat <- (dat[which(dat$MODASE>= 0.9 | dat$SNGASE>= 0.9),])
    tag <- as.vector((paste(dat$contig, dat$position, sep= ":")))
    dat$tag <- tag
    dat_filtered_dna <- dat[ ! dat$tag %in% snps_to_remove_dna, ]
    dat_filtered_rna <- dat_filtered_dna[ ! dat_filtered_dna$tag %in% snps_to_remove_rna, ]
    high_tag <- as.vector((paste(dat_filtered_rna$contig, dat_filtered_rna$position, sep= ":")))
    high_tags <- c(high_tags, high_tag)
    
  }
  length(high_tags)
  high_ase_hybrids <- as.data.frame(table(high_tags))
  hist(high_ase_hybrids$Freq)
  snps_showing_high_ase <- high_ase_hybrids[which(high_ase_hybrids$Freq >= 2),]
  snps_showing_high_ase <- as.vector(snps_showing_high_ase$high_tags)
  #number of snps that show high ase exclusively in hybrids 
  length(snps_showing_high_ase)
  
  #### Overlap snps showing ase in hybrids with genes ####
  
  #library(reshape2)
  hybrid_ase <- as.data.frame(snps_showing_high_ase)
  hybrid_ase <- cbind(hybrid_ase, colsplit(hybrid_ase$snps_showing_high_ase, ":", c("genomic_accession", "position")))
  gene_locations <- read.table("C:/Users/Joseph McGirr Lab/Desktop/Cyprinodon/feature_table.txt", header = TRUE, stringsAsFactors = FALSE, quote = "", fill = TRUE, sep = "\t")
  head(gene_locations)
  head(hybrid_ase)
  
  scaffs <- hybrid_ase$genomic_accession
  
  scaffolds   <- c()
  start_s     <- c()
  stop_s      <- c()
  gene_names <- c()
  snp_pss     <- c()
  features   <- c()
  
  for (scaff in scaffs) 
  {
    
    snp_table   <- hybrid_ase[which(hybrid_ase$genomic_accession == scaff),]
    genes_table <- gene_locations[which(gene_locations$genomic_accession == scaff),]
    gene_strt   <- genes_table$start
    gene_stp    <- genes_table$end
    gene_name   <- genes_table$GeneID
    feature     <- genes_table$feature
    
    run <- c(1:length(gene_strt))
    if ((nrow(genes_table)) > 0 )
    {
      for (ps in snp_table$position)
      {
        for (i in run)
        {
          
          
          if (ps >= gene_strt[i] & ps <= gene_stp[i])
          {
            
            scaffolds <- c(scaffolds, scaff)
            start_s <- c(start_s, gene_strt[i])
            stop_s <- c(stop_s, gene_stp[i])
            gene_names <- c(gene_names, gene_name[i])
            snp_pss <- c(snp_pss, ps)
            features <- c(features, feature[i])
            
          }
          
        }}}}
  
  test5 <- data.frame(genomic_accession = scaffolds, start = start_s, end = stop_s, position = snp_pss, 
                      GeneID = gene_names, feature = features,
                      stringsAsFactors = FALSE)
  test6 <- merge(test5, gene_locations, by = c("genomic_accession", "start", "end", "feature", "GeneID"))
  test6$tag <- paste(test6$genomic_accession, test6$position, test6$start, test6$end,test6$feature, test6$GeneID, sep=":")
  hybrid_ase_genes <- subset(test6, !duplicated(test6[,22]))
  nrow(hybrid_ase_genes)
  length(unique(hybrid_ase_genes$GeneID))
  #write.table(hybrid_ase_genes, "c:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/axm/ASE/MAMBA/posteriors/results/hybrid_mod_to_sngase_embryo_genes.txt", sep = "\t", quote = FALSE, row.names = FALSE)
  
  #### cis regulatory candidates (ase in am hybrids and de am) ####
  
  all_DE <- read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/blast/all_DE_master.csv", header = TRUE, stringsAsFactors = FALSE)
  head(all_DE)
  head(hybrid_ase_genes)
  am_de <- all_DE[which(all_DE$comp == "am_embryo"),]
  am_de$genomic_accession <- am_de$Chr
  cis <- merge(hybrid_ase_genes, am_de, by = c("GeneID", "genomic_accession"))
  head(cis)
  length(cis)
  nrow(cis)
  cis <- subset(cis, !duplicated(cis[,22]))
  nrow(cis)
  #write.table(cis, "c:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/axm/ASE/MAMBA/posteriors/results/hybrid__mod_to_sngase_cis_candidate_embryo_genes.txt", sep = "\t", quote = FALSE, row.names = FALSE)
  
  
  #### Overlap fixed snps with genes showing ase in hybrids ####
  
  fixed_snp_genes <- read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/43_pupfish/fst/43_nomaysim/fixed_vxm_fst_merge.csv", header = TRUE, stringsAsFactors = FALSE)
  ase_genes <- read.table("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/axm/ASE/MAMBA/posteriors/results/hybrid__mod_to_sngase_cis_candidate_embryo_genes.txt", header = TRUE, stringsAsFactors = FALSE, sep = '\t')
  head(fixed_snp_genes)
  head(ase_genes)
  ase_genes$symbol <- ase_genes$symbol.x
  fixed_snp_ase_genes <- merge(fixed_snp_genes, ase_genes, by = c("symbol"))
  head(fixed_snp_ase_genes)
  fixed_snp_ase_genes <- subset(fixed_snp_ase_genes, !duplicated(fixed_snp_ase_genes['tag']))
  nrow(fixed_snp_ase_genes)
  
  
  
  
  
  
  
  
  test_tag <- as.vector((paste(CPAE1$contig, CPAE1$position, sep= ":")))
  CPAE1$tag <- test_tag
  CPAE1_test <- CPAE1[ ! CPAE1$tag %in% snps_to_remove, ]
  nrow(CPAE1_test)
  nrow(CPAE1)
  length(CPAE1$tag) - length(setdiff(test_tag, snps_to_remove))
  length(CPAE1$tag) - length(intersect(test_tag, snps_to_remove))
  #tagt <- as.vector((paste(CPAE1_test$contig, CPAE1_test$position, CPAE1_test$refAllele, sep= ":")))
  
  tag_freqs <- as.data.frame(table(tags))
  head(tag_freqs)
  hist(tag_freqs$Freq)
  sum(tag_freqs$Freq)
  nrow(tag_freqs[which(tag_freqs$Freq >= 10),])
  
  
  
  #
  
  
  
  
  
  

#### MAMBA DRAFTS ####
#### convert gatk ASEReadCounter datatables into mamba format ####

dat <- read.table("c:/Users/Joseph McGirr Lab/Desktop/ubuntushare/MAMBA/allele_counts/CPME1/CPME1_hetvcf_ASE_gatk.txt", header = TRUE, stringsAsFactors = FALSE)

#choose a test gene region
#dat <-  dat[which(dat$contig == "NW_015150787.1" & dat$position > 100000 & dat$position < 109000),]
dat <-  dat[which(dat$refCount > 0 & dat$altCount > 0),]
head(dat)
nrow(dat)
tail(dat)

dat["RSID"] <- paste(dat$contig, dat$position, sep = "_")
dat["INDIVIDUAL"] <- "ind1"
dat["REF_COUNT"] <- dat$refCount
dat["NONREF_COUNT"] <- dat$altCount
dat["TOTAL_COUNT"] <- dat$totalCount
dat["EXON_INFO"] <- "SEVERE_IMPACT=STOP_GAINED"
dat["TISSUE"] <- "TIS1"
keeps <- c("INDIVIDUAL", "RSID", "REF_COUNT", "NONREF_COUNT", "TOTAL_COUNT", "EXON_INFO", "TISSUE")

#create dataframe where tis1 and tis2 are identical for read counts. Only care about ASE in individual and comparing to ASE in non-hybrid
final <- dat[keeps]
#final <- head(final,((nrow(final))/3))
#final <- head(final,1000)
head(final)
final2 <- final
final2["TISSUE"] <- "TIS2"
head(final2)
final3 <- rbind(final,final2)
tail(final3)
head(final3)

#write.table(final3, "c:/Users/Joseph McGirr Lab/Desktop/ubuntushare/MAMBA/allele_counts/CPME1/CPME1_hets.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
nrow(final)
final
nrow(final3)

cpm <- final
cpa <- final
axm1 <- final



############# DRAFTS ##################################################################################################

setwd("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/axm/ASE/")

comp <- read.table("amE_v_axmE.diff.sites_in_files", header = TRUE, stringsAsFactors = FALSE)
head(comp)
nrow(comp)
unique(comp$IN_FILE)

hybrid_het_sites <- comp[which(comp$IN_FILE == "2"),]
nrow(hybrid_het_sites)
head(hybrid_het_sites)


keeps <- c("CHROM", "POS2")
excludes <- hybrid_het_sites[keeps]
head(excludes)
nrow(excludes)
nrow(comp)
nrow(comp) - nrow(excludes)
test_gene <-  excludes[which(excludes$CHROM == "NW_015150551.1"),]
excludes <- tail(test_gene,100)

#write.table(excludes, "axmE_dlx6_het_positions", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


genes <- read.table("C:/Users/Joseph McGirr Lab/Desktop/Cyprinodon/ref_C_variegatus-1.0_scaffolds_no_header.gff3", header = FALSE, stringsAsFactors = FALSE, sep = "\t")
head(genes)
genes <- cbind(genes, colsplit(genes[,9], "gbkey=", c("junk", "typejunk")))
genes <- cbind(genes, colsplit(genes$typejunk, ";", c("type", "junker")))
head(genes)
genes1 <- genes[which(genes$type == "mRNA"),]
nrow(genes1)
genes2 <- cbind(genes1, colsplit(genes1$junker, "gene=", c("junk", "genejunk")))
genes2 <- cbind(genes2, colsplit(genes2$genejunk, ";", c("gene", "junker")))
head(genes2)
keeps <- c("V1", "V4", "V5", "gene")
final <- genes2[keeps]
nrow(final)
head(final)


#write.table(final, "genes_mRNA.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


###########
#file paths#
setwd("C:/Users/Joseph McGirr Lab/Desktop/ubuntushare/MAMBA/allele_counts/")
gatk_files <- c("axmE1_hetvcf_ASE_gatk",
"axmE2_hetvcf_ASE_gatk",
"axmJ1_hetvcf_ASE_gatk",
"axmJ4_hetvcf_ASE_gatk",
"axmJ5_hetvcf_ASE_gatk",
"axmJ6_hetvcf_ASE_gatk",
"CPAE1_hetvcf_ASE_gatk",
"CPAE2_hetvcf_ASE_gatk",
"CPAE3_hetvcf_ASE_gatk",
"CPAJ1_hetvcf_ASE_gatk",
"CPAJ2_hetvcf_ASE_gatk",
"CPAJ3_hetvcf_ASE_gatk",
"CPME1_hetvcf_ASE_gatk",
"CPME2_hetvcf_ASE_gatk",
"CPME3_hetvcf_ASE_gatk",
"CPMJ1_hetvcf_ASE_gatk",
"CPMJ2_hetvcf_ASE_gatk",
"CPMJ3_hetvcf_ASE_gatk",
"CPPE1_hetvcf_ASE_gatk",
"CPPE2_hetvcf_ASE_gatk",
"CPPE3_hetvcf_ASE_gatk",
"CPPJ1_hetvcf_ASE_gatk",
"CPPJ2_hetvcf_ASE_gatk",
"LLAE1_hetvcf_ASE_gatk",
"LLAE2_hetvcf_ASE_gatk",
"LLAE3_hetvcf_ASE_gatk",
"LLAJ1_hetvcf_ASE_gatk",
"LLAJ2_hetvcf_ASE_gatk",
"LLAJ3_hetvcf_ASE_gatk",
"LLME1_hetvcf_ASE_gatk",
"LLME2_hetvcf_ASE_gatk",
"LLME3_hetvcf_ASE_gatk",
"LLMJ1_hetvcf_ASE_gatk",
"LLMJ2_hetvcf_ASE_gatk",
"LLMJ3_hetvcf_ASE_gatk")

for(file in gatk files) 
{
  infile <- paste(file,".txt",sep="")
  outfile <- paste(i,"-edit.txt",sep="")
  
  data <- read.table(infile,header=TRUE,sep="\t",row.names=NULL)
  dat["RSID"] <- paste(dat$contig, dat$position, sep = ":")
  dat["INDIVIDUAL"] <- "axmE1"
  dat["REF_COUNT"] <- dat$refCount
  dat["NONREF_COUNT"] <- dat$altCount
  dat["TOTAL_COUNT"] <- dat$totalCount
  dat["EXON_INFO"] <- "SEVERE_IMPACT=STOP_GAINED"
  dat["TISSUE"] <- "TIS1"
  keeps <- c("INDIVIDUAL", "RSID", "REF_COUNT", "NONREF_COUNT", "TOTAL_COUNT", "EXON_INFO", "TISSUE")
  
  final <- dat[keeps]
  #write.table(final,quote=FALSE,sep=", ",outfile)
}
 

 
dat <- read.table("c:/Users/Joseph McGirr Lab/Desktop/ubuntushare/MAMBA/allele_counts/axmE1_hetvcf_ASE_gatk.txt", header = TRUE, stringsAsFactors = FALSE)
head(dat)
dat["RSID"] <- paste(dat$contig, dat$position, sep = "_")
dat["INDIVIDUAL"] <- "axmE1"
dat["REF_COUNT"] <- dat$refCount
dat["NONREF_COUNT"] <- dat$altCount
dat["TOTAL_COUNT"] <- dat$totalCount
dat["EXON_INFO"] <- "SEVERE_IMPACT=STOP_GAINED"
dat["TISSUE"] <- "TIS1"
keeps <- c("INDIVIDUAL", "RSID", "REF_COUNT", "NONREF_COUNT", "TOTAL_COUNT", "EXON_INFO", "TISSUE")

final <- dat[keeps]
#final <- head(final,((nrow(final))/3))
#final <- head(final,1000)
head(final)
final2 <- final
final2["TISSUE"] <- "TIS2"
head(final2)
final3 <- rbind(final,final2)
head(final3,50)

#write.table(final3, "c:/Users/Joseph McGirr Lab/Desktop/ubuntushare/MAMBA/allele_counts/test/test2/test2_mamba.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
nrow(final)
nrow(final3)

#### make a test file ####

dat <- read.table("c:/Users/Joseph McGirr Lab/Desktop/ubuntushare/MAMBA/allele_counts/axmE1_hetvcf_ASE_gatk.txt", header = TRUE, stringsAsFactors = FALSE)
head(dat)
dat <- dat[which(dat$refCount != 0 & dat$altCount != 0),]
dat["RSID"] <- paste(dat$contig, dat$position, sep = "_")
dat["INDIVIDUAL"] <- "IND1"
dat["REF_COUNT"] <- dat$refCount
dat["NONREF_COUNT"] <- dat$altCount
dat["TOTAL_COUNT"] <- dat$totalCount
dat["EXON_INFO"] <- "SEVERE_IMPACT=STOP_GAINED"
dat["TISSUE"] <- "TIS1"
keeps <- c("INDIVIDUAL", "RSID", "REF_COUNT", "NONREF_COUNT", "TOTAL_COUNT", "EXON_INFO", "TISSUE")

final <- dat[keeps]
nrow(final)


#write.table(final, "c:/Users/Joseph McGirr Lab/Desktop/ubuntushare/MAMBA/allele_counts/test/test2/test2_mamba.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

ase <- read.table("c:/Users/Joseph McGirr Lab/Desktop/ubuntushare/MAMBA/test/asetest.txt", header = TRUE, stringsAsFactors = FALSE)
head(ase)

unique(ase$RSID)
unique(ase$INDIVIDUAL)
unique(ase$REF_COUNT)
unique(ase$NONREF_COUNT)
unique(ase$TOTAL_COUNT)
unique(ase$EXON_INFO)
unique(ase$TISSUE)
