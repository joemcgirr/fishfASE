#------------------------------------------------------------------------------------
#
# mbased.R
#
# This script is used after parental_allele_counts.py
# to estimate snp or gene level alelle specific expression
#
# example input:
# https://github.com/joemcgirr/fishfASE/blob/master/examples/CUT1_parental_counts.txt
# 
# example output:
# https://github.com/joemcgirr/fishfASE/blob/master/examples/CUT1_mbased_ase_gene_level.txt
#
#------------------------------------------------------------------------------------

library("MBASED")
library("dplyr")
library("magrittr")
library("plyranges")

# F1 rna sample name
indiv <- "CUT1" 

# directories
out_dir <- "/pine/scr/j/m/jmcgirr/pupfish_transcriptomes/ase/mbased/gene_output/"
cts_dir <- "/pine/scr/j/m/jmcgirr/pupfish_transcriptomes/ase/allele_counts/"
gene <- read.table("/pine/scr/j/m/jmcgirr/pupfish_transcriptomes/ase/mbased/c_brontotheroides.all.renamed.putative_function.genes_only_reformated_known_final.bed", header = FALSE, stringsAsFactors = FALSE, sep = "	")

# parental counts output by parental_allele_counts.py
parental_cts <- read.table(paste(cts_dir, indiv, "_parental_counts.txt",sep = ""),sep = "	", stringsAsFactors = FALSE, header = TRUE)

parental_cts_ranges <- parental_cts
parental_cts_ranges$seqnames <- parental_cts_ranges$chrom
parental_cts_ranges$start <- parental_cts_ranges$position
parental_cts_ranges$end <- parental_cts_ranges$start +1
parental_cts_ranges <- parental_cts_ranges %>% as_granges()

names(gene) <- c("seqnames", "start","end","GeneID", "strand")
gene <- gene %>% as_granges()
gene <- join_overlap_intersect(gene, parental_cts_ranges) %>% as.data.frame()

gene$geneIndex <- gene$GeneID
gene$alleleOneCount <- gene$refCount
gene$alleleTwoCount <- gene$altCount
gene$aseIndex <- paste(gene$geneIndex, gene$snpIndex, sep = ":")
gene$zero_counts <- paste(gene$refCount, gene$altCount, sep= ":")
gene <- gene %>% filter(zero_counts != "0:0")
cts <- gene[c("geneIndex", "snpIndex","alleleOneCount","alleleTwoCount","refAllele","altAllele", "chrom","position","aseIndex")]

# change aseID=cts$aseIndex to aseID=cts$geneIndex for gene/exon level ase rather than snp level
set.seed(988482)
mySNVs <- GRanges(
  seqnames=cts$chrom,ranges=IRanges(start=cts$position, width=1),
  aseID=cts$geneIndex,allele1=cts$refAllele,allele2=cts$altAllele)
names(mySNVs) <- cts$geneIndex
## create input RangedSummarizedExperiment object
mySample <- SummarizedExperiment(assays=list(
  lociAllele1Counts=matrix(cts$alleleOneCount,ncol=1,dimnames=list(names(mySNVs),
  "mySample")),lociAllele2Counts=matrix(cts$alleleTwoCount,ncol=1,dimnames=list(names(mySNVs),"mySample"))),rowRanges=mySNVs)

ASEresults_1s_haplotypesKnown <- runMBASED(
  ASESummarizedExperiment=mySample,
  isPhased=TRUE,
  numSim=10^5,
  BPPARAM = SerialParam()
)

class(ASEresults_1s_haplotypesKnown)
names(assays(ASEresults_1s_haplotypesKnown))
head(assays(ASEresults_1s_haplotypesKnown)$majorAlleleFrequency)
head(assays(metadata(ASEresults_1s_haplotypesKnown)$locusSpecificResults)$allele1IsMajor)

summarizeASEResults_1s <- function(MBASEDOutput) {
  geneOutputDF <- data.frame(
    majorAlleleFrequency=assays(MBASEDOutput)$majorAlleleFrequency[,1],
    #allele1IsMajor=assays(metadata(MBASEDOutput)$locusSpecificResults)$allele1IsMajor,
    pValueASE=assays(MBASEDOutput)$pValueASE[,1],
    pValueHeterogeneity=assays(MBASEDOutput)$pValueHeterogeneity[,1]
  )
  lociOutputGR <- rowRanges(metadata(MBASEDOutput)$locusSpecificResults)
  lociOutputGR$allele1IsMajor <- assays(metadata(MBASEDOutput)$locusSpecificResults)$allele1IsMajor[,1]
  lociOutputGR$MAF <- assays(metadata(MBASEDOutput)$locusSpecificResults)$MAF[,1]
  lociOutputList <- split(lociOutputGR, factor(lociOutputGR$aseID, levels=unique(lociOutputGR$aseID)))
  return(
    list(
      geneOutput=geneOutputDF,
      locusOutput=lociOutputList
    )
  )
}

ase <- as.data.frame(summarizeASEResults_1s(ASEresults_1s_haplotypesKnown)$geneOutput)
ase$geneID <- rownames(ase)

write.table(ase,paste(out_dir,indiv,"_mbased_ase_phased_",ase_run_type,".txt",sep = ""), quote = FALSE, row.names = FALSE, sep = "	")