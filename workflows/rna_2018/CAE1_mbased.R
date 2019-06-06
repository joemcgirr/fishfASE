ase_dir <- "/pine/scr/j/m/jmcgirr/mbased/"
cts_dir <- "/pine/scr/j/m/jmcgirr/wasp/vcf/allele_counts/maternal_counts/"
indiv <- "CAE1"

#source("http://bioconductor.org/biocLite.R")
#biocLite("MBASED")
library("MBASED")

cts <- read.table(paste(cts_dir, indiv, "maternal_counts_features.txt",sep = ""),sep = "	",fill = TRUE, stringsAsFactors = FALSE, header = TRUE)
nrow(cts)
length(unique(cts$mrnaID))
nrow(cts[which(cts$geneID == ""),])
cts <- cts[which(cts$geneID != ""),]
nrow(cts)
head(cts)
## set gene or transcript isoform as locus
cts$locus <- paste(cts$mrnaID,cts$snpIndex, sep = ";")
#cts$locus <- paste(cts$geneID,snpIndex, sep = ";")

set.seed(988482)
mySNVs <- GRanges(
 seqnames=cts$chrom,ranges=IRanges(start=cts$position, width=1),
 aseID=cts$mrnaID,allele1=cts$refAllele,allele2=cts$altAllele)
names(mySNVs) <- cts$locus
## create input RangedSummarizedExperiment object
mySample <- SummarizedExperiment(assays=list(
  lociAllele1Counts=matrix(cts$refCount,ncol=1,dimnames=list(names(mySNVs),
        'mySample')),lociAllele2Counts=matrix(
      cts$altCount,ncol=1,dimnames=list(names(mySNVs),
        'mySample'))),rowRanges=mySNVs)

ASEresults_1s_haplotypesKnown <- runMBASED(
ASESummarizedExperiment=mySample,
isPhased=TRUE,
numSim=10^6,
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
tail(summarizeASEResults_1s(ASEresults_1s_haplotypesKnown)$locusOutput,5)
ase <- as.data.frame(summarizeASEResults_1s(ASEresults_1s_haplotypesKnown)$geneOutput)
head(ase)
nrow(ase)
ase$mrnaID <- rownames(ase)
write.table(ase,paste(ase_dir,indiv,"_mbased_ase.txt",sep = ""), quote = FALSE, row.names = FALSE, sep = "	")

