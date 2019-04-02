library(DESeq2)
library(reshape2)
library(seqinr)
library(plyr)
library(MASS)
library(AER)
library(rlang)

all_cts <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/ase_data/all_samples_2018_counts_mrna_saf.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")


pops1 <- "OSPA"
pops2 <- "OSPP"

pops1 <- c("OSPA","OSPP")
pops2 <- "OAxOP"

stage <- "8dpf"


setwd("C:/Users/jmcgirr/Documents/all_2018_samples/conditions/")
master <- read.table("D:/Martin Lab/rna_2018/all_2018_samples/table_maker_master_outlier_rm.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
head(master)
pops1_table <- master[master$f1 %in% pops1, ]
pops2_table <- master[master$f1 %in% pops2, ]
pops1_table$species <- "a" 
pops2_table$species <- "b"
comp_table <- rbind(pops1_table,pops2_table)
comp_table <- comp_table[which(comp_table$stage == stage),]
sample_size_a <- nrow(comp_table[which(comp_table$species == "a"),])
sample_size_b <- nrow(comp_table[which(comp_table$species == "b"),])
sample_size_a
sample_size_b
write.table(comp_table, "C:/Users/jmcgirr/Documents/test_delete.txt", row.names = FALSE, quote= FALSE,sep="\t")
#pops1 <- "OSPA_and_OSPP"
comp_file <- paste(pops1, "_vs_", pops2, "_" ,stage, sep = "")
for (i in c(1))
{
  

setwd("C:/Users/jmcgirr/Documents/all_2018_samples/conditions/")
col_data <-             paste("condition_species_", comp_file, ".txt", sep = "")
cts_data <-             paste("DESeq_counts_", comp_file, ".txt", sep = "")
cts_data_genes <-       paste("DE_", comp_file, ".csv", sep = "")
genes_out <-            paste("DE_", comp_file, "_genes.csv", sep = "")

sample_list <- read.table(col_data, header = TRUE, stringsAsFactors = FALSE)
samples_a <- sample_list[which(sample_list$species == "a"),]
samples_b <- sample_list[which(sample_list$species == "b"),]
add_sequencing_round_to_model <- length(unique(samples_a$sequencing_round)) + length(unique(samples_b$sequencing_round))
keeps <- c("Geneid", sample_list$sample)
keeper <- sample_list$sample

cts <- all_cts[keeps]
#head(cts)
write.table(cts, "C:/Users/jmcgirr/Documents/test_delete.txt", row.names = FALSE, quote= FALSE,sep="\t") 

### DIFFERENTIAL EXPRESSION ANALYSIS ###
#vignette("DESeq2")
cts <- as.matrix(read.table("C:/Users/jmcgirr/Documents/test_delete.txt" ,sep = "\t",header = TRUE,row.names=1))
head(cts)
nrow(cts)
colData <- as.matrix(read.table(col_data ,header = TRUE,row.names=1))
head(colData)
ncol(cts)
nrow(colData)


#if (add_sequencing_round_to_model > 3)
#{
# dds <- DESeqDataSetFromMatrix(countData = cts,
#                                colData = colData,
#                                design= ~species+sequencing_round)
#}else 
#{
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = colData,
                              design= ~species)
#}
### count number of unused genes
#unused <- dds[ rowSums(counts(dds)) < 1, ]
#unused
# most stringent -- filter out genes if read count is zero for more than half of individuals.

dds <- estimateSizeFactors(dds)
idx <- rowSums(counts(dds, normalized=TRUE) >= 10 ) >= (nrow(colData))
#length(rowSums(counts(dds, normalized=TRUE) >= 1 ) >= (nrow(colData)/2))
length(idx)
dds <- dds[idx,]
dds <- DESeq(dds)

# evolutionletters filter
#length(rowSums(counts(dds)) > 1)
#dds <- dds[ rowSums(counts(dds)) > 1, ]
#dds <- DESeq(dds)
#res <- results(dds, independentFiltering=FALSE)
#res$pvalue[res$baseMean < 10] <- NA
#res$padj <- p.adjust(res$pvalue, method="BH")

res <- results(dds, alpha=0.05)
#res <- results(dds, contrast=c("f1","CRPA","NCA"))
resOrdered <- res[order(res$padj),]
summary(res)
resLFC <- lfcShrink(dds, coef=2, res=res)

res_ordered <- as.data.frame(resOrdered)
res_ordered$Geneid <- rownames(res_ordered)
nrow(res_ordered)
head(res_ordered)
length(unique(res_ordered$Geneid))
options(scipen=999)
#plotCounts(dds,gene="gene678", intgroup='species')

total_genes <- nrow(res_ordered)
de_total <- nrow(res_ordered[which(res_ordered$padj < 0.05),])
de_up <- (nrow(res_ordered[which(res_ordered$log2FoldChange > 0 & res_ordered$padj < 0.05),]))/total_genes
de_dn <- (nrow(res_ordered[which(res_ordered$log2FoldChange < 0 & res_ordered$padj < 0.05),]))/total_genes
prop_de <- de_total/total_genes
sample_sizes_plot <- paste(sample_size_a, " vs ", sample_size_b, sep = "")
total_genes_plot <- paste(total_genes, "transcripts", sep = " ")
}
de_plot <- paste("C:/Users/jmcgirr/Documents/all_2018_samples/manuscript_figs/",comp_file, "_de_plot.tiff", sep = "")
de_total
#tiff(de_plot, width = 3.5, height = 3.5, units = 'in', res = 1000)
plotMA(resLFC, ylim=c(-4.7,4.3), main = "", colSig = blu, colNonSig = "grey",
       yaxt = 'n', cex = 0.58, xlab = "", ylab = "")#,xaxt = 'n')
#axis(1,c(10,1000,100000))
axis(2,c(-4,-2,0,2,4))
#dev.off()
#legend("topright", legend=c(sample_sizes_plot,total_genes_plot, de_total_plot, prop_de),cex=1.0, bty = 'n')
#legend("bottomleft", legend=c(paste((100*(round(de_dn, digits = 3))), "% DE down", sep = "")),cex=0.8, bty = 'n')
#legend("topleft", legend=c(paste((100*(round(de_up, digits = 3))), "% DE up", sep = "")),cex=0.8, bty = 'n')
#legend('bottomright', legend=c(boot_ci, dns_ci),cex=0.8, bty = 'n')
#dev.off()





###############################################
###############################################
############ visualize ase at cool genes ######
###############################################

library(reshape2)
library(seqinr)

norm_cts <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/norm_cts_all_2018_no_seq_round_control.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
inheritance_comps <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/ase_data/ase_comparisons.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
master <- read.table("D:/Martin Lab/rna_2018/all_2018_samples/table_maker_master_outlier_rm.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
norm_cts <- cbind(norm_cts, colsplit(norm_cts$Geneid, ";", c("related_accession", "gene_name")))
mrna <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/ase_data/mrna.saf", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
mrna <- cbind(mrna, colsplit(mrna$GeneID, ";", c("related_accession", "gene_name")))
maternal_counts_dir <- "C:/Users/jmcgirr/Documents/all_2018_samples/ase_data/maternal_counts/"
gff <- read.delim("D:/Cyprinodon/ref_C_variegatus-1.0_scaffolds_no_header.gff3",header = FALSE, stringsAsFactors = FALSE, sep = "\t")
cool_genes_out_dir <- "C:/Users/jmcgirr/Documents/all_2018_samples/manuscript_figs/can_genes/"

# itga 8 opxom parallel
#cool_gene <- "XM_015370652.1"
# sech24d cmxcp adaptive mis
#cool_gene <- "XM_015397247.1"


cool_genes_comps <- read.table("D:/test_plots/ai_fst_dxy_gemma.txt", header = TRUE, stringsAsFactors = FALSE)
i <- 1
pdf(paste("D:/test_plots/tests_",cool_genes_comps$interesting_because[i],".pdf",sep = ""),width=7,height=5)
for (i in c(1:nrow(cool_genes_comps)))
{
  cool_stage <- cool_genes_comps$stage[i]
  pop_gen_comp <- cool_genes_comps$pops[i]
  cool_p1 <- cool_genes_comps$p1[i]
  cool_p2 <- cool_genes_comps$p2[i]
  cool_h <- cool_genes_comps$h[i]
  cool_cols <- c(gre,blu,grb)
  
  p1_inds <- master[which(master$f1 == cool_p1 & master$stage == cool_stage),]
  p2_inds <- master[which(master$f1 == cool_p2 & master$stage == cool_stage),]
  hy_inds <- master[which(master$f1 == cool_h & master$stage == cool_stage),]
  p1_inds <- p1_inds$sample
  p2_inds <- p2_inds$sample
  hy_inds <- hy_inds$sample
  
  p1 <- strsplit(pop_gen_comp, "x")[[1]][1]
  p2 <- strsplit(pop_gen_comp, "x")[[1]][2]
  dxy_main <- read.csv(paste("D:/Martin Lab/rna_2018/fst/dna/",pop_gen_comp,"_popgen_dna_stats_corr_dxy.csv",sep = ""), header = TRUE, stringsAsFactors = FALSE)

  fst_main <- read.table(paste("C:/Users/jmcgirr/Documents/all_2018_samples/fst_dna/",pop_gen_comp,"_fst",sep = ""), header = TRUE, stringsAsFactors = FALSE)
  fst_main <- fst_main[which(fst_main$WEIR_AND_COCKERHAM_FST >=0 & fst_main$WEIR_AND_COCKERHAM_FST <=1),]
  lake <- strsplit(pop_gen_comp,split ="")[[1]][1]
  taj_m_main <-    read.table(paste("C:/Users/jmcgirr/Documents/all_2018_samples/fst_dna/",lake,"m_taj_d_20kb.txt.Tajima.D", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  taj_a_main <-    read.table(paste("C:/Users/jmcgirr/Documents/all_2018_samples/fst_dna/",lake,"a_taj_d_20kb.txt.Tajima.D", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  taj_p_main <-    read.table(paste("C:/Users/jmcgirr/Documents/all_2018_samples/fst_dna/",lake,"p_taj_d_20kb.txt.Tajima.D", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  taj_m_main <- na.omit(taj_m_main)
  taj_a_main <- na.omit(taj_a_main)
  taj_p_main <- na.omit(taj_p_main)
  sweed_p1_main <-    read.table(paste("C:/Users/jmcgirr/Documents/all_2018_samples/sweed/correct_seq/SweeD_Report.",p1,"_pop_bottle_58_grid_500_sweeps.txt", sep = ""), header = FALSE, stringsAsFactors = FALSE)
  sweed_p2_main <-    read.table(paste("C:/Users/jmcgirr/Documents/all_2018_samples/sweed/correct_seq/SweeD_Report.",p2,"_pop_bottle_58_grid_500_sweeps.txt", sep = ""), header = FALSE, stringsAsFactors = FALSE)
  
  cool_genes <- strsplit(cool_genes_comps$genes[i],";")[[1]]
  for (cool_gene in cool_genes)
  {
  gene_chrom <- mrna[which(mrna$related_accession == cool_gene),]
  sweed_p1 <- sweed_p1_main[which(sweed_p1_main$V1 ==gene_chrom$Chr[1]),]
  sweed_p2 <- sweed_p2_main[which(sweed_p2_main$V1 ==gene_chrom$Chr[1]),]
  fst <- fst_main[which(fst_main$CHROM == gene_chrom$Chr[1]),]
  fixed <- fst[which(fst$WEIR_AND_COCKERHAM_FST == 1),]
  smooth_fst <-   smooth.spline(fst$POS, fst$WEIR_AND_COCKERHAM_FST, spar = .2)
  dxy <- dxy_main[which(dxy_main$scaffold == gene_chrom$Chr[1]),]
  
  smooth_dxy <-   smooth.spline(dxy$start, dxy$corr_dxy, spar = .2)
  smooth_p1_pi <- smooth.spline(dxy$start, dxy[,14], spar = .2)
  smooth_p2_pi <- smooth.spline(dxy$start, dxy[,16], spar = .2)
  
  tj_m <- taj_m_main[which(taj_m_main$CHROM == gene_chrom$Chr[1]),]
  smooth_m <-   smooth.spline(tj_m$BIN_START, tj_m$TajimaD, spar = .1)
  tj_p <- taj_p_main[which(taj_p_main$CHROM == gene_chrom$Chr[1]),]
  smooth_p <-   smooth.spline(tj_p$BIN_START, tj_p$TajimaD, spar = .1)
  tj_a <- taj_a_main[which(taj_a_main$CHROM == gene_chrom$Chr[1]),]
  smooth_a <-   smooth.spline(tj_a$BIN_START, tj_a$TajimaD, spar = .1)
  

  plot(fixed$POS, fixed$WEIR_AND_COCKERHAM_FST,ylim = c(0,1),
        xlab = "", ylab = "fst", yaxt = 'n',xlim =c((gene_chrom$Start[1] -500000),(gene_chrom$End[1] +500000)), main = paste(cool_p1, "vs", cool_p2,cool_gene, sep = " "))
  axis(2,c(0,0.5,1))
  lines(smooth_fst, lwd = 2, col = "black")
  rect(gene_chrom$Start[1], -1000, gene_chrom$End[1], 1000, border = NA, col = col2alpha(blu,0.2))
  plot(dxy$start, dxy$corr_dxy,col = "white",
        xlab = "", ylab = "dxy",xlim =c((gene_chrom$Start[1] -500000),(gene_chrom$End[1] +500000)))#, main = paste(cool_p1, "vs", cool_p2, sep = " "))
  lines(smooth_dxy, lwd = 2, col = "grey")
  rect(gene_chrom$Start[1], -1000, gene_chrom$End[1], 1000, border = NA, col = col2alpha(blu,0.2))
  plot(dxy$start, dxy[,14],col = "white",
       xlab = "", ylab = "pi",xlim =c((gene_chrom$Start[1] -500000),(gene_chrom$End[1] +500000)))#, main = paste(cool_p1, "vs", cool_p2, sep = " "))
  lines(smooth_p1_pi, lwd = 2, col = red)
  lines(smooth_p2_pi, lwd = 2, col = blu)
  rect(gene_chrom$Start[1], -1000, gene_chrom$End[1], 1000, border = NA, col = col2alpha(blu,0.2))

  plot(tj_p$BIN_START, tj_p$TajimaD,ylab = "D",ylim = c(-2.5,2.5),xlab = "",
       xlim = c((gene_chrom$Start[1] -500000),(gene_chrom$End[1] +500000)),col = "white", yaxt = "n")
  axis(2,c(-2,-1,0,1,2))
  lines(smooth_m, lwd = 2, col = gre)
  lines(smooth_p, lwd = 2, col = blu)
  lines(smooth_a, lwd = 2, col = red)
  rect(gene_chrom$Start[1], -1000, gene_chrom$End[1], 1000, border = NA, col = col2alpha(blu,0.2))

  sweed_p1$mid <- sweed_p1$V2 + ((sweed_p1$V2[2]-sweed_p1$V2[1])/2)
  sweed_p2$mid <- sweed_p2$V2 + ((sweed_p2$V2[2]-sweed_p2$V2[1])/2)
  sweed_p1$norm_clr = (sweed_p1$V3-min(sweed_p1$V3))/(max(sweed_p1$V3)-min(sweed_p1$V3))
  sweed_p2$norm_clr = (sweed_p2$V3-min(sweed_p2$V3))/(max(sweed_p2$V3)-min(sweed_p2$V3))
  smooth_sweed_p1 <-   smooth.spline(sweed_p1$mid, sweed_p1$norm_clr, spar = 0.2)
  smooth_sweed_p2 <-   smooth.spline(sweed_p2$mid, sweed_p2$norm_clr, spar = 0.2)
  
  plot(sweed_p1$mid, sweed_p1$norm_clr,xlab = "",ylab = "CLR",col = "white",xlim = c((gene_chrom$Start[1] -500000),(gene_chrom$End[1] +500000)))
  lines(smooth_sweed_p1, lwd = 2, col = red)
  rect(gene_chrom$Start[1], -1000, gene_chrom$End[1], 1000, border = NA, col = col2alpha(blu,0.2))
  rect(gene_chrom$Start[1], -1, gene_chrom$End[1], 1, border = NA, col = col2alpha(blu,0.2))
  plot(sweed_p2$mid, sweed_p2$norm_clr,xlab = "",ylab = "CLR",col = "white",xlim = c((gene_chrom$Start[1] -500000),(gene_chrom$End[1] +500000)))
  lines(smooth_sweed_p2, lwd = 2, col = blu)
  rect(gene_chrom$Start[1], -1000, gene_chrom$End[1], 1000, border = NA, col = col2alpha(blu,0.2))
  rect(gene_chrom$Start[1], -1, gene_chrom$End[1], 1, border = NA, col = col2alpha(blu,0.2))
  
  }
  
}
dev.off()

cool_genes_comps <- read.table("D:/test_plots/ai_fst_dxy.txt", header = TRUE, stringsAsFactors = FALSE)
i <- 1
pdf(paste("D:/test_plots/tests_",cool_genes_comps$interesting_because[i],".pdf",sep = ""),width=7,height=5)
for (i in c(1:nrow(cool_genes_comps)))
{
  cool_stage <- cool_genes_comps$stage[i]
  pop_gen_comp <- cool_genes_comps$pops[i]
  cool_p1 <- cool_genes_comps$p1[i]
  cool_p2 <- cool_genes_comps$p2[i]
  cool_h <- cool_genes_comps$h[i]
  cool_cols <- c(gre,blu,grb)
  
  p1_inds <- master[which(master$f1 == cool_p1 & master$stage == cool_stage),]
  p2_inds <- master[which(master$f1 == cool_p2 & master$stage == cool_stage),]
  hy_inds <- master[which(master$f1 == cool_h & master$stage == cool_stage),]
  p1_inds <- p1_inds$sample
  p2_inds <- p2_inds$sample
  hy_inds <- hy_inds$sample
  
  p1 <- strsplit(pop_gen_comp, "x")[[1]][1]
  p2 <- strsplit(pop_gen_comp, "x")[[1]][2]
  dxy_main <- read.csv(paste("D:/Martin Lab/rna_2018/fst/dna/",pop_gen_comp,"_popgen_dna_stats_corr_dxy.csv",sep = ""), header = TRUE, stringsAsFactors = FALSE)
  
  fst_main <- read.table(paste("C:/Users/jmcgirr/Documents/all_2018_samples/fst_dna/",pop_gen_comp,"_fst",sep = ""), header = TRUE, stringsAsFactors = FALSE)
  fst_main <- fst_main[which(fst_main$WEIR_AND_COCKERHAM_FST >=0 & fst_main$WEIR_AND_COCKERHAM_FST <=1),]
  lake <- strsplit(pop_gen_comp,split ="")[[1]][1]
  taj_m_main <-    read.table(paste("C:/Users/jmcgirr/Documents/all_2018_samples/fst_dna/",lake,"m_taj_d_20kb.txt.Tajima.D", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  taj_a_main <-    read.table(paste("C:/Users/jmcgirr/Documents/all_2018_samples/fst_dna/",lake,"a_taj_d_20kb.txt.Tajima.D", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  taj_p_main <-    read.table(paste("C:/Users/jmcgirr/Documents/all_2018_samples/fst_dna/",lake,"p_taj_d_20kb.txt.Tajima.D", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  taj_m_main <- na.omit(taj_m_main)
  taj_a_main <- na.omit(taj_a_main)
  taj_p_main <- na.omit(taj_p_main)
  sweed_p1_main <-    read.table(paste("C:/Users/jmcgirr/Documents/all_2018_samples/sweed/correct_seq/SweeD_Report.",p1,"_pop_bottle_58_grid_500_sweeps.txt", sep = ""), header = FALSE, stringsAsFactors = FALSE)
  sweed_p2_main <-    read.table(paste("C:/Users/jmcgirr/Documents/all_2018_samples/sweed/correct_seq/SweeD_Report.",p2,"_pop_bottle_58_grid_500_sweeps.txt", sep = ""), header = FALSE, stringsAsFactors = FALSE)
  
  cool_genes <- strsplit(cool_genes_comps$genes[i],";")[[1]]
  for (cool_gene in cool_genes)
  {
    gene_chrom <- mrna[which(mrna$related_accession == cool_gene),]
    sweed_p1 <- sweed_p1_main[which(sweed_p1_main$V1 ==gene_chrom$Chr[1]),]
    sweed_p2 <- sweed_p2_main[which(sweed_p2_main$V1 ==gene_chrom$Chr[1]),]
    fst <- fst_main[which(fst_main$CHROM == gene_chrom$Chr[1]),]
    fixed <- fst[which(fst$WEIR_AND_COCKERHAM_FST == 1),]
    smooth_fst <-   smooth.spline(fst$POS, fst$WEIR_AND_COCKERHAM_FST, spar = .2)
    dxy <- dxy_main[which(dxy_main$scaffold == gene_chrom$Chr[1]),]
    
    smooth_dxy <-   smooth.spline(dxy$start, dxy$corr_dxy, spar = .2)
    smooth_p1_pi <- smooth.spline(dxy$start, dxy[,14], spar = .2)
    smooth_p2_pi <- smooth.spline(dxy$start, dxy[,16], spar = .2)
    
    tj_m <- taj_m_main[which(taj_m_main$CHROM == gene_chrom$Chr[1]),]
    smooth_m <-   smooth.spline(tj_m$BIN_START, tj_m$TajimaD, spar = .1)
    tj_p <- taj_p_main[which(taj_p_main$CHROM == gene_chrom$Chr[1]),]
    smooth_p <-   smooth.spline(tj_p$BIN_START, tj_p$TajimaD, spar = .1)
    tj_a <- taj_a_main[which(taj_a_main$CHROM == gene_chrom$Chr[1]),]
    smooth_a <-   smooth.spline(tj_a$BIN_START, tj_a$TajimaD, spar = .1)
    
    
    plot(fixed$POS, fixed$WEIR_AND_COCKERHAM_FST,ylim = c(0,1),
         xlab = "", ylab = "fst", yaxt = 'n',xlim =c((gene_chrom$Start[1] -500000),(gene_chrom$End[1] +500000)), main = paste(cool_p1, "vs", cool_p2,cool_gene, sep = " "))
    axis(2,c(0,0.5,1))
    lines(smooth_fst, lwd = 2, col = "black")
    rect(gene_chrom$Start[1], -1000, gene_chrom$End[1], 1000, border = NA, col = col2alpha(blu,0.2))
    plot(dxy$start, dxy$corr_dxy,col = "white",
         xlab = "", ylab = "dxy",xlim =c((gene_chrom$Start[1] -500000),(gene_chrom$End[1] +500000)))#, main = paste(cool_p1, "vs", cool_p2, sep = " "))
    lines(smooth_dxy, lwd = 2, col = "grey")
    rect(gene_chrom$Start[1], -1000, gene_chrom$End[1], 1000, border = NA, col = col2alpha(blu,0.2))
    plot(dxy$start, dxy[,14],col = "white",
         xlab = "", ylab = "pi",xlim =c((gene_chrom$Start[1] -500000),(gene_chrom$End[1] +500000)))#, main = paste(cool_p1, "vs", cool_p2, sep = " "))
    lines(smooth_p1_pi, lwd = 2, col = red)
    lines(smooth_p2_pi, lwd = 2, col = blu)
    rect(gene_chrom$Start[1], -1000, gene_chrom$End[1], 1000, border = NA, col = col2alpha(blu,0.2))
    
    plot(tj_p$BIN_START, tj_p$TajimaD,ylab = "D",ylim = c(-2.5,2.5),xlab = "",
         xlim = c((gene_chrom$Start[1] -500000),(gene_chrom$End[1] +500000)),col = "white", yaxt = "n")
    axis(2,c(-2,-1,0,1,2))
    lines(smooth_m, lwd = 2, col = gre)
    lines(smooth_p, lwd = 2, col = blu)
    lines(smooth_a, lwd = 2, col = red)
    rect(gene_chrom$Start[1], -1000, gene_chrom$End[1], 1000, border = NA, col = col2alpha(blu,0.2))
    
    sweed_p1$mid <- sweed_p1$V2 + ((sweed_p1$V2[2]-sweed_p1$V2[1])/2)
    sweed_p2$mid <- sweed_p2$V2 + ((sweed_p2$V2[2]-sweed_p2$V2[1])/2)
    sweed_p1$norm_clr = (sweed_p1$V3-min(sweed_p1$V3))/(max(sweed_p1$V3)-min(sweed_p1$V3))
    sweed_p2$norm_clr = (sweed_p2$V3-min(sweed_p2$V3))/(max(sweed_p2$V3)-min(sweed_p2$V3))
    smooth_sweed_p1 <-   smooth.spline(sweed_p1$mid, sweed_p1$norm_clr, spar = 0.2)
    smooth_sweed_p2 <-   smooth.spline(sweed_p2$mid, sweed_p2$norm_clr, spar = 0.2)
    
    plot(sweed_p1$mid, sweed_p1$norm_clr,xlab = "",ylab = "CLR",col = "white",xlim = c((gene_chrom$Start[1] -500000),(gene_chrom$End[1] +500000)))
    lines(smooth_sweed_p1, lwd = 2, col = red)
    rect(gene_chrom$Start[1], -1000, gene_chrom$End[1], 1000, border = NA, col = col2alpha(blu,0.2))
    rect(gene_chrom$Start[1], -1, gene_chrom$End[1], 1, border = NA, col = col2alpha(blu,0.2))
    plot(sweed_p2$mid, sweed_p2$norm_clr,xlab = "",ylab = "CLR",col = "white",xlim = c((gene_chrom$Start[1] -500000),(gene_chrom$End[1] +500000)))
    lines(smooth_sweed_p2, lwd = 2, col = blu)
    rect(gene_chrom$Start[1], -1000, gene_chrom$End[1], 1000, border = NA, col = col2alpha(blu,0.2))
    rect(gene_chrom$Start[1], -1, gene_chrom$End[1], 1, border = NA, col = col2alpha(blu,0.2))
    
  }
  
}
dev.off()

cool_genes_comps <- read.table("D:/test_plots/pmp_fst_dxy_gemma.txt", header = TRUE, stringsAsFactors = FALSE)
i <- 1
pdf(paste("D:/test_plots/tests_",cool_genes_comps$interesting_because[i],".pdf",sep = ""),width=7,height=5)
for (i in c(1:nrow(cool_genes_comps)))
{
  cool_stage <- cool_genes_comps$stage[i]
  pop_gen_comp <- cool_genes_comps$pops[i]
  cool_p1 <- cool_genes_comps$p1[i]
  cool_p2 <- cool_genes_comps$p2[i]
  cool_h <- cool_genes_comps$h[i]
  cool_cols <- c(gre,blu,grb)
  
  p1_inds <- master[which(master$f1 == cool_p1 & master$stage == cool_stage),]
  p2_inds <- master[which(master$f1 == cool_p2 & master$stage == cool_stage),]
  hy_inds <- master[which(master$f1 == cool_h & master$stage == cool_stage),]
  p1_inds <- p1_inds$sample
  p2_inds <- p2_inds$sample
  hy_inds <- hy_inds$sample
  
  p1 <- strsplit(pop_gen_comp, "x")[[1]][1]
  p2 <- strsplit(pop_gen_comp, "x")[[1]][2]
  dxy_main <- read.csv(paste("D:/Martin Lab/rna_2018/fst/dna/",pop_gen_comp,"_popgen_dna_stats_corr_dxy.csv",sep = ""), header = TRUE, stringsAsFactors = FALSE)
  
  fst_main <- read.table(paste("C:/Users/jmcgirr/Documents/all_2018_samples/fst_dna/",pop_gen_comp,"_fst",sep = ""), header = TRUE, stringsAsFactors = FALSE)
  fst_main <- fst_main[which(fst_main$WEIR_AND_COCKERHAM_FST >=0 & fst_main$WEIR_AND_COCKERHAM_FST <=1),]
  lake <- strsplit(pop_gen_comp,split ="")[[1]][1]
  taj_m_main <-    read.table(paste("C:/Users/jmcgirr/Documents/all_2018_samples/fst_dna/",lake,"m_taj_d_20kb.txt.Tajima.D", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  taj_a_main <-    read.table(paste("C:/Users/jmcgirr/Documents/all_2018_samples/fst_dna/",lake,"a_taj_d_20kb.txt.Tajima.D", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  taj_p_main <-    read.table(paste("C:/Users/jmcgirr/Documents/all_2018_samples/fst_dna/",lake,"p_taj_d_20kb.txt.Tajima.D", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  taj_m_main <- na.omit(taj_m_main)
  taj_a_main <- na.omit(taj_a_main)
  taj_p_main <- na.omit(taj_p_main)
  sweed_p1_main <-    read.table(paste("C:/Users/jmcgirr/Documents/all_2018_samples/sweed/correct_seq/SweeD_Report.",p1,"_pop_bottle_58_grid_500_sweeps.txt", sep = ""), header = FALSE, stringsAsFactors = FALSE)
  sweed_p2_main <-    read.table(paste("C:/Users/jmcgirr/Documents/all_2018_samples/sweed/correct_seq/SweeD_Report.",p2,"_pop_bottle_58_grid_500_sweeps.txt", sep = ""), header = FALSE, stringsAsFactors = FALSE)
  
  cool_genes <- strsplit(cool_genes_comps$genes[i],";")[[1]]
  for (cool_gene in cool_genes)
  {
    gene_chrom <- mrna[which(mrna$related_accession == cool_gene),]
    sweed_p1 <- sweed_p1_main[which(sweed_p1_main$V1 ==gene_chrom$Chr[1]),]
    sweed_p2 <- sweed_p2_main[which(sweed_p2_main$V1 ==gene_chrom$Chr[1]),]
    fst <- fst_main[which(fst_main$CHROM == gene_chrom$Chr[1]),]
    fixed <- fst[which(fst$WEIR_AND_COCKERHAM_FST == 1),]
    smooth_fst <-   smooth.spline(fst$POS, fst$WEIR_AND_COCKERHAM_FST, spar = .2)
    dxy <- dxy_main[which(dxy_main$scaffold == gene_chrom$Chr[1]),]
    
    smooth_dxy <-   smooth.spline(dxy$start, dxy$corr_dxy, spar = .2)
    smooth_p1_pi <- smooth.spline(dxy$start, dxy[,14], spar = .2)
    smooth_p2_pi <- smooth.spline(dxy$start, dxy[,16], spar = .2)
    
    tj_m <- taj_m_main[which(taj_m_main$CHROM == gene_chrom$Chr[1]),]
    smooth_m <-   smooth.spline(tj_m$BIN_START, tj_m$TajimaD, spar = .1)
    tj_p <- taj_p_main[which(taj_p_main$CHROM == gene_chrom$Chr[1]),]
    smooth_p <-   smooth.spline(tj_p$BIN_START, tj_p$TajimaD, spar = .1)
    tj_a <- taj_a_main[which(taj_a_main$CHROM == gene_chrom$Chr[1]),]
    smooth_a <-   smooth.spline(tj_a$BIN_START, tj_a$TajimaD, spar = .1)
    
    
    plot(fixed$POS, fixed$WEIR_AND_COCKERHAM_FST,ylim = c(0,1),
         xlab = "", ylab = "fst", yaxt = 'n',xlim =c((gene_chrom$Start[1] -500000),(gene_chrom$End[1] +500000)), main = paste(cool_p1, "vs", cool_p2,cool_gene, sep = " "))
    axis(2,c(0,0.5,1))
    lines(smooth_fst, lwd = 2, col = "black")
    rect(gene_chrom$Start[1], -1000, gene_chrom$End[1], 1000, border = NA, col = col2alpha(blu,0.2))
    plot(dxy$start, dxy$corr_dxy,col = "white",
         xlab = "", ylab = "dxy",xlim =c((gene_chrom$Start[1] -500000),(gene_chrom$End[1] +500000)))#, main = paste(cool_p1, "vs", cool_p2, sep = " "))
    lines(smooth_dxy, lwd = 2, col = "grey")
    rect(gene_chrom$Start[1], -1000, gene_chrom$End[1], 1000, border = NA, col = col2alpha(blu,0.2))
    plot(dxy$start, dxy[,14],col = "white",
         xlab = "", ylab = "pi",xlim =c((gene_chrom$Start[1] -500000),(gene_chrom$End[1] +500000)))#, main = paste(cool_p1, "vs", cool_p2, sep = " "))
    lines(smooth_p1_pi, lwd = 2, col = red)
    lines(smooth_p2_pi, lwd = 2, col = blu)
    rect(gene_chrom$Start[1], -1000, gene_chrom$End[1], 1000, border = NA, col = col2alpha(blu,0.2))
    
    plot(tj_p$BIN_START, tj_p$TajimaD,ylab = "D",ylim = c(-2.5,2.5),xlab = "",
         xlim = c((gene_chrom$Start[1] -500000),(gene_chrom$End[1] +500000)),col = "white", yaxt = "n")
    axis(2,c(-2,-1,0,1,2))
    lines(smooth_m, lwd = 2, col = gre)
    lines(smooth_p, lwd = 2, col = blu)
    lines(smooth_a, lwd = 2, col = red)
    rect(gene_chrom$Start[1], -1000, gene_chrom$End[1], 1000, border = NA, col = col2alpha(blu,0.2))
    
    sweed_p1$mid <- sweed_p1$V2 + ((sweed_p1$V2[2]-sweed_p1$V2[1])/2)
    sweed_p2$mid <- sweed_p2$V2 + ((sweed_p2$V2[2]-sweed_p2$V2[1])/2)
    sweed_p1$norm_clr = (sweed_p1$V3-min(sweed_p1$V3))/(max(sweed_p1$V3)-min(sweed_p1$V3))
    sweed_p2$norm_clr = (sweed_p2$V3-min(sweed_p2$V3))/(max(sweed_p2$V3)-min(sweed_p2$V3))
    smooth_sweed_p1 <-   smooth.spline(sweed_p1$mid, sweed_p1$norm_clr, spar = 0.2)
    smooth_sweed_p2 <-   smooth.spline(sweed_p2$mid, sweed_p2$norm_clr, spar = 0.2)
    
    plot(sweed_p1$mid, sweed_p1$norm_clr,xlab = "",ylab = "CLR",col = "white",xlim = c((gene_chrom$Start[1] -500000),(gene_chrom$End[1] +500000)))
    lines(smooth_sweed_p1, lwd = 2, col = red)
    rect(gene_chrom$Start[1], -1000, gene_chrom$End[1], 1000, border = NA, col = col2alpha(blu,0.2))
    rect(gene_chrom$Start[1], -1, gene_chrom$End[1], 1, border = NA, col = col2alpha(blu,0.2))
    plot(sweed_p2$mid, sweed_p2$norm_clr,xlab = "",ylab = "CLR",col = "white",xlim = c((gene_chrom$Start[1] -500000),(gene_chrom$End[1] +500000)))
    lines(smooth_sweed_p2, lwd = 2, col = blu)
    rect(gene_chrom$Start[1], -1000, gene_chrom$End[1], 1000, border = NA, col = col2alpha(blu,0.2))
    rect(gene_chrom$Start[1], -1, gene_chrom$End[1], 1, border = NA, col = col2alpha(blu,0.2))
    
  }
  
}
dev.off()

cool_genes_comps <- read.table("D:/test_plots/pmp_fst_dxy.txt", header = TRUE, stringsAsFactors = FALSE)
i <- 1
pdf(paste("D:/test_plots/tests_",cool_genes_comps$interesting_because[i],".pdf",sep = ""),width=7,height=5)
for (i in c(1:nrow(cool_genes_comps)))
{
  cool_stage <- cool_genes_comps$stage[i]
  pop_gen_comp <- cool_genes_comps$pops[i]
  cool_p1 <- cool_genes_comps$p1[i]
  cool_p2 <- cool_genes_comps$p2[i]
  cool_h <- cool_genes_comps$h[i]
  cool_cols <- c(gre,blu,grb)
  
  p1_inds <- master[which(master$f1 == cool_p1 & master$stage == cool_stage),]
  p2_inds <- master[which(master$f1 == cool_p2 & master$stage == cool_stage),]
  hy_inds <- master[which(master$f1 == cool_h & master$stage == cool_stage),]
  p1_inds <- p1_inds$sample
  p2_inds <- p2_inds$sample
  hy_inds <- hy_inds$sample
  
  p1 <- strsplit(pop_gen_comp, "x")[[1]][1]
  p2 <- strsplit(pop_gen_comp, "x")[[1]][2]
  dxy_main <- read.csv(paste("D:/Martin Lab/rna_2018/fst/dna/",pop_gen_comp,"_popgen_dna_stats_corr_dxy.csv",sep = ""), header = TRUE, stringsAsFactors = FALSE)
  
  fst_main <- read.table(paste("C:/Users/jmcgirr/Documents/all_2018_samples/fst_dna/",pop_gen_comp,"_fst",sep = ""), header = TRUE, stringsAsFactors = FALSE)
  fst_main <- fst_main[which(fst_main$WEIR_AND_COCKERHAM_FST >=0 & fst_main$WEIR_AND_COCKERHAM_FST <=1),]
  lake <- strsplit(pop_gen_comp,split ="")[[1]][1]
  taj_m_main <-    read.table(paste("C:/Users/jmcgirr/Documents/all_2018_samples/fst_dna/",lake,"m_taj_d_20kb.txt.Tajima.D", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  taj_a_main <-    read.table(paste("C:/Users/jmcgirr/Documents/all_2018_samples/fst_dna/",lake,"a_taj_d_20kb.txt.Tajima.D", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  taj_p_main <-    read.table(paste("C:/Users/jmcgirr/Documents/all_2018_samples/fst_dna/",lake,"p_taj_d_20kb.txt.Tajima.D", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  taj_m_main <- na.omit(taj_m_main)
  taj_a_main <- na.omit(taj_a_main)
  taj_p_main <- na.omit(taj_p_main)
  sweed_p1_main <-    read.table(paste("C:/Users/jmcgirr/Documents/all_2018_samples/sweed/correct_seq/SweeD_Report.",p1,"_pop_bottle_58_grid_500_sweeps.txt", sep = ""), header = FALSE, stringsAsFactors = FALSE)
  sweed_p2_main <-    read.table(paste("C:/Users/jmcgirr/Documents/all_2018_samples/sweed/correct_seq/SweeD_Report.",p2,"_pop_bottle_58_grid_500_sweeps.txt", sep = ""), header = FALSE, stringsAsFactors = FALSE)
  
  cool_genes <- strsplit(cool_genes_comps$genes[i],";")[[1]]
  for (cool_gene in cool_genes)
  {
    gene_chrom <- mrna[which(mrna$related_accession == cool_gene),]
    sweed_p1 <- sweed_p1_main[which(sweed_p1_main$V1 ==gene_chrom$Chr[1]),]
    sweed_p2 <- sweed_p2_main[which(sweed_p2_main$V1 ==gene_chrom$Chr[1]),]
    fst <- fst_main[which(fst_main$CHROM == gene_chrom$Chr[1]),]
    fixed <- fst[which(fst$WEIR_AND_COCKERHAM_FST == 1),]
    smooth_fst <-   smooth.spline(fst$POS, fst$WEIR_AND_COCKERHAM_FST, spar = .2)
    dxy <- dxy_main[which(dxy_main$scaffold == gene_chrom$Chr[1]),]
    
    smooth_dxy <-   smooth.spline(dxy$start, dxy$corr_dxy, spar = .2)
    smooth_p1_pi <- smooth.spline(dxy$start, dxy[,14], spar = .2)
    smooth_p2_pi <- smooth.spline(dxy$start, dxy[,16], spar = .2)
    
    tj_m <- taj_m_main[which(taj_m_main$CHROM == gene_chrom$Chr[1]),]
    smooth_m <-   smooth.spline(tj_m$BIN_START, tj_m$TajimaD, spar = .1)
    tj_p <- taj_p_main[which(taj_p_main$CHROM == gene_chrom$Chr[1]),]
    smooth_p <-   smooth.spline(tj_p$BIN_START, tj_p$TajimaD, spar = .1)
    tj_a <- taj_a_main[which(taj_a_main$CHROM == gene_chrom$Chr[1]),]
    smooth_a <-   smooth.spline(tj_a$BIN_START, tj_a$TajimaD, spar = .1)
    
    
    plot(fixed$POS, fixed$WEIR_AND_COCKERHAM_FST,ylim = c(0,1),
         xlab = "", ylab = "fst", yaxt = 'n',xlim =c((gene_chrom$Start[1] -500000),(gene_chrom$End[1] +500000)), main = paste(cool_p1, "vs", cool_p2,cool_gene, sep = " "))
    axis(2,c(0,0.5,1))
    lines(smooth_fst, lwd = 2, col = "black")
    rect(gene_chrom$Start[1], -1000, gene_chrom$End[1], 1000, border = NA, col = col2alpha(blu,0.2))
    plot(dxy$start, dxy$corr_dxy,col = "white",
         xlab = "", ylab = "dxy",xlim =c((gene_chrom$Start[1] -500000),(gene_chrom$End[1] +500000)))#, main = paste(cool_p1, "vs", cool_p2, sep = " "))
    lines(smooth_dxy, lwd = 2, col = "grey")
    rect(gene_chrom$Start[1], -1000, gene_chrom$End[1], 1000, border = NA, col = col2alpha(blu,0.2))
    plot(dxy$start, dxy[,14],col = "white",
         xlab = "", ylab = "pi",xlim =c((gene_chrom$Start[1] -500000),(gene_chrom$End[1] +500000)))#, main = paste(cool_p1, "vs", cool_p2, sep = " "))
    lines(smooth_p1_pi, lwd = 2, col = red)
    lines(smooth_p2_pi, lwd = 2, col = blu)
    rect(gene_chrom$Start[1], -1000, gene_chrom$End[1], 1000, border = NA, col = col2alpha(blu,0.2))
    
    plot(tj_p$BIN_START, tj_p$TajimaD,ylab = "D",ylim = c(-2.5,2.5),xlab = "",
         xlim = c((gene_chrom$Start[1] -500000),(gene_chrom$End[1] +500000)),col = "white", yaxt = "n")
    axis(2,c(-2,-1,0,1,2))
    lines(smooth_m, lwd = 2, col = gre)
    lines(smooth_p, lwd = 2, col = blu)
    lines(smooth_a, lwd = 2, col = red)
    rect(gene_chrom$Start[1], -1000, gene_chrom$End[1], 1000, border = NA, col = col2alpha(blu,0.2))
    
    sweed_p1$mid <- sweed_p1$V2 + ((sweed_p1$V2[2]-sweed_p1$V2[1])/2)
    sweed_p2$mid <- sweed_p2$V2 + ((sweed_p2$V2[2]-sweed_p2$V2[1])/2)
    sweed_p1$norm_clr = (sweed_p1$V3-min(sweed_p1$V3))/(max(sweed_p1$V3)-min(sweed_p1$V3))
    sweed_p2$norm_clr = (sweed_p2$V3-min(sweed_p2$V3))/(max(sweed_p2$V3)-min(sweed_p2$V3))
    smooth_sweed_p1 <-   smooth.spline(sweed_p1$mid, sweed_p1$norm_clr, spar = 0.2)
    smooth_sweed_p2 <-   smooth.spline(sweed_p2$mid, sweed_p2$norm_clr, spar = 0.2)
    
    plot(sweed_p1$mid, sweed_p1$norm_clr,xlab = "",ylab = "CLR",col = "white",xlim = c((gene_chrom$Start[1] -500000),(gene_chrom$End[1] +500000)))
    lines(smooth_sweed_p1, lwd = 2, col = red)
    rect(gene_chrom$Start[1], -1000, gene_chrom$End[1], 1000, border = NA, col = col2alpha(blu,0.2))
    rect(gene_chrom$Start[1], -1, gene_chrom$End[1], 1, border = NA, col = col2alpha(blu,0.2))
    plot(sweed_p2$mid, sweed_p2$norm_clr,xlab = "",ylab = "CLR",col = "white",xlim = c((gene_chrom$Start[1] -500000),(gene_chrom$End[1] +500000)))
    lines(smooth_sweed_p2, lwd = 2, col = blu)
    rect(gene_chrom$Start[1], -1000, gene_chrom$End[1], 1000, border = NA, col = col2alpha(blu,0.2))
    rect(gene_chrom$Start[1], -1, gene_chrom$End[1], 1, border = NA, col = col2alpha(blu,0.2))
    
  }
  
}
dev.off()










cool_genes_out_dir <- "C:/Users/jmcgirr/Documents/all_2018_samples/manuscript_figs/fig_5_can_genes/"

cool_stage <- "8dpf"
cool_gene <- "XM_015371309.1"
pop_gen_comp <- "caxcp"
p1_col <- red
p2_col <- blu

#plot counts
{
  #cool_p1 <- "CRPA"
  #cool_p2 <- "CRPM"
  #cool_h <- "CRPP"
  #p3 <- "CAxCM"
  #p4 <- "CAxCP"
  #p5 <- "CMxCP"
  #p6 <- "UPxUA"
  #p7 <- "NCA"
  #cool_stage <- "8dpf"
  #cool_cols <- c(red,gre,blu,yel,pur,grb,"white","black")
  
  cool_p1 <- "OSPA"
  cool_p2 <- "OSPP"
  cool_h <- "OAxOP"
  p3 <- "OSPM"
  p4 <- "OAxOM"
  p5 <- "OPxOM"
  
  cool_p1 <- "CRPA"
  cool_p2 <- "CRPP"
  cool_h <- "CAxCP"
  p3 <- "CRPM"
  p4 <- "CAxCM"
  p5 <- "CMxCP"
  
  p6 <- "UPxUA"
  p7 <- "NCA"
  
  cool_cols <- c(red,blu,pur,gre,yel,grb,"white","black")
  
  p1_inds <- master[which(master$f1 == cool_p1 & master$stage == cool_stage),]
  p2_inds <- master[which(master$f1 == cool_p2 & master$stage == cool_stage),]
  hy_inds <- master[which(master$f1 == cool_h & master$stage == cool_stage),]
  p3_inds <- master[which(master$f1 == p3 & master$stage == cool_stage),]
  p4_inds <- master[which(master$f1 == p4 & master$stage == cool_stage),]
  p5_inds <- master[which(master$f1 == p5 & master$stage == cool_stage),]
  p6_inds <- master[which(master$f1 == p6 & master$stage == cool_stage),]
  p7_inds <- master[which(master$f1 == p7 & master$stage == cool_stage),]
  p1_inds <- p1_inds$sample
  p2_inds <- p2_inds$sample
  hy_inds <- hy_inds$sample
  p3_inds <-p3_inds$sample
  p4_inds <-p4_inds$sample
  p5_inds <-p5_inds$sample
  p6_inds <-p6_inds$sample
  p7_inds <- p7_inds$sample
  norm_cts_gene <- norm_cts[which(norm_cts$related_accession == cool_gene),]
  ph_cts <- norm_cts_gene[c(p1_inds,p2_inds,hy_inds,p3_inds,p4_inds,p5_inds,p6_inds,p7_inds)]
  ph_cts  <- as.numeric(ph_cts[1,])
  inds <- c(p1_inds,p2_inds,hy_inds,p3_inds,p4_inds,p5_inds,p6_inds,p7_inds)
  type_inds <- c(rep(cool_p1,length(p1_inds)),rep(cool_p2,length(p2_inds)),rep(cool_h,length(hy_inds)),
                 rep(p3,length(p3_inds)),rep(p4,length(p4_inds)),rep(p5,length(p5_inds)),rep(p6,length(p6_inds)),rep(p7,length(p7_inds)))
  plot_counts <- data.frame(type_ind = type_inds, sample = inds, cts = ph_cts)
  plot_counts$type_ind <- as.character(plot_counts$type_ind)
  plot_counts$type_ind <- factor(plot_counts$type_ind, levels=unique(plot_counts$type_ind))
  plot_counts$type_ind <- factor(plot_counts$type_ind, level = c(cool_p1,cool_h,cool_p2, p3,p4,p5,p6,p7))
  plot_title <- paste(cool_gene,"\n",norm_cts_gene$gene_name[1], sep = "")
  cool_gene_out <- paste(cool_genes_out_dir,norm_cts_gene$gene_name[1],"_",cool_h,"_",cool_stage, sep = "")
  
  #tiff(paste(cool_gene_out, "_gene_counts.tiff", sep = ""), width = 6.2, height = 5, units = 'in', res = 1000)
  plot(plot_counts$type_ind,plot_counts$cts, ylab = "normalized counts", border = "white",
       cex.axis = 1.1, cex.names = 1.5,cex.lab = 1.1)#, main = plot_title)
  p1_cts <- plot_counts[which(plot_counts$type_ind == cool_p1),]
  points(jitter(rep(1,nrow(p1_cts)),10), p1_cts$cts, pch = 21,bg =cool_cols[1],  col = "black", cex = 2)
  p2_cts <- plot_counts[which(plot_counts$type_ind == cool_p2),]
  points(jitter(rep(3,nrow(p2_cts)),2), p2_cts$cts, pch = 21,bg =cool_cols[2],  col = "black", cex = 2)
  hy_cts <- plot_counts[which(plot_counts$type_ind == cool_h),]
  points(jitter(rep(2,nrow(hy_cts)),9), hy_cts$cts, pch = 21,bg =cool_cols[3],  col = "black", cex = 2)
  
  p3_cts <- plot_counts[which(plot_counts$type_ind == p3),]
  points(jitter(rep(4,nrow(p3_cts)),3), p3_cts$cts, pch = 21,bg =cool_cols[4],  col = "black", cex = 2)
  p4_cts <- plot_counts[which(plot_counts$type_ind == p4),]
  points(jitter(rep(5,nrow(p4_cts)),3), p4_cts$cts, pch = 21,bg =cool_cols[5],  col = "black", cex = 2)
  p5_cts <- plot_counts[which(plot_counts$type_ind == p5),]
  points(jitter(rep(6,nrow(p5_cts)),2), p5_cts$cts, pch = 21,bg =cool_cols[6],  col = "black", cex = 2)
  p6_cts <- plot_counts[which(plot_counts$type_ind == p6),]
  points(jitter(rep(7,nrow(p6_cts)),3), p6_cts$cts, pch = 21,bg =cool_cols[7],  col = "black", cex = 2)
  p7_cts <- plot_counts[which(plot_counts$type_ind == p7),]
  points(jitter(rep(8,nrow(p7_cts)),2), p7_cts$cts, pch = 21,bg =cool_cols[8],  col = "black", cex = 2)
  #dev.off()
}

#load pop comp
{

p1_inds <- master[which(master$f1 == cool_p1 & master$stage == cool_stage),]
p2_inds <- master[which(master$f1 == cool_p2 & master$stage == cool_stage),]
hy_inds <- master[which(master$f1 == cool_h & master$stage == cool_stage),]
p1_inds <- p1_inds$sample
p2_inds <- p2_inds$sample
hy_inds <- hy_inds$sample

p1 <- strsplit(pop_gen_comp, "x")[[1]][1]
p2 <- strsplit(pop_gen_comp, "x")[[1]][2]
dxy_main <- read.csv(paste("D:/Martin Lab/rna_2018/fst/dna/",pop_gen_comp,"_popgen_dna_stats_corr_dxy.csv",sep = ""), header = TRUE, stringsAsFactors = FALSE)

fst_main <- read.table(paste("C:/Users/jmcgirr/Documents/all_2018_samples/fst_dna/",pop_gen_comp,"_fst",sep = ""), header = TRUE, stringsAsFactors = FALSE)
fst_main <- fst_main[which(fst_main$WEIR_AND_COCKERHAM_FST >=0 & fst_main$WEIR_AND_COCKERHAM_FST <=1),]
lake <- strsplit(pop_gen_comp,split ="")[[1]][1]
taj_m_main <-    read.table(paste("C:/Users/jmcgirr/Documents/all_2018_samples/fst_dna/",lake,"m_taj_d_20kb.txt.Tajima.D", sep = ""), header = TRUE, stringsAsFactors = FALSE)
taj_a_main <-    read.table(paste("C:/Users/jmcgirr/Documents/all_2018_samples/fst_dna/",lake,"a_taj_d_20kb.txt.Tajima.D", sep = ""), header = TRUE, stringsAsFactors = FALSE)
taj_p_main <-    read.table(paste("C:/Users/jmcgirr/Documents/all_2018_samples/fst_dna/",lake,"p_taj_d_20kb.txt.Tajima.D", sep = ""), header = TRUE, stringsAsFactors = FALSE)
taj_m_main <- na.omit(taj_m_main)
taj_a_main <- na.omit(taj_a_main)
taj_p_main <- na.omit(taj_p_main)
sweed_p1_main <-    read.table(paste("C:/Users/jmcgirr/Documents/all_2018_samples/sweed/",p1,"_pop_bottle_58_sweeps.txt", sep = ""), header = FALSE, stringsAsFactors = FALSE)
sweed_p2_main <-    read.table(paste("C:/Users/jmcgirr/Documents/all_2018_samples/sweed/",p2,"_pop_bottle_58_sweeps.txt", sep = ""), header = FALSE, stringsAsFactors = FALSE)
}

xrange <- c(450000,800000)
par(mfrow = c(7,1),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
#plot
{
gene_chrom <- mrna[which(mrna$related_accession == cool_gene),]
sweed_p1 <- sweed_p1_main[which(sweed_p1_main$V1 ==gene_chrom$Chr[1]),]
sweed_p2 <- sweed_p2_main[which(sweed_p2_main$V1 ==gene_chrom$Chr[1]),]
fst <- fst_main[which(fst_main$CHROM == gene_chrom$Chr[1]),]
fixed <- fst[which(fst$WEIR_AND_COCKERHAM_FST == 1),]
smooth_fst <-   smooth.spline(fst$POS, fst$WEIR_AND_COCKERHAM_FST, spar = .5)
dxy <- dxy_main[which(dxy_main$scaffold == gene_chrom$Chr[1]),]

smooth_dxy <-   smooth.spline(dxy$start, dxy$corr_dxy, spar = .5)
smooth_p1_pi <- smooth.spline(dxy$start, dxy[,14], spar = .5)
smooth_p2_pi <- smooth.spline(dxy$start, dxy[,16], spar = .5)

tj_m <- taj_m_main[which(taj_m_main$CHROM == gene_chrom$Chr[1]),]
smooth_m <-   smooth.spline(tj_m$BIN_START, tj_m$TajimaD, spar = .5)
tj_p <- taj_p_main[which(taj_p_main$CHROM == gene_chrom$Chr[1]),]
smooth_p <-   smooth.spline(tj_p$BIN_START, tj_p$TajimaD, spar = .5)
tj_a <- taj_a_main[which(taj_a_main$CHROM == gene_chrom$Chr[1]),]
smooth_a <-   smooth.spline(tj_a$BIN_START, tj_a$TajimaD, spar = .5)


plot(fixed$POS, fixed$WEIR_AND_COCKERHAM_FST,ylim = c(0,1),
     xlab = "", ylab = "fst", yaxt = 'n',xlim =xrange, xaxt = "n")
axis(2,c(0,0.5,1),cex.axis=1.1)
lines(smooth_fst, lwd = 2, col = "black")
rect(gene_chrom$Start[1], -1000, gene_chrom$End[1], 1000, border = NA, col = col2alpha(blu,0.2))
plot(dxy$start, dxy$corr_dxy,col = "white",
     xlab = "", ylab = "dxy",xlim =xrange, xaxt = "n", yaxt="n")
axis(2, at=c(0,0.005,0.01, 0.015), labels=c(0,"",0.01,""),cex.axis=1.1)
lines(smooth_dxy, lwd = 2, col = "grey")
rect(gene_chrom$Start[1], -1000, gene_chrom$End[1], 1000, border = NA, col = col2alpha(blu,0.2))

plot(dxy$start, dxy[,14],col = "white",
     xlab = "", ylab = "pi",xlim =xrange, xaxt = "n", yaxt="n")
axis(2, at=c(0,0.005,0.01, 0.015), labels=c(0,"",0.01,""),cex.axis=1.1)
lines(smooth_p1_pi, lwd = 2, col = p1_col)
lines(smooth_p2_pi, lwd = 2, col = p2_col)
rect(gene_chrom$Start[1], -1000, gene_chrom$End[1], 1000, border = NA, col = col2alpha(blu,0.2))

plot(tj_p$BIN_START, tj_p$TajimaD,ylab = "D",ylim = c(-1.5,2.5),xlab = "",
     xlim = xrange,col = "white", yaxt = "n", xaxt = "n")
axis(2,c(-1,0,1,2),cex.axis=1.1)
lines(smooth_m, lwd = 2, col = gre)
lines(smooth_p, lwd = 2, col = blu)
lines(smooth_a, lwd = 2, col = red)
rect(gene_chrom$Start[1], -1000, gene_chrom$End[1], 1000, border = NA, col = col2alpha(blu,0.2))

smooth_sweed_p1 <-   smooth.spline(sweed_p1$V2, sweed_p1$V4, spar = 0.3)
smooth_sweed_p2 <-   smooth.spline(sweed_p2$V2, sweed_p2$V4, spar = 0.3)
plot(sweed_p1$V2, sweed_p1$V4,xlab = "",ylab = "CLR",xlim = xrange, xaxt = "n", yaxt="n")
axis(2, at=c(0,0.015,0.03), labels=c(0,0.015,0.03),cex.axis=1.1)
lines(smooth_sweed_p1, lwd = 2, col = p1_col)
rect(gene_chrom$Start[1], -1000, gene_chrom$End[1], 1000, border = NA, col = col2alpha(blu,0.2))
rect(gene_chrom$Start[1], -1, gene_chrom$End[1], 1, border = NA, col = col2alpha(blu,0.2))
plot(sweed_p2$V2, sweed_p2$V4,xlab = "",ylab = "CLR",xlim = xrange, xaxt = "n", yaxt="n")
axis(2, at=c(0,0.03,0.06), labels=c(0,0.03,0.06),cex.axis=1.1)
lines(smooth_sweed_p2, lwd = 2, col = p2_col)
rect(gene_chrom$Start[1], -1000, gene_chrom$End[1], 1000, border = NA, col = col2alpha(blu,0.2))
rect(gene_chrom$Start[1], -1, gene_chrom$End[1], 1, border = NA, col = col2alpha(blu,0.2))
  
gem <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/GEMMA/output_58/mean_gammas_up_jaw_summed_20kb.txt", header = FALSE, stringsAsFactors = FALSE)
head(gem)
gem <- cbind(gem, colsplit(gem[,3], "\\|", c("end", "pip")))
gem <- gem[which(gem$V1 ==gene_chrom$Chr[1]),]
smooth_gem <-   smooth.spline(gem$V2, gem$pip, spar = .1)
#tiff(paste(cool_gene_out, "_gemma.tiff", sep = ""), width = 5, height = 3, units = 'in', res = 1000)
plot(gem$V2, gem$pip,xlab = "",ylab = "PIP",xlim = xrange, 
     xaxt = 'n', col = "white", yaxt = "n")
lines(smooth_gem, lwd = 2, col = pur)
axis(2, at=c(0,0.001,0.002), labels=c(0,"",0.002),cex.axis=1.1)
#axis(2,c(0,0.002))
rect(gene_chrom$Start[1], -1000, gene_chrom$End[1], 1000, border = NA, col = col2alpha(blu,0.2))
#dev.off()
#plot(gem$V2, gem$pip,xlab = "",ylab = "PIP", 
#      col = "white")#, yaxt = "n")
#lines(smooth_gem, lwd = 2, col = pur)

}

















cool_stage <- "8dpf"
cool_gene <- "XM_015395742.1"
pop_gen_comp <- "cmxcp"
#plot counts
{
  #cool_p1 <- "CRPA"
  #cool_p2 <- "CRPM"
  #cool_h <- "CRPP"
  #p3 <- "CAxCM"
  #p4 <- "CAxCP"
  #p5 <- "CMxCP"
  #p6 <- "UPxUA"
  #p7 <- "NCA"
  #cool_stage <- "8dpf"
  #cool_cols <- c(red,gre,blu,yel,pur,grb,"white","black")
  
  cool_p1 <- "OSPA"
  cool_p2 <- "OSPP"
  cool_h <- "OAxOP"
  p3 <- "OSPM"
  p4 <- "OAxOM"
  p5 <- "OPxOM"

  cool_p1 <- "CRPA"
  cool_p2 <- "CRPP"
  cool_h <- "CAxCP"
  p3 <- "CRPM"
  p4 <- "CAxCM"
  p5 <- "CMxCP"
  
  p6 <- "UPxUA"
  p7 <- "NCA"

  cool_cols <- c(red,blu,pur,gre,yel,grb,"white","black")
  
  p1_inds <- master[which(master$f1 == cool_p1 & master$stage == cool_stage),]
  p2_inds <- master[which(master$f1 == cool_p2 & master$stage == cool_stage),]
  hy_inds <- master[which(master$f1 == cool_h & master$stage == cool_stage),]
  p3_inds <- master[which(master$f1 == p3 & master$stage == cool_stage),]
  p4_inds <- master[which(master$f1 == p4 & master$stage == cool_stage),]
  p5_inds <- master[which(master$f1 == p5 & master$stage == cool_stage),]
  p6_inds <- master[which(master$f1 == p6 & master$stage == cool_stage),]
  p7_inds <- master[which(master$f1 == p7 & master$stage == cool_stage),]
  p1_inds <- p1_inds$sample
  p2_inds <- p2_inds$sample
  hy_inds <- hy_inds$sample
  p3_inds <-p3_inds$sample
  p4_inds <-p4_inds$sample
  p5_inds <-p5_inds$sample
  p6_inds <-p6_inds$sample
  p7_inds <- p7_inds$sample
  norm_cts_gene <- norm_cts[which(norm_cts$related_accession == cool_gene),]
  ph_cts <- norm_cts_gene[c(p1_inds,p2_inds,hy_inds,p3_inds,p4_inds,p5_inds,p6_inds,p7_inds)]
  ph_cts  <- as.numeric(ph_cts[1,])
  inds <- c(p1_inds,p2_inds,hy_inds,p3_inds,p4_inds,p5_inds,p6_inds,p7_inds)
  type_inds <- c(rep(cool_p1,length(p1_inds)),rep(cool_p2,length(p2_inds)),rep(cool_h,length(hy_inds)),
                 rep(p3,length(p3_inds)),rep(p4,length(p4_inds)),rep(p5,length(p5_inds)),rep(p6,length(p6_inds)),rep(p7,length(p7_inds)))
  plot_counts <- data.frame(type_ind = type_inds, sample = inds, cts = ph_cts)
  plot_counts$type_ind <- as.character(plot_counts$type_ind)
  plot_counts$type_ind <- factor(plot_counts$type_ind, levels=unique(plot_counts$type_ind))
  plot_counts$type_ind <- factor(plot_counts$type_ind, level = c(cool_p1,cool_h,cool_p2, p3,p4,p5,p6,p7))
  plot_title <- paste(cool_gene,"\n",norm_cts_gene$gene_name[1], sep = "")
  cool_gene_out <- paste(cool_genes_out_dir,norm_cts_gene$gene_name[1],"_",cool_h,"_",cool_stage, sep = "")
  
  #tiff(paste(cool_gene_out, "_gene_counts.tiff", sep = ""), width = 6.2, height = 5, units = 'in', res = 1000)
  plot(plot_counts$type_ind,plot_counts$cts, ylab = "normalized counts", border = "white",
       cex.axis = 1.1, cex.names = 1.5,cex.lab = 1.1)#, main = plot_title)
  p1_cts <- plot_counts[which(plot_counts$type_ind == cool_p1),]
  points(jitter(rep(1,nrow(p1_cts)),10), p1_cts$cts, pch = 21,bg =cool_cols[1],  col = "black", cex = 2)
  p2_cts <- plot_counts[which(plot_counts$type_ind == cool_p2),]
  points(jitter(rep(3,nrow(p2_cts)),2), p2_cts$cts, pch = 21,bg =cool_cols[2],  col = "black", cex = 2)
  hy_cts <- plot_counts[which(plot_counts$type_ind == cool_h),]
  points(jitter(rep(2,nrow(hy_cts)),9), hy_cts$cts, pch = 21,bg =cool_cols[3],  col = "black", cex = 2)
  
  p3_cts <- plot_counts[which(plot_counts$type_ind == p3),]
  points(jitter(rep(4,nrow(p3_cts)),3), p3_cts$cts, pch = 21,bg =cool_cols[4],  col = "black", cex = 2)
  p4_cts <- plot_counts[which(plot_counts$type_ind == p4),]
  points(jitter(rep(5,nrow(p4_cts)),3), p4_cts$cts, pch = 21,bg =cool_cols[5],  col = "black", cex = 2)
  p5_cts <- plot_counts[which(plot_counts$type_ind == p5),]
  points(jitter(rep(6,nrow(p5_cts)),2), p5_cts$cts, pch = 21,bg =cool_cols[6],  col = "black", cex = 2)
  p6_cts <- plot_counts[which(plot_counts$type_ind == p6),]
  points(jitter(rep(7,nrow(p6_cts)),3), p6_cts$cts, pch = 21,bg =cool_cols[7],  col = "black", cex = 2)
  p7_cts <- plot_counts[which(plot_counts$type_ind == p7),]
  points(jitter(rep(8,nrow(p7_cts)),2), p7_cts$cts, pch = 21,bg =cool_cols[8],  col = "black", cex = 2)
  #dev.off()
}

cool_p1 <- "CRPM"
cool_p2 <- "CRPP"
cool_h <- "CMxCP"
cool_cols <- c(gre,blu,grb)
#plot pop gen
{
  p1_inds <- master[which(master$f1 == cool_p1 & master$stage == cool_stage),]
  p2_inds <- master[which(master$f1 == cool_p2 & master$stage == cool_stage),]
  hy_inds <- master[which(master$f1 == cool_h & master$stage == cool_stage),]
  p1_inds <- p1_inds$sample
  p2_inds <- p2_inds$sample
  hy_inds <- hy_inds$sample
  
  p1 <- strsplit(pop_gen_comp, "x")[[1]][1]
  p2 <- strsplit(pop_gen_comp, "x")[[1]][2]
  dxy <- read.csv(paste("D:/Martin Lab/rna_2018/fst/dna/",pop_gen_comp,"_popgen_dna_stats_corr_dxy.csv",sep = ""), header = TRUE, stringsAsFactors = FALSE)
  #pi_p1 <- read.table(paste("C:/Users/jmcgirr/Documents/all_2018_samples/fst_dna/",p1,"_pi_20kb.txt.windowed.pi",sep = ""), header = TRUE, stringsAsFactors = FALSE)
  #pi_p2 <- read.table(paste("C:/Users/jmcgirr/Documents/all_2018_samples/fst_dna/",p2,"_pi_20kb.txt.windowed.pi",sep = ""), header = TRUE, stringsAsFactors = FALSE)

  fst <- read.table(paste("C:/Users/jmcgirr/Documents/all_2018_samples/fst_dna/",pop_gen_comp,"_fst",sep = ""), header = TRUE, stringsAsFactors = FALSE)
  #pop_gen_dna[pop_gen_dna<0] <- 0
  #pop_gen_dna <-  na.omit(pop_gen_dna)
  lake <- strsplit(pop_gen_comp,split ="")[[1]][1]
  
  taj_m <-    read.table(paste("C:/Users/jmcgirr/Documents/all_2018_samples/fst_dna/",lake,"m_taj_d_20kb.txt.Tajima.D", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  taj_a <-    read.table(paste("C:/Users/jmcgirr/Documents/all_2018_samples/fst_dna/",lake,"a_taj_d_20kb.txt.Tajima.D", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  taj_p <-    read.table(paste("C:/Users/jmcgirr/Documents/all_2018_samples/fst_dna/",lake,"p_taj_d_20kb.txt.Tajima.D", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  taj_m <- na.omit(taj_m)
  taj_a <- na.omit(taj_a)
  taj_p <- na.omit(taj_p)
  sweed_p1 <-    read.table(paste("C:/Users/jmcgirr/Documents/all_2018_samples/sweed/",p1,"_pop_bottle_58_sweeps.txt", sep = ""), header = FALSE, stringsAsFactors = FALSE)
  sweed_p2 <-    read.table(paste("C:/Users/jmcgirr/Documents/all_2018_samples/sweed/",p2,"_pop_bottle_58_sweeps.txt", sep = ""), header = FALSE, stringsAsFactors = FALSE)
  
  gene_chrom <- mrna[which(mrna$related_accession == cool_gene),]
  sweed_p1 <- sweed_p1[which(sweed_p1$V1 ==gene_chrom$Chr[1]),]
  sweed_p2 <- sweed_p2[which(sweed_p2$V1 ==gene_chrom$Chr[1]),]
  fst <- fst[which(fst$CHROM == gene_chrom$Chr[1]),]
  fixed <- fst[which(fst$WEIR_AND_COCKERHAM_FST == 1),]
  smooth_fst <-   smooth.spline(fst$POS, fst$WEIR_AND_COCKERHAM_FST, spar = .3)
  dxy <- dxy[which(dxy$scaffold == gene_chrom$Chr[1]),]
  smooth_dxy <-   smooth.spline(dxy$start, dxy$corr_dxy, spar = .01)
  #smooth_p1_pi <- smooth.spline(pi_p1$BIN_START, pi_p1$PI, spar = .1)
  #smooth_p2_pi <- smooth.spline(pi_p2$BIN_START, pi_p2$PI, spar = .1)
  smooth_p1_pi <- smooth.spline(dxy$start, dxy[,14], spar = .1)
  smooth_p2_pi <- smooth.spline(dxy$start, dxy[,16], spar = .1)
  
  tj_m <- taj_m[which(taj_m$CHROM == gene_chrom$Chr[1]),]
  smooth_m <-   smooth.spline(tj_m$BIN_START, tj_m$TajimaD, spar = .2)
  tj_p <- taj_p[which(taj_p$CHROM == gene_chrom$Chr[1]),]
  smooth_p <-   smooth.spline(tj_p$BIN_START, tj_p$TajimaD, spar = .2)
  tj_a <- taj_a[which(taj_a$CHROM == gene_chrom$Chr[1]),]
  smooth_a <-   smooth.spline(tj_a$BIN_START, tj_a$TajimaD, spar = .2)
  
  #tiff(paste(cool_gene_out, "_fst_dna.tiff", sep = ""), width = 5, height = 3, units = 'in', res = 1000)
  plot(fixed$POS, fixed$WEIR_AND_COCKERHAM_FST,ylim = c(0,1),
       xaxt = 'n', xlab = "", ylab = "", yaxt = 'n',xlim =c((gene_chrom$Start[1] -100000),(gene_chrom$End[1] +100000)))#, main = paste(cool_p1, "vs", cool_p2, sep = " "))
  axis(2,c(0,0.5,1))
  lines(smooth_fst, lwd = 2, col = "black")
  rect(gene_chrom$Start[1], -1000, gene_chrom$End[1], 1000, border = NA, col = col2alpha(blu,0.2))
  #dev.off()
  #tiff(paste(cool_gene_out, "_dxy_dna.tiff", sep = ""), width = 5, height = 3, units = 'in', res = 1000)
  plot(dxy$start, dxy$corr_dxy,col = "white",
       xaxt = 'n', xlab = "", ylab = "",xlim =c((gene_chrom$Start[1] -100000),(gene_chrom$End[1] +100000)))#, main = paste(cool_p1, "vs", cool_p2, sep = " "))
  #axis(2,c(0.004,0.005))
  lines(smooth_dxy, lwd = 2, col = "grey")
  rect(gene_chrom$Start[1], -1000, gene_chrom$End[1], 1000, border = NA, col = col2alpha(blu,0.2))
  #dev.off()
  #tiff(paste(cool_gene_out, "_pi_dna.tiff", sep = ""), width = 5, height = 3, units = 'in', res = 1000)
  plot(dxy$start, dxy[,14],col = "white",xaxt = "n",
        xlab = "", ylab = "",xlim =c((gene_chrom$Start[1] -100000),(gene_chrom$End[1] +100000)))#, main = paste(cool_p1, "vs", cool_p2, sep = " "))
  #axis(2,c(0.001,0.005))
  lines(smooth_p1_pi, lwd = 2, col = red)
  lines(smooth_p2_pi, lwd = 2, col = blu)
  rect(gene_chrom$Start[1], -1000, gene_chrom$End[1], 1000, border = NA, col = col2alpha(blu,0.2))
  #dev.off()
  
  #tiff(paste(cool_gene_out, "_tajima_dna.tiff", sep = ""), width = 5, height = 3, units = 'in', res = 1000)
  plot(tj_p$BIN_START, tj_p$TajimaD,xlab = round((gene_chrom$Start[1] +100000)-(gene_chrom$End[1] -100000),digits = -3),ylab = "D",ylim = c(-1,2.5),
       xlim = c((gene_chrom$Start[1] -100000),(gene_chrom$End[1] +100000)), xaxt = 'n',col = "white", yaxt = "n")
  axis(2,c(-1,0,1,2))
  lines(smooth_m, lwd = 2, col = gre)
  lines(smooth_p, lwd = 2, col = blu)
  lines(smooth_a, lwd = 2, col = red)
  rect(gene_chrom$Start[1], -1000, gene_chrom$End[1], 1000, border = NA, col = col2alpha(blu,0.2))
  #dev.off()
  
  smooth_sweed_p1 <-   smooth.spline(sweed_p1$V2, sweed_p1$V4, spar = .3)
  smooth_sweed_p2 <-   smooth.spline(sweed_p2$V2, sweed_p1$V4, spar = .3)
  #tiff(paste(cool_gene_out, "_sweed.tiff", sep = ""), width = 5, height = 3.5, units = 'in', res = 1000)
  plot(sweed_p1$V2, sweed_p1$V4,xlab = "",ylab = "CLR",xlim = c((gene_chrom$Start[1] -100000),(gene_chrom$End[1] +100000)), xaxt = 'n',col = "white")
  lines(smooth_sweed_p1, lwd = 2, col = red)
  lines(smooth_sweed_p2, lwd = 2, col = blu)
  rect(gene_chrom$Start[1], -1000, gene_chrom$End[1], 1000, border = NA, col = col2alpha(blu,0.2))
  
  #dev.off()
}

#plot gemma
{
gem <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/GEMMA/output_58/mean_gammas_up_jaw_summed_20kb.txt", header = FALSE, stringsAsFactors = FALSE)
head(gem)
gem <- cbind(gem, colsplit(gem[,3], "\\|", c("end", "pip")))
gem <- gem[which(gem$V1 ==gene_chrom$Chr[1]),]
smooth_gem <-   smooth.spline(gem$V2, gem$pip, spar = .1)
#tiff(paste(cool_gene_out, "_gemma.tiff", sep = ""), width = 5, height = 3, units = 'in', res = 1000)
plot(gem$V2, gem$pip,xlab = "",ylab = "PIP",xlim = c((gene_chrom$Start[1] -100000),(gene_chrom$End[1] +100000)), 
     xaxt = 'n', col = "white")#, yaxt = "n")
lines(smooth_gem, lwd = 2, col = pur)
#axis(2,c(0,0.002))
rect(gene_chrom$Start[1], -1000, gene_chrom$End[1], 1000, border = NA, col = col2alpha(blu,0.2))
#dev.off()
}
# plot ase
{
  ind <- hy_inds[[1]]
  cts <- read.csv(paste(maternal_counts_dir,ind, "_parental_counts_features.txt",sep = ""),na.strings=c("","NA"), header = TRUE, stringsAsFactors = FALSE, sep = "\t") 
  cts <- cbind(cts, colsplit(cts$mrnaID, ";", c("related_accession", "gene_name")))
  cts <- cts[which(cts$chrom == gene_chrom$Chr[1]),]# & cts$related_accession == cool_gene),]
  keeps <- c("chrom","mrnaID","position")
  allele_cts_ps <- cts[keeps]
  exons <- gff[grep(cool_gene, gff$V9),]
  exons <- exons[which(exons$V3 == "exon"),]
  exons <- exons[c("V1","V4","V5")]
  for (ind in hy_inds)
  {
    cts <- read.csv(paste(maternal_counts_dir,ind, "_parental_counts_features.txt",sep = ""),na.strings=c("","NA"), header = TRUE, stringsAsFactors = FALSE, sep = "\t") 
    cts <- cbind(cts, colsplit(cts$mrnaID, ";", c("related_accession", "gene_name")))
    cts <- cts[which(cts$chrom == gene_chrom$Chr[1] & cts$related_accession == cool_gene),]
    colnames(cts)[colnames(cts)=="momCount"] <- paste(ind,"_momCount", sep = "")
    colnames(cts)[colnames(cts)=="dadCount"] <- paste(ind,"_dadCount", sep = "")
    keeps <- c("chrom","mrnaID","position",paste(ind,"_momCount", sep = ""),paste(ind,"_dadCount", sep = ""))
    allele_cts <- merge(allele_cts_ps, cts[keeps],all = TRUE, by = c("chrom","mrnaID","position"))
    hy_cts <- na.omit(as.vector(as.matrix(allele_cts[,c(4,5)])))
    mid_point <- max(hy_cts) - ((max(hy_cts) - min(hy_cts))/2)
    y_ax_range <- c(min(hy_cts) -((max(hy_cts)-mid_point)*0.2),max(hy_cts) +((max(hy_cts)-mid_point)*0.1))
    
    #tiff(paste(cool_gene_out,"_",ind, "_allele_counts.tiff", sep = ""), width = 5.5, height = 3.5, units = 'in', res = 1000)
    plot(allele_cts$position, c(1:length(allele_cts$position)),ylim = y_ax_range,axes = FALSE,
         xlim = c(min(exons$V4)-300,max(exons$V5)+300),col = "white",xlab = "" ,ylab = "", xaxt = 'n')#, main = plot_title)
    axis(2, cex.axis = 1.05)
    axis(1, labels = FALSE, lwd.tick=0)
    title(xlab = paste(round(((max(exons$V5)+300)-(min(exons$V4)+300)), digits = -3)," bp", sep = ""), line=0.3, cex.lab=1)
    abline(h =mid_point,col = col2alpha(blu,0.2), lwd = 2)
    for(i in c(1:nrow(exons)))
    {
      rect(exons$V4[i], mid_point - ((y_ax_range[2]-y_ax_range[1]) /40), exons$V5[i], mid_point + ((y_ax_range[2]-y_ax_range[1])/40), border = NA, col = col2alpha("darkblue",0.2))
    }
    points(allele_cts$position, allele_cts[,4], pch = 21,bg =cool_cols[1],  col = "black", cex = 1.1)
    points(allele_cts$position, allele_cts[,5], pch = 21,bg =cool_cols[2],  col = "black", cex = 1.1)
    #dev.off()
  }
  
}



  
  mpp1 <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/ase_data/mpp1_8dpf_avg_caxcp.txt", header= TRUE, stringsAsFactors = FALSE)
  head(mpp1)
  
  #tiff(paste(cool_gene_out,"_",ind, "_allele_counts.tiff", sep = ""), width = 5.5, height = 4.5, units = 'in', res = 1000)
  
  plot(mpp1$position, mpp1$dad_avg,ylim = c(min(mpp1$mom_avg)-30,max(mpp1$dad_avg)),axes = FALSE,
       xlim = c(min(exons$V4)-300,max(exons$V5)+300),col = "white",xlab = "" ,ylab = "", xaxt = 'n')#, main = plot_title)
  axis(2, cex.axis = 1.05)
  #axis(1, labels = FALSE, lwd.tick=0)
  points(mpp1$position, mpp1$mom_avg, pch = 21,bg =cool_cols[1],  col = "black", cex = 1.3)
  points(mpp1$position, mpp1$dad_avg, pch = 21,bg =cool_cols[2],  col = "black", cex = 1.3)
  abline(h =130,col = col2alpha(blu,0.2), lwd = 2)
  
  for(i in c(1:nrow(exons)))
  {
    rect(exons$V4[i], 120, exons$V5[i], 140, border = NA, col = col2alpha(blu,0.2))
  }
  #dev.off()
  



mpp1 <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/ase_data/mpp1_8dpf_avg_caxcp.txt", header= TRUE, stringsAsFactors = FALSE)
head(mpp1)

#tiff(paste(cool_gene_out,"_",ind, "_allele_counts.tiff", sep = ""), width = 7, height = 6, units = 'in', res = 1000)

plot(mpp1$position, mpp1$dad_avg,ylim = c(min(mpp1$mom_avg)-30,max(mpp1$dad_avg)),axes = FALSE,
     xlim = c(min(exons$V4)-300,max(exons$V5)+300),col = "white",xlab = "" ,ylab = "", xaxt = 'n')#, main = plot_title)
axis(2, cex.axis = 1.05)
#axis(1, labels = FALSE, lwd.tick=0)
points(mpp1$position, mpp1$mom_avg, pch = 21,bg =cool_cols[1],  col = "black", cex = 1.8)
points(mpp1$position, mpp1$dad_avg, pch = 21,bg =cool_cols[2],  col = "black", cex = 1.8)
abline(h =125,col = col2alpha(blu,0.2), lwd = 2)

for(i in c(1:nrow(exons)))
{
  rect(exons$V4[i], 115, exons$V5[i], 130, border = NA, col = col2alpha(blu,0.2))
}
#dev.off()



#####



