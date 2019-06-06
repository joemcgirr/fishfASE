# source("http://bioconductor.org/biocLite.R")
# biocLite("DESeq2")
library("DESeq2")
library(reshape2)
library(seqinr)
library(plyr)
library(MASS)
library(AER)
#RnaSeqSampleSize for power

#### counts ####
#setwd("D:/Martin Lab/RNA-seq/axm/post_reviews/")
#dpf17_parental <- read.table("D:/Martin Lab/RNA-seq/axm/post_reviews/17dpf_parental_am_counts_my_gtf_geneid_final", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
#dpf17_axm <- read.table("D:/Martin Lab/RNA-seq/axm/post_reviews/17dpf_axm_counts_my_gtf_geneid_final", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
#dpf8 <- read.table("D:/Martin Lab/RNA-seq/mapping_correction/DESeq_counts_am_embryo.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
#round3_and_4 <- read.table("D:/Martin Lab/rna_2018/all_2018_samples/conditions/all_48hpf_and_8dpf_counts_my_gtf_geneid_final", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
#head(dpf8)
#head(round3_and_4)
#dpf8$Chr <- NULL
#dpf8$Start <- NULL
#dpf8$End <- NULL
#dpf8$Strand <- NULL
#dpf8$Length <- NULL
#dpf17_parental$Chr <- NULL
#dpf17_parental$Start <- NULL
#dpf17_parental$End <- NULL
#dpf17_parental$Strand <- NULL
#dpf17_parental$Length <- NULL
#dpf17_axm$Chr <- NULL
#dpf17_axm$Start <- NULL
#dpf17_axm$End <- NULL
#dpf17_axm$Strand <- NULL
#dpf17_axm$Length <- NULL
#round3_and_4$Chr <- NULL
#round3_and_4$Start <- NULL
#round3_and_4$End <- NULL
#round3_and_4$Strand <- NULL
#round3_and_4$Length <- NULL
#head(dpf8)
#head(dpf17_axm)
#df1 <- merge(dpf17_axm, dpf17_parental, by = ("Geneid"))
#df2 <- merge(round3_and_4, df1, by = ("Geneid"))
#df3 <- merge(dpf8, df2, by = ("Geneid"))
#
#head(df3)
#length(names(df3))
#df3 = df3[-1,]
#write.table(df3, "D:/Martin Lab/RNA-seq/axm/post_reviews/all_rna_4_rounds_counts.txt", row.names = FALSE, quote = FALSE, sep = "\t")

all_comps <- read.table("D:/Martin Lab/RNA-seq/axm/post_reviews/comparisons_mirror_design.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
head(all_comps)
final_features <- read.table("D:/Martin Lab/rna_2018/all_2018_samples/features_gff.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t", quote = "")
all_cts <- read.table("D:/Martin Lab/RNA-seq/axm/post_reviews/all_rna_4_rounds_counts.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
head(all_cts)
blast_key <- read.table("D:/Martin Lab/rna_2018/all_2018_samples/blast/cyprinodon_to_danio_two_way_best_hit_symbols.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
head(blast_key)
#i <- 15
#subset_comps <- c(1,2,4,15)

for (i in (1:nrow(all_comps)))
  #for (i in subset_comps)
{
  pops1_name <- all_comps$pops1[i]
  pops2_name <- all_comps$pops2[i]
  pops1 <- strsplit(pops1_name,"_&_")[[1]]
  pops2 <- strsplit(pops2_name,"_&_")[[1]]
  stage <- all_comps$stage[i]
  pops1_name <- paste("(", pops1_name, ")", sep = "")
  pops2_name <- paste("(", pops2_name, ")", sep = "")
  comp_name <- paste("condition_species_",pops1_name, "_vs_", pops2_name, "_" ,stage,".txt", sep = "")
  single_run <- paste(pops1_name, "_vs_", pops2_name, "_" ,stage, sep = "")
  
  setwd("D:/Martin Lab/RNA-seq/axm/post_reviews/conditions")
  master <- read.table("D:/Martin Lab/RNA-seq/axm/post_reviews/table_maker_mirror_design.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
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
  write.table(comp_table, comp_name, row.names = FALSE, quote= FALSE,sep="\t")
  
  comp_file <- single_run
  
  setwd("D:/Martin Lab/RNA-seq/axm/post_reviews/conditions")
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
  write.table(cts, cts_data, row.names = FALSE, quote= FALSE,sep="\t") 
  
  ### DIFFERENTIAL EXPRESSION ANALYSIS ###
  #vignette("DESeq2")
  cts <- as.matrix(read.table(cts_data ,sep = "\t",header = TRUE,row.names=1))
  head(cts)
  nrow(cts)
  colData <- as.matrix(read.table(col_data ,header = TRUE,row.names=1))
  head(colData)
  ncol(cts)
  nrow(colData)
  
  
  #if ((add_sequencing_round_to_model > 3) & unique(samples_a$sequencing_round) == unique(samples_b$sequencing_round))
  #{
  #  dds <- DESeqDataSetFromMatrix(countData = cts,
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
  # most stringent -- filter out genes if read count is less than 10 for any individual.
  
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
  de_total_plot <- paste(de_total, "DE", sep = " ")
  prop_de <- paste((100*(round(prop_de, digits = 3))), "% DE", sep = "")
  plotMA(resLFC, ylim=c(-6,5), main = comp_file)
  legend("bottomright", legend=c(sample_sizes_plot,total_genes_plot, de_total_plot, prop_de),cex=1.0, bty = 'n')
  legend("bottomleft", legend=c(paste((100*(round(de_dn, digits = 3))), "% DE down", sep = "")),cex=0.8, bty = 'n')
  legend("topleft", legend=c(paste((100*(round(de_up, digits = 3))), "% DE up", sep = "")),cex=0.8, bty = 'n')
  
  #head(resLFC)
  de_plot <- paste("D:/Martin Lab/RNA-seq/axm/post_reviews/de_plots/",comp_file, "_de_plot.tiff", sep = "")
  tiff(de_plot, width = 7, height = 6, units = 'in', res = 1000)
  plotMA(resLFC, ylim=c(-6,5), main = comp_file)
  legend("bottomright", legend=c(sample_sizes_plot,total_genes_plot, de_total_plot, prop_de),cex=1.0, bty = 'n')
  legend("bottomleft", legend=c(paste((100*(round(de_dn, digits = 3))), "% DE down", sep = "")),cex=0.8, bty = 'n')
  legend("topleft", legend=c(paste((100*(round(de_up, digits = 3))), "% DE up", sep = "")),cex=0.8, bty = 'n')
  dev.off()
  
  #### overlap with genes ###
  
  write.csv(res_ordered, file= cts_data_genes)
  genes <- read.csv(cts_data_genes, header = TRUE, stringsAsFactors = FALSE)
  final <- merge(genes, final_features, by = c("Geneid"))
  length(unique(final$Geneid))
  length(unique(final$Geneid)) / length(unique(res_ordered$Geneid))
  head(final)
  nrow(final)
  head(blast_key)
  final <- merge(final, blast_key,all.x = TRUE, by = c("product_accession"))
  final <- final[order(final$padj, decreasing = FALSE),]
  write.csv(final, genes_out, row.names = FALSE)
  
}






#####
########################################
########################################
##### inheritance plots ################
########################################

#### make plot like McManus 2010 4B

inheritance_comps <- read.table("D:/Martin Lab/RNA-seq/axm/post_reviews/inheritance_patterns/inheritance_comparisons.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
DE_genes_dir <- "D:/Martin Lab/RNA-seq/axm/post_reviews/conditions/"
inheritance_plots_dir <- "D:/Martin Lab/RNA-seq/axm/post_reviews/inheritance_patterns/plots/"
inheritance_dir <- "D:/Martin Lab/RNA-seq/axm/post_reviews/inheritance_patterns/"
blast_key <- read.table("D:/Martin Lab/rna_2018/all_2018_samples/blast/cyprinodon_to_danio_two_way_best_hit_symbols.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")

#i <- 3

for (i in (1:nrow(inheritance_comps)))
{
  parent_pops1_name <- inheritance_comps$parent_pops1[i]
  parent_pops2_name <- inheritance_comps$parent_pops2[i]
  hy_name <- inheritance_comps$hybrids[i]
  stage <- inheritance_comps$stage[i]
  parent_pops_v_hy  <- paste(DE_genes_dir, "DE_(", parent_pops1_name,"_&_",parent_pops2_name, ")_vs_(",hy_name, ")_",stage,"_genes",".csv",sep = "")
  parent_pops1_v_hy <- paste(DE_genes_dir, "DE_(", parent_pops1_name, ")_vs_(",hy_name, ")_",stage,"_genes",".csv",sep = "")
  parent_pops2_v_hy <- paste(DE_genes_dir, "DE_(", parent_pops2_name, ")_vs_(",hy_name, ")_",stage,"_genes",".csv",sep = "")
  parent_pops1_v_parent_pops2 <- paste(DE_genes_dir, "DE_(", parent_pops1_name, ")_vs_(",parent_pops2_name, ")_",stage,"_genes",".csv",sep = "")
  outfile_plot  <- paste(inheritance_plots_dir, hy_name, "_",stage,"_misexpression",".tiff",sep = "")
  plot_title  <- paste(hy_name, " ",stage," misexpression",sep = "")
  
  ph    <- read.csv(parent_pops_v_hy  , header = TRUE, stringsAsFactors = FALSE)
  p1h   <- read.csv(parent_pops1_v_hy , header = TRUE, stringsAsFactors = FALSE)
  p2h   <- read.csv(parent_pops2_v_hy , header = TRUE, stringsAsFactors = FALSE)
  p1p2  <- read.csv(parent_pops1_v_parent_pops2, header = TRUE, stringsAsFactors = FALSE)
  
  head(ph)
  
  ph$lfc_ph     <- ph$log2FoldChange
  p1h$lfc_p1h   <- p1h$log2FoldChange 
  p2h$lfc_p2h   <- p2h$log2FoldChange  
  p1p2$lfc_p1p2 <- p1p2$log2FoldChange
  ph$padj_ph     <- ph$padj
  p1h$padj_p1h   <- p1h$padj 
  p2h$padj_p2h   <- p2h$padj  
  p1p2$padj_p1p2 <- p1p2$padj
  
  head(ph)
  
  keeps <- c("tag","lfc_ph","padj_ph","lfc_p1h","padj_p1h")
  all_genes <- merge(ph, p1h, by = ("tag"))
  all_genes <- all_genes[keeps]
  all_genes <- merge(all_genes, p2h, by = ("tag"))
  keeps <- c("tag","lfc_ph","padj_ph","lfc_p1h","padj_p1h","lfc_p2h", "padj_p2h")
  all_genes <- all_genes[keeps]
  all_genes <- merge(all_genes, p1p2, by = ("tag"))
  keeps <- c("tag","lfc_ph","padj_ph","lfc_p1h","padj_p1h","lfc_p2h", "padj_p2h", "lfc_p1p2", "padj_p1p2", "symbol", "related_accession", "product_accession", "name")
  all_genes <- all_genes[keeps]
  head(all_genes)
  nrow(all_genes)
  length(unique(all_genes$tag))
  
  
  total_genes <- nrow(all_genes)
  con <- all_genes[which(all_genes$padj_ph > 0.05 & all_genes$padj_p1h > 0.05 & all_genes$padj_p2h > 0.05 & all_genes$padj_p1p2 > 0.05),]
  con$inheritance <- "conserved"
  con$inheritance <- 1
  nrow(con)
  add <- all_genes[which(all_genes$padj_p1h > 0.05 & all_genes$padj_p2h > 0.05 & all_genes$padj_p1p2 < 0.05),]
  nrow(add)
  add1 <- add[which(add$lfc_p1h < 0 & add$lfc_p2h > 0),]
  nrow(add1)
  add2 <- add[which(add$lfc_p1h > 0 & add$lfc_p2h < 0),]
  add <- rbind(add1, add2)
  add$inheritance <- "additive"
  add$inheritance <- 2
  nrow(add)
  p1_dom <- all_genes[which(all_genes$padj_p1h > 0.05 & all_genes$padj_p2h < 0.05 & all_genes$padj_p1p2 < 0.05),]
  p1_dom$inheritance <- "p1_dominant"
  p1_dom$inheritance <- 3
  nrow(p1_dom)
  p2_dom <- all_genes[which(all_genes$padj_p1h < 0.05 & all_genes$padj_p2h > 0.05 & all_genes$padj_p1p2 < 0.05),]
  p2_dom$inheritance <- "p2_dominant"
  p2_dom$inheritance <- 4
  nrow(p2_dom)
  over_dom <- all_genes[which(all_genes$padj_ph < 0.05 & all_genes$lfc_p1h > 0 & all_genes$lfc_p2h > 0),]
  over_dom$inheritance <- "over_dominant"
  over_dom$inheritance <- 5
  nrow(over_dom)
  under_dom <- all_genes[which(all_genes$padj_ph < 0.05 & all_genes$lfc_p1h < 0 & all_genes$lfc_p2h < 0),]
  under_dom$inheritance <- "under_dominant"
  under_dom$inheritance <- 6
  nrow(under_dom)
  
  prop_con      <-paste(round((nrow(con)/total_genes * 100), digits = 2),"%", sep = "") 
  prop_add      <-paste(round((nrow(add)/total_genes * 100), digits = 2),"%", sep = "")
  prop_p1_dom   <-paste(round((nrow(p1_dom)/total_genes * 100), digits = 2),"%", sep = "")      
  prop_p2_dom   <-paste(round((nrow(p2_dom)/total_genes * 100), digits = 2),"%", sep = "")      
  prop_overdom  <-paste(round((nrow(over_dom)/total_genes * 100), digits = 2),"%", sep = "")         
  prop_underdom <-paste(round((nrow(under_dom)/total_genes * 100), digits = 2),"%", sep = "")      
  
  inheritance <- rbind(con, add, p1_dom, p2_dom, over_dom, under_dom)
  head(inheritance)
  nrow(inheritance)
  length(unique(inheritance$tag))
  inheritance <- merge(inheritance, blast_key,all.x = TRUE, by = c("product_accession"))
  inheritance_table_out  <- paste(inheritance_dir, hy_name, "_",stage,"_inheritance.txt",sep = "")
  #write.table(inheritance, inheritance_table_out, row.names = FALSE, quote= FALSE,sep="\t") 
  
  
  
  cols <- 	c("#f6c700", "#E77200", "#FF00CC","#000000", "#bd1e24", "#0067a7")
  
  #con = yellow ="#f6c700"
  #add = orange = "#E77200"
  #p1_dom = pink = "#FF00CC"
  #p2_dom = black = "#000000"
  #over_dom = red = #bd1e24
  #under_dom = blue = #0067a7
  
  inheritance <- inheritance[order(inheritance$inheritance),]
  cols_in <- cols[inheritance$inheritance]
  x_name <- paste(parent_pops1_name, hy_name, sep = " vs. ")
  y_name <- paste(parent_pops2_name, hy_name, sep = " vs. ")
  
  tiff(outfile_plot, width = 6, height = 6, units = 'in', res = 1000)
  plot(inheritance$lfc_p1h, inheritance$lfc_p2h, pch =16, col = cols_in, 
       ylab = "log2 fold change molluscivores vs. hybrids",
       xlab = "log2 fold change generalists vs. hybrids",cex.axis=1.5,cex.lab = 1.5,
       main = "" )
  abline(v=0, col = "black", lty = 3, lwd = 1.8)
  abline(h=0, col = "black", lty = 3, lwd = 1.8)
  #legend("topleft", inset=0,
  #       c(paste("conserved", prop_con),
  #         paste("additive", prop_add),
  #         paste("generalist dominant", prop_p1_dom),
  #         paste("molluscivore dominant", prop_p2_dom), 
  #         paste("overdominant", prop_overdom), 
  #         paste("underdominant", prop_underdom)),
  #       cex = 1.1, fill=cols, horiz=FALSE, bty="n")
  dev.off()

  cols_b <- c("#f6c700", "#0067a7", "#bd1e24","#000000", "#FF00CC","#E77200")
  labs <- c("additive","generalist dominant","molluscivore dominant","overdominant","underdominant", "conserved")
  bars <- c((nrow(con)/total_genes * 100), (nrow(under_dom)/total_genes * 100), (nrow(over_dom)/total_genes * 100), (nrow(p2_dom)/total_genes * 100), (nrow(p1_dom)/total_genes * 100),(nrow(add)/total_genes * 100))
  
  #tiff(paste(outfile_plot, "bar.tiff", sep = ""), width = 6, height = 4, units = 'in', res = 1000)
  #b1 <- barplot(bars, las =1, horiz = TRUE, font.axis = 2, cex.axis=1.5, xlim = c(0,55), col = cols_b)
  #barplot(bars, las =1, horiz = TRUE, font.axis = 2, cex.axis=1.5, xlim = c(0,55), col = cols_b)
  #
  #offset <- 4
  #text(y=b1[6,], x=(nrow(add)/total_genes * 100)+offset,prop_add, font = 2, cex = 1)
  #text(y=b1[5,], x=(nrow(p1_dom)/total_genes * 100)+offset,prop_p1_dom, font = 2, cex = 1)
  #text(y=b1[4,], x=(nrow(p2_dom)/total_genes * 100)+offset,prop_p2_dom, font = 2, cex = 1)
  #text(y=b1[3,], x=(nrow(over_dom)/total_genes * 100)+offset,"25.83 %", font = 2, cex = 1)
  #text(y=b1[2,], x=(nrow(under_dom)/total_genes * 100)+offset,"25.77 %", font = 2, cex = 1)
  #text(y=b1[1,], x=(nrow(con)/total_genes * 100)+offset,prop_con, font = 2, cex = 1)
  #dev.off()
  
  tiff(paste(outfile_plot, "bar.tiff", sep = ""), width = 6, height = 4, units = 'in', res = 1000)
  b1 <- barplot(bars, las =1, horiz = TRUE, font.axis = 2, cex.axis=1.5, xlim = c(0,100), col = cols_b)
  barplot(bars, las =1, horiz = TRUE, font.axis = 2, cex.axis=1.5, xlim = c(0,100), col = cols_b)
  
  offset <- 6
  text(y=b1[6,], x=(nrow(add)/total_genes * 100)+offset,prop_add, font = 2, cex = 1)
  text(y=b1[5,], x=(nrow(p1_dom)/total_genes * 100)+offset,prop_p1_dom, font = 2, cex = 1)
  text(y=b1[4,], x=(nrow(p2_dom)/total_genes * 100)+offset,prop_p2_dom, font = 2, cex = 1)
  text(y=b1[3,], x=(nrow(over_dom)/total_genes * 100)+offset,prop_overdom, font = 2, cex = 1)
  text(y=b1[2,], x=(nrow(under_dom)/total_genes * 100)+offset,prop_underdom, font = 2, cex = 1)
  text(y=b1[1,], x=(nrow(con)/total_genes * 100)+offset,prop_con, font = 2, cex = 1)
  dev.off()
}  

#Hybrid inheritance was considered additive if gene expression was intermediate between
#generalists and molluscivores with significant differential expression between generalists and
#molluscivores. Inheritance was dominant if expression was intermediate between parental
#species and hybrid expression was significantly different from one parent but not the other.
#Genes showing misexpression in hybrids showed transgressive inheritance, where hybrid gene
#expression was significantly higher (overdominant) or lower (underdominant) than parental
#populations.

# schematic
cols = c("#bd1e24","#800080","#0067a7")
par(mfrow= c(1,1))
for (i in collapser)
{

tiff("D:/Martin Lab/RNA-seq/axm/post_reviews/inheritance_patterns/plots/conserved.tiff", width = 5, height = 6, units = 'in', res = 600)
par(lwd=3)
barplot(c(12.5,12.5,12.5), ylim = c(0,29),col = "grey89", border = cols, names.arg = c("G", "H", "M"),xlab="",ylab="",yaxt="n",lwd = 20, cex.names =1.5, font = 2)
axis(2, tcl=0, labels=FALSE, lwd = 3)
par(lwd=1)
#legend(1,28, "conserved", bg = "#f6c700",adj = 0.13, text.font =2, cex = 1.5)
dev.off()

tiff("D:/Martin Lab/RNA-seq/axm/post_reviews/inheritance_patterns/plots/add_1.tiff", width = 5, height = 6, units = 'in', res = 600)
par(lwd=3)
barplot(c(5,12.5,20), ylim = c(0,29),col = "grey89", border = cols, names.arg = c("G", "H", "M"),xlab="",ylab="",yaxt="n",lwd = 20, cex.names =1.5, font = 2)
x <- barplot(c(5,12.5,20), ylim = c(0,29),col = "grey89", border = cols, names.arg = c("G", "H", "M"),xlab="",ylab="",yaxt="n",lwd = 20, cex.names =1.5, font = 2)
axis(2, tcl=0, labels=FALSE)
y <- 22
offset <- 1.2
#legend(1.15,28, "additive", bg = "#E77200",adj = 0.15, text.font =2, cex = 1.5)
lines(x[c(1,3)],c(y, y))
lines(x[c(1,1)],c(y, y-offset))
lines(x[c(3,3)],c(y, y-offset))
text(x[1]+((x[3]-x[1])/2),y+offset,"*", cex = 3, font = 2)
dev.off()

tiff("D:/Martin Lab/RNA-seq/axm/post_reviews/inheritance_patterns/plots/add_2.tiff", width = 5, height = 6, units = 'in', res = 600)
par(lwd=3)
barplot(c(20,12.5,5), ylim = c(0,29),col = "grey89", border = cols, names.arg = c("G", "H", "M"),xlab="",ylab="",yaxt="n",lwd = 20, cex.names =1.5, font = 2)
x <- barplot(c(20,12.5,5), ylim = c(0,29),col = "grey89", border = cols, names.arg = c("G", "H", "M"),xlab="",ylab="",yaxt="n",lwd = 20, cex.names =1.5, font = 2)
axis(2, tcl=0, labels=FALSE)
par(lwd=1)
y <- 22
offset <- 0.5
legend(1.15,28, "additive", bg = "#E77200",adj = 0.15, text.font =2, cex = 1.5)
lines(x[c(1,3)],c(y, y))
lines(x[c(1,1)],c(y, y-offset))
lines(x[c(3,3)],c(y, y-offset))
text(x[1]+((x[3]-x[1])/2),y+offset,"*", cex = 2)
dev.off()

tiff("D:/Martin Lab/RNA-seq/axm/post_reviews/inheritance_patterns/plots/g_dom.tiff", width = 5, height = 6, units = 'in', res = 600)
par(lwd=3)
barplot(c(18,18,5), ylim = c(0,29),col = "grey89", border = cols, names.arg = c("G", "H", "M"),xlab="",ylab="",yaxt="n",lwd = 20, cex.names =1.5, font = 2)
x <- barplot(c(18,18,5), ylim = c(0,29),col = "grey89", border = cols, names.arg = c("G", "H", "M"),xlab="",ylab="",yaxt="n",lwd = 20, cex.names =1.5, font = 2)
axis(2, tcl=0, labels=FALSE)
par(lwd=1)
y <- 22
offset <- 0.5
legend(0.5,28, "generalist dominant", bg = "#FF00CC",adj = 0.08, text.font =2, cex = 1.5)
lines(x[c(1,3)],c(y, y))
lines(x[c(1,1)],c(y, y-offset))
lines(x[c(3,3)],c(y, y-offset))
text(x[1]+((x[3]-x[1])/2),y+offset,"*", cex = 2)
y <- 20
offset <- 0.5
lines(x[c(2,3)],c(y, y))
lines(x[c(2,2)],c(y, y-offset))
lines(x[c(3,3)],c(y, y-offset))
text(x[2]+((x[3]-x[2])/2),y+offset,"*", cex = 2)
dev.off()

tiff("D:/Martin Lab/RNA-seq/axm/post_reviews/inheritance_patterns/plots/g_dom2.tiff", width = 5, height = 6, units = 'in', res = 600)
par(lwd=3)
barplot(c(5,5,18), ylim = c(0,29),col = "grey89", border = cols, names.arg = c("G", "H", "M"),xlab="",ylab="",yaxt="n",lwd = 20, cex.names =1.5, font = 2)
x <- barplot(c(5,5,18), ylim = c(0,29),col = "grey89", border = cols, names.arg = c("G", "H", "M"),xlab="",ylab="",yaxt="n",lwd = 20, cex.names =1.5, font = 2)
axis(2, tcl=0, labels=FALSE)
par(lwd=1)
y <- 22
offset <- 0.5
legend(0.5,28, "generalist dominant", bg = "#FF00CC",adj = 0.08, text.font =2, cex = 1.5)
lines(x[c(1,3)],c(y, y))
lines(x[c(1,1)],c(y, y-offset))
lines(x[c(3,3)],c(y, y-offset))
text(x[1]+((x[3]-x[1])/2),y+offset,"*", cex = 2)
y <- 20
offset <- 0.5
lines(x[c(2,3)],c(y, y))
lines(x[c(2,2)],c(y, y-offset))
lines(x[c(3,3)],c(y, y-offset))
text(x[2]+((x[3]-x[2])/2),y+offset,"*", cex = 2)
dev.off()

tiff("D:/Martin Lab/RNA-seq/axm/post_reviews/inheritance_patterns/plots/m_dom.tiff", width = 5, height = 6, units = 'in', res = 600)
par(lwd=3)
barplot(c(5,18,18), ylim = c(0,29),col = "grey89", border = cols, names.arg = c("G", "H", "M"),xlab="",ylab="",yaxt="n",lwd = 20, cex.names =1.5, font = 2)
x <- barplot(c(5,18,18), ylim = c(0,29),col = "grey89", border = cols, names.arg = c("G", "H", "M"),xlab="",ylab="",yaxt="n",lwd = 20, cex.names =1.5, font = 2)
axis(2, tcl=0, labels=FALSE)
par(lwd=1)
y <- 22
offset <- 0.5
legend(0.35,28, "molluscivore dominant", bg = "#000000",adj = 0.08, text.font =2, text.col = "white", cex = 1.5)
lines(x[c(1,3)],c(y, y))
lines(x[c(1,1)],c(y, y-offset))
lines(x[c(3,3)],c(y, y-offset))
text(x[1]+((x[3]-x[1])/2),y+offset,"*", cex = 2)
y <- 20
offset <- 0.5
lines(x[c(1,2)],c(y, y))
lines(x[c(2,2)],c(y, y-offset))
lines(x[c(1,1)],c(y, y-offset))
text(x[1]+((x[2]-x[1])/2),y+offset,"*", cex = 2)
dev.off()

tiff("D:/Martin Lab/RNA-seq/axm/post_reviews/inheritance_patterns/plots/m_dom2.tiff", width = 5, height = 6, units = 'in', res = 600)
par(lwd=3)
barplot(c(18,5,5), ylim = c(0,29),col = "grey89", border = cols, names.arg = c("G", "H", "M"),xlab="",ylab="",yaxt="n",lwd = 20, cex.names =1.5, font = 2)
x <- barplot(c(18,5,5), ylim = c(0,29),col = "grey89", border = cols, names.arg = c("G", "H", "M"),xlab="",ylab="",yaxt="n",lwd = 20, cex.names =1.5, font = 2)
axis(2, tcl=0, labels=FALSE)
par(lwd=1)
y <- 22
offset <- 0.5
legend(0.35,28, "molluscivore dominant", bg = "#000000",adj = 0.08, text.font =2, text.col = "white", cex = 1.5)
lines(x[c(1,3)],c(y, y))
lines(x[c(1,1)],c(y, y-offset))
lines(x[c(3,3)],c(y, y-offset))
text(x[1]+((x[3]-x[1])/2),y+offset,"*", cex = 2)
y <- 20
offset <- 0.5
lines(x[c(1,2)],c(y, y))
lines(x[c(2,2)],c(y, y-offset))
lines(x[c(1,1)],c(y, y-offset))
text(x[1]+((x[2]-x[1])/2),y+offset,"*", cex = 2)
dev.off()

tiff("D:/Martin Lab/RNA-seq/axm/post_reviews/inheritance_patterns/plots/over.tiff", width = 5, height = 6, units = 'in', res = 600)
par(lwd=3)
barplot(c(5,18,5), ylim = c(0,29),col = "grey89", border = cols, names.arg = c("G", "H", "M"),xlab="",ylab="",yaxt="n",lwd = 20, cex.names =1.5, font = 2)
x <- barplot(c(5,18,5), ylim = c(0,29),col = "grey89", border = cols, names.arg = c("G", "H", "M"),xlab="",ylab="",yaxt="n",lwd = 20, cex.names =1.5, font = 2)
axis(2, tcl=0, labels=FALSE)
par(lwd=1)
y <- 21
offset <- 0.5
legend(0.8,28, "overdominant", bg = "#bd1e24",adj = 0.08, text.font =2, text.col = "white", cex = 1.5)
lines(c(0.7,1.8),c(y, y))
lines(c(1.8,1.8),c(y, y-offset))
lines(x[c(1,1)],c(y, y-offset))
text(x[1]+((1.8-x[1])/2),y+offset,"*", cex = 2)
y <- 21
offset <- 0.5
lines(c(2,3.1),c(y, y))
lines(c(2,2),c(y, y-offset))
lines(c(3.1,3.1),c(y, y-offset))
text(x[2]+((3.1-x[2])/2),y+offset,"*", cex = 2)
dev.off()

tiff("D:/Martin Lab/RNA-seq/axm/post_reviews/inheritance_patterns/plots/under.tiff", width = 5, height = 6, units = 'in', res = 600)
par(lwd=3)
barplot(c(18,5,18), ylim = c(0,29),col = "grey89", border = cols, names.arg = c("G", "H", "M"),xlab="",ylab="",yaxt="n",lwd = 20, cex.names =1.5, font = 2)
x <- barplot(c(18,5,18), ylim = c(0,29),col = "grey89", border = cols, names.arg = c("G", "H", "M"),xlab="",ylab="",yaxt="n",lwd = 20, cex.names =1.5, font = 2)
axis(2, tcl=0, labels=FALSE)
par(lwd=1)
y <- 21
offset <- 0.5
legend(0.75,28, "underdominant", bg = "#0067a7",adj = 0.08, text.font =2, text.col = "white", cex = 1.5)
lines(c(0.7,1.8),c(y, y))
lines(c(1.8,1.8),c(y, y-offset))
lines(x[c(1,1)],c(y, y-offset))
text(x[1]+((1.8-x[1])/2),y+offset,"*", cex = 2)
y <- 21
offset <- 0.5
lines(c(2,3.1),c(y, y))
lines(c(2,2),c(y, y-offset))
lines(c(3.1,3.1),c(y, y-offset))
text(x[2]+((3.1-x[2])/2),y+offset,"*", cex = 2)
dev.off()
}

for (i in collapser)
{
  
  tiff("D:/Martin Lab/RNA-seq/axm/post_reviews/inheritance_patterns/plots/conserved.tiff", width = 4, height = 6.5, units = 'in', res = 600)
  par(lwd=3)
  barplot(c(12.5,12.5,12.5), ylim = c(0,29),col = "grey89", border = cols, names.arg = c("G", "H", "M"),xlab="",ylab="",yaxt="n",lwd = 20, cex.names =1.5, font = 2)
  axis(2, tcl=0, labels=FALSE, lwd = 3)
  
  #legend(1,28, "conserved", bg = "#f6c700",adj = 0.13, text.font =2, cex = 1.5)
  dev.off()
  
  tiff("D:/Martin Lab/RNA-seq/axm/post_reviews/inheritance_patterns/plots/add_1.tiff", width = 4, height = 6.5, units = 'in', res = 600)
  par(lwd=3)
  barplot(c(5,12.5,20), ylim = c(0,29),col = "grey89", border = cols, names.arg = c("G", "H", "M"),xlab="",ylab="",yaxt="n",lwd = 20, cex.names =1.5, font = 2)
  x <- barplot(c(5,12.5,20), ylim = c(0,29),col = "grey89", border = cols, names.arg = c("G", "H", "M"),xlab="",ylab="",yaxt="n",lwd = 20, cex.names =1.5, font = 2)
  axis(2, tcl=0, labels=FALSE, lwd = 3)
  y <- 22
  offset <- 1
  #legend(1.15,28, "additive", bg = "#E77200",adj = 0.15, text.font =2, cex = 1.5)
  lines(x[c(1,3)],c(y, y))
  lines(x[c(1,1)],c(y, y-offset))
  lines(x[c(3,3)],c(y, y-offset))
  text(x[1]+((x[3]-x[1])/2),y+offset,"*", cex = 3, font = 2)
  dev.off()
  
  tiff("D:/Martin Lab/RNA-seq/axm/post_reviews/inheritance_patterns/plots/add_2.tiff", width = 4, height = 6.5, units = 'in', res = 600)
  par(lwd=3)
  barplot(c(20,12.5,5), ylim = c(0,29),col = "grey89", border = cols, names.arg = c("G", "H", "M"),xlab="",ylab="",yaxt="n",lwd = 20, cex.names =1.5, font = 2)
  x <- barplot(c(20,12.5,5), ylim = c(0,29),col = "grey89", border = cols, names.arg = c("G", "H", "M"),xlab="",ylab="",yaxt="n",lwd = 20, cex.names =1.5, font = 2)
  axis(2, tcl=0, labels=FALSE, lwd = 3)
  
  y <- 22
  offset <- 1
  #legend(1.15,28, "additive", bg = "#E77200",adj = 0.15, text.font =2, cex = 1.5)
  lines(x[c(1,3)],c(y, y))
  lines(x[c(1,1)],c(y, y-offset))
  lines(x[c(3,3)],c(y, y-offset))
  text(x[1]+((x[3]-x[1])/2),y+offset,"*", cex = 3, font = 2)
  dev.off()
  
  tiff("D:/Martin Lab/RNA-seq/axm/post_reviews/inheritance_patterns/plots/g_dom.tiff", width = 4, height = 6.5, units = 'in', res = 600)
  par(lwd=3)
  barplot(c(18,18,5), ylim = c(0,29),col = "grey89", border = cols, names.arg = c("G", "H", "M"),xlab="",ylab="",yaxt="n",lwd = 20, cex.names =1.5, font = 2)
  x <- barplot(c(18,18,5), ylim = c(0,29),col = "grey89", border = cols, names.arg = c("G", "H", "M"),xlab="",ylab="",yaxt="n",lwd = 20, cex.names =1.5, font = 2)
  axis(2, tcl=0, labels=FALSE, lwd = 3)
  
  y <- 22
  offset <- 1
  #legend(0.5,28, "generalist dominant", bg = "#FF00CC",adj = 0.08, text.font =2, cex = 1.5)
  lines(x[c(1,3)],c(y, y))
  lines(x[c(1,1)],c(y, y-offset))
  lines(x[c(3,3)],c(y, y-offset))
  text(x[1]+((x[3]-x[1])/2),y+offset,"*", cex = 3, font = 2)
  y <- 20
  offset <- 1
  lines(x[c(2,3)],c(y, y))
  lines(x[c(2,2)],c(y, y-offset))
  lines(x[c(3,3)],c(y, y-offset))
  text(x[2]+((x[3]-x[2])/2),y+offset,"*", cex = 3, font = 2)
  dev.off()
  
  tiff("D:/Martin Lab/RNA-seq/axm/post_reviews/inheritance_patterns/plots/g_dom2.tiff", width = 4, height = 6.5, units = 'in', res = 600)
  par(lwd=3)
  barplot(c(5,5,18), ylim = c(0,29),col = "grey89", border = cols, names.arg = c("G", "H", "M"),xlab="",ylab="",yaxt="n",lwd = 20, cex.names =1.5, font = 2)
  x <- barplot(c(5,5,18), ylim = c(0,29),col = "grey89", border = cols, names.arg = c("G", "H", "M"),xlab="",ylab="",yaxt="n",lwd = 20, cex.names =1.5, font = 2)
  axis(2, tcl=0, labels=FALSE, lwd = 3)
  
  y <- 22
  offset <- 1
  #legend(0.5,28, "generalist dominant", bg = "#FF00CC",adj = 0.08, text.font =2, cex = 1.5)
  lines(x[c(1,3)],c(y, y))
  lines(x[c(1,1)],c(y, y-offset))
  lines(x[c(3,3)],c(y, y-offset))
  text(x[1]+((x[3]-x[1])/2),y+offset,"*", cex = 3, font = 2)
  y <- 20
  offset <- 1
  lines(x[c(2,3)],c(y, y))
  lines(x[c(2,2)],c(y, y-offset))
  lines(x[c(3,3)],c(y, y-offset))
  text(x[2]+((x[3]-x[2])/2),y+offset,"*", cex = 3, font = 2)
  dev.off()
  
  tiff("D:/Martin Lab/RNA-seq/axm/post_reviews/inheritance_patterns/plots/m_dom.tiff", width = 4, height = 6.5, units = 'in', res = 600)
  par(lwd=3)
  barplot(c(5,18,18), ylim = c(0,29),col = "grey89", border = cols, names.arg = c("G", "H", "M"),xlab="",ylab="",yaxt="n",lwd = 20, cex.names =1.5, font = 2)
  x <- barplot(c(5,18,18), ylim = c(0,29),col = "grey89", border = cols, names.arg = c("G", "H", "M"),xlab="",ylab="",yaxt="n",lwd = 20, cex.names =1.5, font = 2)
  axis(2, tcl=0, labels=FALSE, lwd = 3)
  
  y <- 22
  offset <- 1
  #legend(0.35,28, "molluscivore dominant", bg = "#000000",adj = 0.08, text.font =2, text.col = "white", cex = 1.5)
  lines(x[c(1,3)],c(y, y))
  lines(x[c(1,1)],c(y, y-offset))
  lines(x[c(3,3)],c(y, y-offset))
  text(x[1]+((x[3]-x[1])/2),y+offset,"*", cex = 3, font = 2)
  y <- 20
  offset <- 1
  lines(x[c(1,2)],c(y, y))
  lines(x[c(2,2)],c(y, y-offset))
  lines(x[c(1,1)],c(y, y-offset))
  text(x[1]+((x[2]-x[1])/2),y+offset,"*", cex = 3, font = 2)
  dev.off()
  
  tiff("D:/Martin Lab/RNA-seq/axm/post_reviews/inheritance_patterns/plots/m_dom2.tiff", width = 4, height = 6.5, units = 'in', res = 600)
  par(lwd=3)
  barplot(c(18,5,5), ylim = c(0,29),col = "grey89", border = cols, names.arg = c("G", "H", "M"),xlab="",ylab="",yaxt="n",lwd = 20, cex.names =1.5, font = 2)
  x <- barplot(c(18,5,5), ylim = c(0,29),col = "grey89", border = cols, names.arg = c("G", "H", "M"),xlab="",ylab="",yaxt="n",lwd = 20, cex.names =1.5, font = 2)
  axis(2, tcl=0, labels=FALSE, lwd = 3)
  
  y <- 22
  offset <- 1
  #legend(0.35,28, "molluscivore dominant", bg = "#000000",adj = 0.08, text.font =2, text.col = "white", cex = 1.5)
  lines(x[c(1,3)],c(y, y))
  lines(x[c(1,1)],c(y, y-offset))
  lines(x[c(3,3)],c(y, y-offset))
  text(x[1]+((x[3]-x[1])/2),y+offset,"*", cex = 3, font = 2)
  y <- 20
  offset <- 1
  lines(x[c(1,2)],c(y, y))
  lines(x[c(2,2)],c(y, y-offset))
  lines(x[c(1,1)],c(y, y-offset))
  text(x[1]+((x[2]-x[1])/2),y+offset,"*", cex = 3, font = 2)
  dev.off()
  
  tiff("D:/Martin Lab/RNA-seq/axm/post_reviews/inheritance_patterns/plots/over.tiff", width = 4, height = 6.5, units = 'in', res = 600)
  par(lwd=3)
  barplot(c(5,18,5), ylim = c(0,29),col = "grey89", border = cols, names.arg = c("G", "H", "M"),xlab="",ylab="",yaxt="n",lwd = 20, cex.names =1.5, font = 2)
  x <- barplot(c(5,18,5), ylim = c(0,29),col = "grey89", border = cols, names.arg = c("G", "H", "M"),xlab="",ylab="",yaxt="n",lwd = 20, cex.names =1.5, font = 2)
  axis(2, tcl=0, labels=FALSE, lwd = 3)
  
  y <- 21
  offset <- 1
  #legend(0.8,28, "overdominant", bg = "#bd1e24",adj = 0.08, text.font =2, text.col = "white", cex = 1.5)
  lines(c(0.7,1.8),c(y, y))
  lines(c(1.8,1.8),c(y, y-offset))
  lines(x[c(1,1)],c(y, y-offset))
  text(x[1]+((1.8-x[1])/2),y+offset,"*", cex = 3, font = 2)
  y <- 21
  offset <- 1
  lines(c(2,3.1),c(y, y))
  lines(c(2,2),c(y, y-offset))
  lines(c(3.1,3.1),c(y, y-offset))
  text(x[2]+((3.1-x[2])/2),y+offset,"*", cex = 3, font = 2)
  dev.off()
  
  tiff("D:/Martin Lab/RNA-seq/axm/post_reviews/inheritance_patterns/plots/under.tiff", width = 4, height = 6.5, units = 'in', res = 600)
  par(lwd=3)
  barplot(c(18,5,18), ylim = c(0,29),col = "grey89", border = cols, names.arg = c("G", "H", "M"),xlab="",ylab="",yaxt="n",lwd = 20, cex.names =1.5, font = 2)
  x <- barplot(c(18,5,18), ylim = c(0,29),col = "grey89", border = cols, names.arg = c("G", "H", "M"),xlab="",ylab="",yaxt="n",lwd = 20, cex.names =1.5, font = 2)
  axis(2, tcl=0, labels=FALSE, lwd = 3)
  
  y <- 21
  offset <- 1
  #legend(0.75,28, "underdominant", bg = "#0067a7",adj = 0.08, text.font =2, text.col = "white", cex = 1.5)
  lines(c(0.7,1.8),c(y, y))
  lines(c(1.8,1.8),c(y, y-offset))
  lines(x[c(1,1)],c(y, y-offset))
  text(x[1]+((1.8-x[1])/2),y+offset,"*", cex = 3, font = 2)
  y <- 21
  offset <- 1
  lines(c(2,3.1),c(y, y))
  lines(c(2,2),c(y, y-offset))
  lines(c(3.1,3.1),c(y, y-offset))
  text(x[2]+((3.1-x[2])/2),y+offset,"*", cex = 3, font = 2)
  dev.off()
}



#####
########################################
########################################
##### ASE #############################
########################################

library(reshape2)

setwd("D:/Martin Lab/RNA-seq/axm/post_reviews/ASE")
gatk_files <- c("axmJ1", "axmJ4", "axmJ5", "axmJ6", "CPAJ1", "CPAJ2", "CPAJ3", "CPMJ1", "CPMJ2", "CPMJ3", "LLAJ1", "LLAJ2", "LLAJ3", "LLMJ1", "LLMJ2", "LLMJ3", "OAE2", "OAE3", "OAE4", "OME1", "OME2", "OME3", "OUE1", "OUE3", "OUE4", "CAE1", "CAE2", "CAE3", "CME1", "CME2", "CME5")
# make haplotype file

for(file in gatk_files) 
{
  
  infile <- paste(file,"_snp_table.txt",sep="")
  outfile <- paste(file,"_snp_table_haplotypes.file",sep="")
  
  #infile <- "axmJ1_snp_table.txt"
  
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

# make haplotype counts file
for(file in gatk_files) 
{
  
  
  infile1 <- paste(file,"_snp_table_haplotypes.file",sep="")
  infile2 <- paste(file,"_allele_counts.table",sep="")
  outfile <- paste(file,"_haplotype_counts.txt",sep="")
  
  #infile1 <- "axmJ1_snp_table_haplotypes.file"
  #infile2 <- "axmJ1_allele_counts.table"
  
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
  
  write.table(haps_counts,quote=FALSE, row.names = FALSE, col.names = TRUE, outfile)
}

### set minimum coverage
for(file in gatk_files) 
{
  
  setwd("D:/Martin Lab/RNA-seq/axm/post_reviews/ASE")
  infile1 <- paste(file,"_gene_haplotype_counts.txt",sep="")
  outfile <- paste(file,"_min_cov_20.txt",sep="")
  
  haps <- read.table(infile1, header = TRUE, stringsAsFactors = FALSE)
  haps <- haps[which(haps$alleleOneCount >= 10 & haps$alleleTwoCount >= 10),]
  
  setwd("D:/Martin Lab/RNA-seq/axm/post_reviews/ASE/min_cov_20")
  write.table(haps, outfile, quote=FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  
}

# NOTE! use python scripts before moving forward
# merge genes with haplotype counts file (use python snp_genes_overlap_haplotype_counts)

# determine total number of genes that contain a snp with at least
# 10x ref and 10x alt coverage across all samples (total number of genes
# analyzed for ase)

first_samp <- "axmJ1"
first_samp <- "OUE1"
setwd("D:/Martin Lab/RNA-seq/axm/post_reviews/ASE/min_cov_20")
infile1 <- paste(first_samp,"_min_cov_20.txt",sep="")
df1 <- read.table(infile1, header = TRUE, stringsAsFactors = FALSE)
#set_a <- paste(df1$geneIndex, df1$snpIndex, sep = ":")
informative_genes <- df1$geneIndex
informative_genes <- df1$snpIndex

all_samps <- c("axmJ4","axmJ5","axmJ6","CPAJ1","CPAJ2","CPAJ3","CPMJ1","CPMJ2","CPMJ3","LLAJ1","LLAJ2","LLAJ3","LLMJ1","LLMJ2","LLMJ3","OAE2","OAE3","OAE4","OME1","OME2","OME3","OUE1","OUE3","OUE4","CAE1","CAE2","CAE3","CME1","CME2","CME5")
dpf_8 <-  c("OUE3", "OUE4","OAE2", "OAE3", "OAE4", "OME1", "OME2", "OME3", "CAE1", "CAE2", "CAE3", "CME1", "CME2", "CME5")
dpf_17 <- c("axmJ1", "axmJ4", "axmJ5", "axmJ6","CPAJ1", "CPAJ2", "CPAJ3", "CPMJ1", "CPMJ2", "CPMJ3", "LLAJ1", "LLAJ2", "LLAJ3", "LLMJ1", "LLMJ2", "LLMJ3")

for(file in all_samps)
{
setwd("D:/Martin Lab/RNA-seq/axm/post_reviews/ASE/min_cov_20")
infile1 <- paste(file,"_min_cov_20.txt",sep="")
df1 <- read.table(infile1, header = TRUE, stringsAsFactors = FALSE)
set_b <- df1$geneIndex
set_a <- intersect(set_a, set_b)
}
length(set_a)

for(file in all_samps)
{
  setwd("D:/Martin Lab/RNA-seq/axm/post_reviews/ASE/min_cov_20")
  infile1 <- paste(file,"_min_cov_20.txt",sep="")
  df1 <- read.table(infile1, header = TRUE, stringsAsFactors = FALSE)
  #informative_genes <- c(informative_genes, df1$geneIndex)
  informative_genes <- c(informative_genes, df1$snpIndex)
  
}
length(unique(informative_genes))


#ase results w/ binomial comparison

hybrids_8dpf <-  c("OUE1", "OUE3", "OUE4")
hybrids_17dpf <- c("axmJ1", "axmJ4", "axmJ5", "axmJ6")
parents_8dpf <-  c("OAE2", "OAE3", "OAE4", "OME1", "OME2", "OME3", "CAE1", "CAE2", "CAE3", "CME1", "CME2", "CME5")
parents_17dpf <- c("CPAJ1", "CPAJ2", "CPAJ3", "CPMJ1", "CPMJ2", "CPMJ3", "LLAJ1", "LLAJ2", "LLAJ3", "LLMJ1", "LLMJ2", "LLMJ3")

#parents 8-10 dpf e letters indivs (ase called with dedup bams)
#parents_8dpf <-  c("CPAE1","CPAE2","CPAE3","CPME1","CPME2","CPME3","LLAE1","LLAE2","LLAE3","LLME1","LLME2","LLME3")

sum_stats <- c()
parental_ase_genes <- c()
hybrid_ase_genes <- c()
no_hybrid_ase_gene_names <- c()
no_parent_ase_gene_names <- c()
#file <- "axmJ1"

for(file in hybrids_8dpf)
{
  
  setwd("D:/Martin Lab/RNA-seq/axm/post_reviews/ASE/min_cov_20")
  infile1 <- paste(file,"_min_cov_20.txt",sep="")
  outfile <- paste(file,"_binomial.txt",sep="")
  haps <- read.table(infile1, header = TRUE, stringsAsFactors = FALSE)
  haps$totalCount <- haps$alleleOneCount + haps$alleleTwoCount
  haps <- haps[which(haps$totalCount != 0),]
  haps <- haps[which(haps$alleleOneCount >= 100 & haps$ alleleTwoCount >= 100),]
  binom_tests <- c()
  
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
  
  for (j in unique(geneIndex))
  {
    
    geneIndex_table <- haps[which(haps$geneIndex == j),]
    p_vals <- geneIndex_table$binom_test_p_value
    
    # all snps in gene show ase in same haplotype
    #a1s <- c()
    #a2s <- c()
    #for (snp_in_gene in nrow(geneIndex_table))
    #{
    #a1 <- geneIndex_table$alleleOneCount[snp_in_gene]
    #a2 <- geneIndex_table$alleleTwoCount[snp_in_gene]
    #all_snps_a1 <- as.character(a1 > a2)
    #all_snps_a2 <- as.character(a2 < a1)
    #
    #}
    #a1s <- c(a1s, all_snps_a1)
    #a2s <- c(a2s, all_snps_a2)
    #all_ase <- all(a1s == "TRUE") | all(a2s == "TRUE")
    
    if (all(p_vals <= 0.05))# & all_ase == "TRUE")
    {
      counter <- counter + 1
      ase_gene_names <- c(ase_gene_names, j)
    }
    else
    {
      no_hybrid_ase_gene_names <- c(no_hybrid_ase_gene_names, j)
    }
    
  }
  
  ase_genes <- counter
  total_genes <- length(unique(haps$geneIndex))
  sum_stat <- c(file, ase_genes,total_genes,length(unique(haps$snpIndex)))
  sum_stats <- c(sum_stats, sum_stat)
  #parental_ase_genes <- c(parental_ase_genes, ase_gene_names)
  hybrid_ase_genes <- c(hybrid_ase_genes, ase_gene_names)
  print(file)
  print(ase_genes)
  print(total_genes)
  print(ase_genes/total_genes)
  haps <- haps[haps$geneIndex %in% ase_gene_names, ]
  #write.table(haps, outfile, quote=FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  
} 
hybrid_ase_genes_8 <- hybrid_ase_genes
no_hybrid_ase_gene_8 <- no_hybrid_ase_gene_names
for(file in parents_8dpf)
{
  
  setwd("D:/Martin Lab/RNA-seq/axm/post_reviews/ASE/min_cov_20")
  infile1 <- paste(file,"_min_cov_20.txt",sep="")
  outfile <- paste(file,"_binomial.txt",sep="")
  haps <- read.table(infile1, header = TRUE, stringsAsFactors = FALSE)
  haps$totalCount <- haps$alleleOneCount + haps$alleleTwoCount
  haps <- haps[which(haps$totalCount != 0),]
  haps <- haps[which(haps$alleleOneCount >= 100 & haps$ alleleTwoCount >= 100),]
  binom_tests <- c()
  
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
  
  for (j in unique(geneIndex))
  {
    
    geneIndex_table <- haps[which(haps$geneIndex == j),]
    p_vals <- geneIndex_table$binom_test_p_value
    
    # all snps in gene show ase in same haplotype
    #a1s <- c()
    #a2s <- c()
    #for (snp_in_gene in nrow(geneIndex_table))
    #{
    #  a1 <- geneIndex_table$alleleOneCount[snp_in_gene]
    #  a2 <- geneIndex_table$alleleTwoCount[snp_in_gene]
    #  all_snps_a1 <- as.character(a1 > a2)
    #  all_snps_a2 <- as.character(a2 < a1)
    #  
    #}
    #a1s <- c(a1s, all_snps_a1)
    #a2s <- c(a2s, all_snps_a2)
    #all_ase <- all(a1s == "TRUE") | all(a2s == "TRUE")
    
    if (all(p_vals <= 0.05))# & all_ase == "TRUE")
    {
      counter <- counter + 1
      ase_gene_names <- c(ase_gene_names, j)
    }
    else
    {
      no_parent_ase_gene_names <- c(no_parent_ase_gene_names, j)
    }
    
  }
  
  ase_genes <- counter
  total_genes <- length(unique(haps$geneIndex))
  sum_stat <- c(file, ase_genes,total_genes,length(unique(haps$snpIndex)))
  sum_stats <- c(sum_stats, sum_stat)
  parental_ase_genes <- c(parental_ase_genes, ase_gene_names)
  #hybrid_ase_genes <- c(hybrid_ase_genes, ase_gene_names)
  print(file)
  print(ase_genes)
  print(total_genes)
  print(ase_genes/total_genes)
  haps <- haps[haps$geneIndex %in% ase_gene_names, ]
  #write.table(haps, outfile, quote=FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  
} 
parental_ase_genes_8 <- parental_ase_genes
no_parent_ase_gene_8 <- no_parent_ase_gene_names

parental_ase_genes <- c()
hybrid_ase_genes <- c()
no_hybrid_ase_gene_names <- c()
no_parent_ase_gene_names <- c()

for(file in hybrids_17dpf)
{
  
  setwd("D:/Martin Lab/RNA-seq/axm/post_reviews/ASE/min_cov_20")
  infile1 <- paste(file,"_min_cov_20.txt",sep="")
  outfile <- paste(file,"_binomial.txt",sep="")
  haps <- read.table(infile1, header = TRUE, stringsAsFactors = FALSE)
  haps$totalCount <- haps$alleleOneCount + haps$alleleTwoCount
  haps <- haps[which(haps$totalCount != 0),]
  haps <- haps[which(haps$alleleOneCount >= 100 & haps$ alleleTwoCount >= 100),]
  binom_tests <- c()
  
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
  
  for (j in unique(geneIndex))
  {
    
    geneIndex_table <- haps[which(haps$geneIndex == j),]
    p_vals <- geneIndex_table$binom_test_p_value
    
      # all snps in gene show ase in same haplotype
      #a1s <- c()
      #a2s <- c()
      #for (snp_in_gene in nrow(geneIndex_table))
      #{
      #a1 <- geneIndex_table$alleleOneCount[snp_in_gene]
      #a2 <- geneIndex_table$alleleTwoCount[snp_in_gene]
      #all_snps_a1 <- as.character(a1 > a2)
      #all_snps_a2 <- as.character(a2 < a1)
      #
      #}
      #a1s <- c(a1s, all_snps_a1)
      #a2s <- c(a2s, all_snps_a2)
      #all_ase <- all(a1s == "TRUE") | all(a2s == "TRUE")
    
    if (all(p_vals <= 0.05))# & all_ase == "TRUE")
    {
      counter <- counter + 1
      ase_gene_names <- c(ase_gene_names, j)
    }
    else
    {
      no_hybrid_ase_gene_names <- c(no_hybrid_ase_gene_names, j)
    }
    
  }
  
  ase_genes <- counter
  total_genes <- length(unique(haps$geneIndex))
  sum_stat <- c(file, ase_genes,total_genes,length(unique(haps$snpIndex)))
  sum_stats <- c(sum_stats, sum_stat)
  #parental_ase_genes <- c(parental_ase_genes, ase_gene_names)
  hybrid_ase_genes <- c(hybrid_ase_genes, ase_gene_names)
  print(file)
  print(ase_genes)
  print(total_genes)
  print(ase_genes/total_genes)
  haps <- haps[haps$geneIndex %in% ase_gene_names, ]
  #write.table(haps, outfile, quote=FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  
} 
hybrid_ase_genes_17 <- hybrid_ase_genes
no_hybrid_ase_gene_17 <- no_hybrid_ase_gene_names
for(file in parents_17dpf)
{
  
  setwd("D:/Martin Lab/RNA-seq/axm/post_reviews/ASE/min_cov_20")
  infile1 <- paste(file,"_min_cov_20.txt",sep="")
  outfile <- paste(file,"_binomial.txt",sep="")
  haps <- read.table(infile1, header = TRUE, stringsAsFactors = FALSE)
  haps$totalCount <- haps$alleleOneCount + haps$alleleTwoCount
  haps <- haps[which(haps$totalCount != 0),]
  haps <- haps[which(haps$alleleOneCount >= 100 & haps$ alleleTwoCount >= 100),]
  binom_tests <- c()
  
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
  
  for (j in unique(geneIndex))
  {
    
    geneIndex_table <- haps[which(haps$geneIndex == j),]
    p_vals <- geneIndex_table$binom_test_p_value
    
    # all snps in gene show ase in same haplotype
    #a1s <- c()
    #a2s <- c()
    #for (snp_in_gene in nrow(geneIndex_table))
    #{
    #  a1 <- geneIndex_table$alleleOneCount[snp_in_gene]
    #  a2 <- geneIndex_table$alleleTwoCount[snp_in_gene]
    #  all_snps_a1 <- as.character(a1 > a2)
    #  all_snps_a2 <- as.character(a2 < a1)
    #  
    #}
    #a1s <- c(a1s, all_snps_a1)
    #a2s <- c(a2s, all_snps_a2)
    #all_ase <- all(a1s == "TRUE") | all(a2s == "TRUE")
    
    if (all(p_vals <= 0.05))# & all_ase == "TRUE")
    {
      counter <- counter + 1
      ase_gene_names <- c(ase_gene_names, j)
    }
    else
    {
      no_parent_ase_gene_names <- c(no_parent_ase_gene_names, j)
    }
    
  }
  
  ase_genes <- counter
  total_genes <- length(unique(haps$geneIndex))
  sum_stat <- c(file, ase_genes,total_genes,length(unique(haps$snpIndex)))
  sum_stats <- c(sum_stats, sum_stat)
  parental_ase_genes <- c(parental_ase_genes, ase_gene_names)
  #hybrid_ase_genes <- c(hybrid_ase_genes, ase_gene_names)
  print(file)
  print(ase_genes)
  print(total_genes)
  print(ase_genes/total_genes)
  haps <- haps[haps$geneIndex %in% ase_gene_names, ]
  #write.table(haps, outfile, quote=FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  
} 
parental_ase_genes_17 <- parental_ase_genes
no_parent_ase_gene_17 <- no_parent_ase_gene_names

length(unique(parental_ase_genes_8))
length(unique(no_parent_ase_gene_8))
length(unique(hybrid_ase_genes_8))
length(unique(no_hybrid_ase_gene_8))
length(unique(parental_ase_genes_17))
length(unique(no_parent_ase_gene_17))
length(unique(hybrid_ase_genes_17))
length(unique(no_hybrid_ase_gene_17))

length(setdiff(hybrid_ase_genes_8, parental_ase_genes_8))
length(setdiff(hybrid_ase_genes_17, parental_ase_genes_17))

ase_in_hybrids_only_8   <- as.data.frame(setdiff(hybrid_ase_genes_8, parental_ase_genes_8))
ase_in_hybrids_only_17  <- as.data.frame(setdiff(hybrid_ase_genes_17, parental_ase_genes_17))
no_ase_or_parent_ase_8  <- as.data.frame(unique(setdiff(c(no_hybrid_ase_gene_8, parental_ase_genes_8),hybrid_ase_genes_8)))
no_ase_or_parent_ase_17 <- as.data.frame(unique(setdiff(c(no_hybrid_ase_gene_17, parental_ase_genes_17),hybrid_ase_genes_17)))
ase_in_hybrids_only_8$tag   <- ase_in_hybrids_only_8$`setdiff(hybrid_ase_genes_8, parental_ase_genes_8)`
ase_in_hybrids_only_17$tag  <- ase_in_hybrids_only_17$`setdiff(hybrid_ase_genes_17, parental_ase_genes_17)`
no_ase_or_parent_ase_8$tag  <- no_ase_or_parent_ase_8$`unique(setdiff(c(no_hybrid_ase_gene_8, parental_ase_genes_8), hybrid_ase_genes_8))`
no_ase_or_parent_ase_17$tag <- no_ase_or_parent_ase_17$`unique(setdiff(c(no_hybrid_ase_gene_17, parental_ase_genes_17), hybrid_ase_genes_17))`
ase_in_hybrids_only_8$`setdiff(hybrid_ase_genes_8, parental_ase_genes_8)` <- NULL
ase_in_hybrids_only_17$`setdiff(hybrid_ase_genes_17, parental_ase_genes_17)`<- NULL
no_ase_or_parent_ase_8$`unique(setdiff(c(no_hybrid_ase_gene_8, parental_ase_genes_8), hybrid_ase_genes_8))`<- NULL
no_ase_or_parent_ase_17$`unique(setdiff(c(no_hybrid_ase_gene_17, parental_ase_genes_17), hybrid_ase_genes_17))`<- NULL
ase_in_hybrids_only_8$ase_in_hybrids_not_parents   <- "yes"
ase_in_hybrids_only_17$ase_in_hybrids_not_parents  <- "yes"
no_ase_or_parent_ase_8$ase_in_hybrids_not_parents  <- "no"
no_ase_or_parent_ase_17$ase_in_hybrids_not_parents <-"no"
ase_8 <- rbind(ase_in_hybrids_only_8,no_ase_or_parent_ase_8)
ase_17 <- rbind(ase_in_hybrids_only_17,no_ase_or_parent_ase_17)

inheritance_8dpf  <- read.csv("D:/Martin Lab/RNA-seq/axm/post_reviews/inheritance_patterns/OAxOM_8dpf_inheritance.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
inheritance_17dpf <- read.csv("D:/Martin Lab/RNA-seq/axm/post_reviews/inheritance_patterns/LAxLM_17dpf_inheritance.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
head(inheritance_17dpf)
nrow(inheritance_17dpf)
nrow(inheritance_8dpf)

inheritance_and_ase_8 <- merge(inheritance_8dpf, ase_8, all.x = TRUE,by = c("tag"))
inheritance_and_ase_8[is.na(inheritance_and_ase_8)] <- "not_measured"
head(inheritance_and_ase_8)
nrow(inheritance_and_ase_8[which(inheritance_and_ase_8$ase_in_hybrids_not_parents =="yes"),])

inheritance_and_ase_8$tag[duplicated(inheritance_and_ase_8$tag)]

#write.table(inheritance_and_ase_8, "D:/Martin Lab/RNA-seq/axm/post_reviews/inheritance_patterns/inheritance_and_ase_8.txt", quote=FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


inheritance_and_ase_17 <- merge(inheritance_17dpf, ase_17, all.x = TRUE,by = c("tag"))
inheritance_and_ase_17[is.na(inheritance_and_ase_17)] <- "not_measured"
head(inheritance_and_ase_17)
#write.table(inheritance_and_ase_17, "D:/Martin Lab/RNA-seq/axm/post_reviews/inheritance_patterns/inheritance_and_ase_17.txt", quote=FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

inheritance_and_ase_8dpf  <- read.csv("D:/Martin Lab/RNA-seq/axm/post_reviews/inheritance_patterns/inheritance_and_ase_8.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
inheritance_and_ase_17dpf  <- read.csv("D:/Martin Lab/RNA-seq/axm/post_reviews/inheritance_patterns/inheritance_and_ase_17.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")

inheritance_and_ase <- inheritance_and_ase_17dpf
#de_am_8dpf <- read.csv("D:/Martin Lab/RNA-seq/axm/post_reviews/conditions/DE_(CRPA2_&_OSPA)_vs_(CRPM2_&_OSPM)_8dpf_genes.csv", header = TRUE, stringsAsFactors = FALSE)
#de_am_17dpf <- read.csv("D:/Martin Lab/RNA-seq/axm/post_reviews/conditions/DE_(CRPA1_&_LILA)_vs_(CRPM1_&_LILM)_17dpf_genes.csv", header = TRUE, stringsAsFactors = FALSE)
#inheritance_and_ase <- subset(inheritance_and_ase, !duplicated(inheritance_and_ase[,1]))
nrow(inheritance_and_ase[which(inheritance_and_ase$padj_ph < 0.05 & inheritance_and_ase$padj_p1p2 < 0.05),])
#inheritance_and_ase <- inheritance_and_ase[which(inheritance_and_ase$padj_p1p2 < 0.05),]
nrow(inheritance_and_ase[which(inheritance_and_ase$ase_in_hybrids_not_parents == "yes" | inheritance_and_ase$ase_in_hybrids_not_parents == "no"),])

trans <- inheritance_and_ase[which(inheritance_and_ase$padj_p1p2 <= 0.05 & inheritance_and_ase$ase_in_hybrids_not_parents == "no"),]
nrow(trans)
length(unique(trans$tag))
trans$ase_type <- 1
trans$ase_shape <- 1
mis_over <- inheritance_and_ase[which(inheritance_and_ase$padj_ph <= 0.05 & inheritance_and_ase$lfc_ph > 0 & inheritance_and_ase$ase_in_hybrids_not_parents == "no"),]
mis_over$ase_type <- 2
mis_over$ase_shape <- 1
nrow(mis_over)
mis_under <- inheritance_and_ase[which(inheritance_and_ase$padj_ph <= 0.05 & inheritance_and_ase$lfc_ph < 0 & inheritance_and_ase$ase_in_hybrids_not_parents == "no"),]
mis_under$ase_type <- 3
mis_under$ase_shape <- 1
nrow(mis_under)
mis_ase_over <- inheritance_and_ase[which(inheritance_and_ase$padj_p1p2 > 0.05 & inheritance_and_ase$padj_ph <= 0.05 & inheritance_and_ase$lfc_ph > 0 & inheritance_and_ase$ase_in_hybrids_not_parents == "yes"),]
nrow(mis_ase_over)
mis_ase_over$ase_type <- 5
mis_ase_over$ase_shape <- 3
mis_ase_under <- inheritance_and_ase[which(inheritance_and_ase$padj_p1p2 > 0.05 & inheritance_and_ase$padj_ph <= 0.05 & inheritance_and_ase$lfc_ph < 0 & inheritance_and_ase$ase_in_hybrids_not_parents == "yes"),]
nrow(mis_ase_under)
mis_ase_under$ase_type <- 6
mis_ase_under$ase_shape <- 3
#mis_ase_over <- inheritance_and_ase[which(inheritance_and_ase$padj_p1p2 > 0.05  & inheritance_and_ase$lfc_ph > 0 & inheritance_and_ase$ase_in_hybrids_not_parents == "yes"),]
#nrow(mis_ase_over)
#mis_ase_over$ase_type <- 6
#mis_ase_under <- inheritance_and_ase[which(inheritance_and_ase$padj_p1p2 > 0.05 & inheritance_and_ase$lfc_ph < 0 & inheritance_and_ase$ase_in_hybrids_not_parents == "yes"),]
#nrow(mis_ase_under)
#mis_ase_under$ase_type <- 6

comp <- inheritance_and_ase[which(inheritance_and_ase$padj_p1p2 > 0.05 & inheritance_and_ase$ase_in_hybrids_not_parents == "yes"),]
comp$ase_type <- 4
comp$ase_shape <- 2
nrow(comp)

con <- inheritance_and_ase[which(inheritance_and_ase$padj_ph > 0.05 & inheritance_and_ase$padj_p1p2 > 0.05 & inheritance_and_ase$ase_in_hybrids_not_parents == "no"),]
nrow(con)
con$ase_type <- 7
con$ase_shape <- 1
cis <- inheritance_and_ase[which(inheritance_and_ase$padj_p1p2 <= 0.05 & inheritance_and_ase$ase_in_hybrids_not_parents == "yes"),]
nrow(cis)
cis$ase_type <- 8
cis$ase_shape <- 2

ase_type <- rbind(con, mis_over, mis_under,cis, trans,comp, mis_ase_over, mis_ase_under)
head(ase_type)
nrow(ase_type)
unique(ase_type$ase_type)
length(unique(ase_type$tag))

#yellow ="#f6c700"
#orange = "#E77200"
#pink = "#FF00CC"
#black = "#000000"
#red = #bd1e24
#blue = #0067a7
#

cols <- 	c("#f6c700", "#bd1e24", "#0067a7", "#FF00CC","#E77200", "#E77200","#f6c700","#f6c700")

cols_ase <- cols[ase_type$ase_type]

total_genes <- nrow(ase_type)
nrow(ase_type[which(ase_type$padj_p1p2 > 0.05),])
nrow(ase_type[which(ase_type$padj_p1p2 <= 0.05),])
con <- paste(round(((nrow(ase_type[which(ase_type$ase_type == 1),]) / total_genes) *100), digits = 2), "%", sep = "")
mis_over <-  paste(round(((nrow(ase_type[which(ase_type$ase_type == 2),]) / total_genes) *100), digits = 2), "%", sep = "")
mis_under <- paste(round(((nrow(ase_type[which(ase_type$ase_type == 3),]) / total_genes) *100), digits = 2), "%", sep = "")
cis <-   paste(round(((nrow(ase_type[which(ase_type$ase_type == 4),]) / total_genes) *100), digits = 2), "%", sep = "")
trans <-     paste(round(((nrow(ase_type[which(ase_type$ase_type == 5),]) / total_genes) *100), digits = 2), "%", sep = "")
comp <-   paste(round(((nrow(ase_type[which(ase_type$ase_type == 6),]) / total_genes) *100), digits = 2), "%", sep = "")
cols_l <-	c("#000000","#E77200","#f6c700","#FF00CC","#bd1e24","#0067a7")


tiff("D:/Martin Lab/RNA-seq/axm/post_reviews/ASE/plots/compensatory_17dpf2.tiff", width = 6, height = 6, units = 'in', res = 1000)
plot(ase_type$lfc_p1p2, ase_type$lfc_ph, col = cols_ase, 
     pch = c(16,17,17)[as.factor(ase_type$ase_shape)],
     cex = c(1, 1, 1, 1,1.1,1.1,1,1)[ase_type$ase_type],cex.axis=1.5, ylab = "", xlab = "")#, ylim = c(-2,1.6),xlim = c(-2,1.6))
     #ylab = "log2 fold change parental species vs. hybrids",
     #xlab = "log2 fold change generalists vs. molluscivores")
#legend("bottomleft", inset=0,
#       c(paste("cis", cis),
#         paste("trans", trans),
#         paste("conserved", con),
#         paste("compensatory", comp),
#         paste("overdominant", mis_over),
#         paste("underdominant", mis_under)),
#         cex = 1.4, fill=cols_l, horiz=FALSE, bty="n")
abline(v=0, col = "black", lty = 3, lwd = 1.8)
abline(h=0, col = "black", lty = 3, lwd = 1.8)
dev.off()

cols_b <- c("#f6c700", "#FF00CC","#0067a7", "#bd1e24","#E77200","#000000")
labs <- c("additive","generalist dominant","molluscivore dominant","overdominant","underdominant", "conserved")
bars <- c(((nrow(ase_type[which(ase_type$ase_type == 1),]) / total_genes) *100), ((nrow(ase_type[which(ase_type$ase_type == 6),]) / total_genes) *100), ((nrow(ase_type[which(ase_type$ase_type == 3),]) / total_genes) *100),((nrow(ase_type[which(ase_type$ase_type == 2),]) / total_genes) *100),((nrow(ase_type[which(ase_type$ase_type == 5),]) / total_genes) *100),((nrow(ase_type[which(ase_type$ase_type == 4),]) / total_genes) *100))

#tiff("D:/Martin Lab/RNA-seq/axm/post_reviews/ASE/plots/17dpf_bar.tiff", width = 6, height = 4, units = 'in', res = 1000)
b1 <- barplot(bars, las =1, horiz = TRUE, font.axis = 2, cex.axis=1.5, xlim = c(0,60), col = cols_b)
barplot(bars, las =1, horiz = TRUE, font.axis = 2, cex.axis=1.5, xlim = c(0,60), col = cols_b)

offset <- 7
text(y=b1[1,], x=((nrow(ase_type[which(ase_type$ase_type == 1),]) / total_genes) *100)+offset,con, font = 2, cex = 1)
text(y=b1[2,], x=((nrow(ase_type[which(ase_type$ase_type == 6),]) / total_genes) *100)+offset,comp, font = 2, cex = 1)
text(y=b1[3,], x=((nrow(ase_type[which(ase_type$ase_type == 3),]) / total_genes) *100)+offset,mis_under, font = 2, cex = 1)
text(y=b1[4,], x=((nrow(ase_type[which(ase_type$ase_type == 2),]) / total_genes) *100)+offset,mis_over, font = 2, cex = 1)
text(y=b1[5,], x=((nrow(ase_type[which(ase_type$ase_type == 5),]) / total_genes) *100)+offset,trans, font = 2, cex = 1)
text(y=b1[6,], x=((nrow(ase_type[which(ase_type$ase_type == 4),]) / total_genes) *100)+offset,cis, font = 2, cex = 1)
#dev.off()

head(ase_type)
cis <- ase_type[which(ase_type$inheritance == 4),]
trans <- ase_type[which(ase_type$inheritance == 5),]
both <- rbind(cis, trans)
both$lfc_p1p2 <- abs(both$lfc_p1p2)
head(both)
boxplot(both$lfc_p1p2~as.character(both$inheritance), main = "")
wilcox.test(both$lfc_p1p2~as.character(both$inheritance))



prop_ase  <- read.table("D:/Martin Lab/RNA-seq/axm/post_reviews/ASE/min_cov_20/proportion_of_ase_genes.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
#prop_ase <- prop_ase[which(prop_ase$sample != "axmJ1" &prop_ase$sample != "axmJ4" & prop_ase$sample != "axmJ5" & prop_ase$sample != "axmJ6"),]
prop_ase_17 <- prop_ase[which(prop_ase$stage == "17dpf"),]
prop_ase_8 <- prop_ase[which(prop_ase$stage == "8dpf"),]
head(prop_ase)
par(mai=c(3,1,1,1))
cols = c("#bd1e24","#800080","#0067a7")[as.factor(prop_ase$species)]
dens <- c(30,1000)[as.factor(prop_ase$stage)]

tiff("D:/Martin Lab/RNA-seq/axm/post_reviews/ASE/plots/prop_ase_cov_200.tiff", width = 7.5, height = 3.5, units = 'in', res = 1000)
barplot(prop_ase$prop_ase_genes_100,density = dens,xlab = "",col = cols, cex.names=0.7, main="", ylim = c(0,0.6),
        names.arg=prop_ase$lake, ylab = "proportion of genes showing ASE")
title(xlab="sample", line=2)
l_cols <- c("black", "black", "#bd1e24","#0067a7","#800080")
#legend( "topright",c("8 dpf","17-20 dpf", "generalists", "molluscivores", "hybrids"),fill=l_cols, bty = "n", density = c(30,30, 30, 30, 30))
dev.off()
tiff("D:/Martin Lab/RNA-seq/axm/post_reviews/ASE/plots/prop_ase_legend_2.tiff", width = 7.5, height = 3.5, units = 'in', res = 1000)
barplot(prop_ase$prop_ase_genes,density = dens,xlab = "",col = cols, cex.names=0.7, main="",ylim = c(0,0.6), 
        names.arg=prop_ase$lake, ylab = "proportion of genes showing ASE")
title(xlab="sample", line=2)
l_cols <- c("black", "black", "#bd1e24","#0067a7","#800080")
#legend( "topright",c("8 dpf","17-20 dpf", "generalists", "molluscivores", "hybrids"),fill=l_cols, bty = "n", density = c(1000,1000, 1000, 1000, 1000))
dev.off()


hy <- c("axmJ1","axmJ4","axmJ5","axmJ6")
prop_ase <- prop_ase[!(prop_ase$sample %in% hy),]
fit <- lm(prop_ase$prop_ase_genes~prop_ase$tin_median)
plot(prop_ase$tin_median, prop_ase$prop_ase_genes, xlim = c(30,90), ylim = c(0,.75))
abline(fit)
summary.lm(fit)
predict(fit, data.frame(tin_median= c(30,40,50)))
range(predict(fit, data.frame(tin_median= 32.68)))

cols <- 	c("#000000","#bd1e24")
options(scipen=999)
options(scipen=10)
cols_in <- cols[as.factor(prop_ase$kit)]

resp <- prop_ase$prop_ase_genes
trms <- prop_ase$tin_median
x_lab <- "TIN"
y_lab <- "normalized counts"
out <-   "tin_prop_ase_all"
qcplot <- paste("D:/Martin Lab/RNA-seq/axm/post_reviews/qc_stats/plots/", out, ".tiff", sep = "")

fit <- lm(prop_ase$prop_ase_genes~prop_ase$tin_median)
summary(fit)
modsum <- summary(fit)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]

#tiff(qcplot, width = 5.5, height = 4.5, units = 'in', res = 1000)
par(mai=c(1,1,1,1))
plot    (trms, resp,pch = c(19, 17, 15)[as.factor(prop_ase$species)], col = cols_in,
         xlab = x_lab, ylab = y_lab, las = 1, ylim = c(0,0.6), xlim = c(30,90))
abline(fit)
mylabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))
text(x = 19, y = 2.5, labels = mylabel)
legend('topright', legend = rp, bty = 'n')#,inset=c(-0.03,0))
#dev.off()


ggplot(data = prop_ase, aes(x = tin_median, y = prop_ase_genes)) + geom_point() + 
  xlim(30,90) + 
  ylim(0,.75) +
  stat_smooth(method = "lm", col = "dodgerblue3") +
  theme(panel.background = element_rect(fill = "white"),
        axis.line.x=element_line(),
        axis.line.y=element_line()) 

fit <- aov(prop_ase_genes_100 ~ species + stage, data=prop_ase)
summary.lm(fit)
fit <- aov(prop_ase_genes ~ species, data=prop_ase_17)
summary.lm(fit)
fit <- aov(prop_ase_genes ~ species, data=prop_ase_8)
summary.lm(fit)



inheritance_8dpf  <- read.csv("D:/Martin Lab/RNA-seq/axm/post_reviews/inheritance_patterns/OAxOM_8dpf_inheritance.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
inheritance_17dpf <- read.csv("D:/Martin Lab/RNA-seq/axm/post_reviews/inheritance_patterns/LAxLM_17dpf_inheritance.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
inheritance <- inheritance_17dpf


down <- inheritance[which(inheritance$inheritance == 6),]
down$type <- "underdominant"
#down$type <- 2
down$lfc_ph <- abs(down$lfc_ph)
up <- inheritance[which(inheritance$inheritance == 5),]
up$type <- "overdominant"
#up$type <- 1
mags <- rbind(up, down)
head(mags)
nrow(mags)
boxplot(mags$lfc_ph~mags$inheritance, main = "")
wilcox.test(mags$lfc_ph~mags$inheritance)


#tiff("D:/Martin Lab/RNA-seq/axm/post_reviews/over_under_magnitude_17dpf.tif", width = 5.5, height = 6, units = 'in', res = 600)
p <- ggplot(mags, aes(x=type, y=lfc_ph, fill = type)) + 
  geom_violin(trim = TRUE)
p+ scale_fill_manual(values=c("#bd1e24", "#0067a7")) + theme(axis.text.y = element_text(size=18, color = '#000000'),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                             panel.background = element_blank(),axis.line = element_line(colour = "black"),
                                                             axis.text.x = element_text(size=18, color = '#000000'), legend.position="none",
                                                             axis.title.y=element_blank(),axis.title.x=element_blank()) + stat_summary(fun.data=data_summary, col = "black")
#dev.off()



#####
########################################
########################################
##### PCA and regressions ##############
########################################

#### all samples dds ###
library("DESeq2")
library(pheatmap)
library(vegan)
library(rgl)
library(ape)
library(ggplot2)
library(reshape2)
library(gplots)
library(RColorBrewer)

#mirror_design <- c("Geneid", "axmJ1", "axmJ4", "axmJ5", "axmJ6", "CPAJ1", "CPAJ2", "CPAJ3", "CPMJ1", "CPMJ2", "CPMJ3", "LLAJ1", "LLAJ2", "LLAJ3", "LLMJ1", "LLMJ2", "LLMJ3", "CAE1", "CAE2", "CAE3", "CME1", "CME2", "CME5", "OAE2", "OAE3", "OAE4", "OME1", "OME2", "OME3", "OUE1", "OUE3", "OUE4", "CPAE1", "CPAE2", "CPAE3", "CPME1", "CPME2", "CPME3", "LLAE1", "LLAE2", "LLAE3", "LLME1", "LLME2", "LLME3")
#cts <- read.table("D:/Martin Lab/RNA-seq/axm/post_reviews/all_rna_4_rounds_counts.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
#cts <- cts[mirror_design]
#write.table(cts, "D:/Martin Lab/RNA-seq/axm/post_reviews/mirror_design_counts.txt", row.names = FALSE, quote = FALSE, sep = "\t")
cts <- as.matrix(read.table("D:/Martin Lab/RNA-seq/axm/post_reviews/mirror_design_counts.txt",header = TRUE,row.names=1))
head(cts)
nrow(cts)
colData <- as.matrix(read.table("D:/Martin Lab/RNA-seq/axm/post_reviews/table_maker_mirror_design.txt",header = TRUE,row.names=1))
head(colData)
ncol(cts)
nrow(colData)
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = colData,
                              design= ~f1)
dds <- estimateSizeFactors(dds)
idx <- rowSums(counts(dds, normalized=TRUE) >= 10 ) >= (nrow(colData))
length(rowSums(counts(dds, normalized=TRUE) >= 10 ) >= (nrow(colData)))
length(idx) 
dds <- dds[idx,]
dds <- DESeq(dds)
norm_cts <- counts(dds, normalized=TRUE)
head(norm_cts)
norm_cts <- as.data.frame(norm_cts)
#write.table(norm_cts, "D:/Martin Lab/RNA-seq/axm/post_reviews/normalized_counts.txt", , sep = "\t", quote = FALSE, row.names = FALSE)

# PCA
rld <- rlog(dds)
dists <- dist(t(assay(rld)))
plot(hclust(dists))

#dds.sub <- dds[ , dds$stage %in% c("17dpf") ]
dds.sub <- dds[ , dds$stage %in% c("8dpf") ]
rld.sub <- rld[ , rld$stage %in% c("8dpf") ]
plotPCA(rld, intgroup=c("sequencing_date"))+scale_shape_manual(values=c(25))
geom_point(aes(shape='species'))


cols_in <- brewer.pal(10, "Paired")
cols_in <- c("#bd1e24","#800080","#0067a7")
#tiff("D:/Martin Lab/RNA-seq/axm/post_reviews/supp_plots/pca_8dpf.tiff", width = 7, height = 6, units = 'in', res = 1000)
pcaData <- plotPCA(rld, intgroup=c("species",'sequencing_round_num'), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=species, shape=sequencing_round_num)) +
  scale_color_manual(values = cols_in)+
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.key=element_blank()) 
#dev.off()

# Plot counts for a gene
#17
cols_plot_counts <- c("blue4","blue4","blue4","blue4", "red","red","red","green4","green4","green4", "red","red","red","green4","green4","green4")
#8
cols_plot_counts <- c("red","red","red","green4","green4","green4", "red","red","red","green4","green4","green4","blue4","blue4","blue4","red","red","red","green4","green4","green4")

plotCounts(dds.sub, gene="gene7154", intgroup= "species", pch = 16, main = "",col = cols_plot_counts)

# regressions
cts_table <- read.table("D:/Martin Lab/RNA-seq/axm/post_reviews/mirror_design_counts.txt",header = TRUE, stringsAsFactors = FALSE)
info <- read.table("D:/Martin Lab/RNA-seq/axm/post_reviews/table_maker_mirror_design.txt", header = TRUE, stringsAsFactors = FALSE)
norm_cts <- read.table("D:/Martin Lab/RNA-seq/axm/post_reviews/normalized_counts.txt", header = TRUE, stringsAsFactors = FALSE)
head(cts_table)
head(info)
cts_table$Geneid <- NULL
totals <- data.frame(sample = info$sample,total_raw_cts=colSums(cts_table))
info <- merge(info, totals, by = c("sample"))
totals <- data.frame(sample = info$sample,total_norm_cts=colSums(norm_cts))
info <- merge(info, totals, by = c("sample"))
head(info)
gc_meds <- c()
dups <- c()
reads <- c()
depths <- c()
tins <- c()

#i <- "axmJ1"
for (i in info$sample)
{
  
  gc_file <- paste("D:/Martin Lab/RNA-seq/axm/post_reviews/qc_stats/gc_content/",i, "_gc.GC.xls", sep = "")
  gc <- read.table(gc_file, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
  gc_meds <- c(gc_meds, median(gc$GC.))
  
  tin_file <- paste("D:/Martin Lab/RNA-seq/axm/post_reviews/qc_stats/tins/",i, ".rna.sort.summary.txt", sep = "")
  tin <- read.table(tin_file, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
  tins <- c(tins, tin$TIN.median.)
  
  dup_file <- paste("D:/Martin Lab/RNA-seq/axm/post_reviews/qc_stats/dups/",i, ".rna.sort.bam.seq.DupRate.xls", sep = "")
  dup <- read.table(dup_file, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
  head(dup)
  non_dup <- dup[1,2]
  dup <- dup[-1,]
  non_dup - sum(dup$UniqReadNumber)
  sum(dup$UniqReadNumber) / non_dup
  sum(dup$UniqReadNumber)
  non_dup + sum(dup$UniqReadNumber)
  dups <- c(dups, (sum(dup$UniqReadNumber)))
  
  read_file <- paste("D:/Martin Lab/RNA-seq/axm/post_reviews/qc_stats/stats/",i, "stats.txt", sep = "")
  read <- read.table(read_file, header = FALSE, stringsAsFactors = FALSE, sep = "\t", nrows = 1)
  head(read)
  read <- cbind(read, colsplit(read$V1, "\\+", c("reads", "junk")))
  reads <- c(reads, read$reads)
  
  depth_file <- paste("D:/Martin Lab/RNA-seq/axm/post_reviews/qc_stats/depth/",i, "_avg_depth_across_features.txt", sep = "")
  depth <- read.table(depth_file, header = FALSE, stringsAsFactors = FALSE, sep = "\t", nrows = 1)
  head(depth)
  depths <- c(depths, depth$V1)
  
}
info$median_gc_content <- gc_meds
info$duplicate_reads <- dups
info$reads_mapped <- reads
info$prop_duplicate_reads <- info$duplicate_reads/info$reads_mapped
info$avg_depth_across_features <- depths
info$raw_fastq_reads <- (info$raw_fastq_lines / 4) * 2
#info$prop_reads_mapped <- info$reads_mapped / info$raw_fastq_reads
info$tin_median <- tins
head(info)
info$type_stage <- paste(info$cross_type, info$stage)
#write.table(info, "D:/Martin Lab/RNA-seq/axm/post_reviews/qc_stats/all_qc_stats.txt", quote = FALSE, row.names = FALSE, sep = "\t")


cols <- 	c("#f6c700", "#E77200", "#FF00CC","#000000", "#bd1e24", "#0067a7")
cols <- 	c("#000000","#bd1e24")
#b,c,d,a
#con = yellow ="#f6c700"
#add = orange = "#E77200"
#p1_dom = pink = "#FF00CC"
#p2_dom = black = "#000000"
#over_dom = red = #bd1e24
#under_dom = blue = #0067a7
options(scipen=999)
options(scipen=0)
cols_in <- cols[as.factor(info$kit)]


resp <- info$total_norm_cts
trms <- info$tin_median
x_lab <- "TIN"
y_lab <- "normalized counts"
out <-   "tin"
qcplot <- paste("D:/Martin Lab/RNA-seq/axm/post_reviews/qc_stats/plots/", out, ".tiff", sep = "")

lm <- lm(resp~trms)
summary(lm)
modsum <- summary(lm)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]

#tiff(qcplot, width = 4.5, height = 4.5, units = 'in', res = 600)
par(mai=c(1,1,1,1))
plot    (trms, resp,pch = c(19, 17, 15)[as.factor(info$species)], col = cols_in,
         xlab = x_lab, ylab = "", las = 1)
abline(lm)
mylabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))
text(x = 19, y = 2.5, labels = mylabel)
legend('topleft', legend = rp, bty = 'n',inset=c(-0.03,0))
#legend("topright", title = c("Library Prep Kit"),
#       legend = c("TruSeq", "KAPA"), 
#       inset=c(-0.25,0),bty = "n", xpd=TRUE, mar(c(7,7,7,10)), cex = 0.85, 
#       fill = c("#bd1e24","#000000"))
#legend("bottomright", title = c("Species"),
#       legend = c("hybrid", "generalist", "molluscivore"), 
#       inset=c(-0.28,0),bty = "n", xpd=TRUE, mar(c(7,7,7,10)), cex = 0.9, 
#       pch = c(17, 19, 15), col = "black")
#dev.off()


library(ggplot2)
library(ggpubr)
require(gridExtra)
#c('****' = 1e-04, '***' = 0.001, '**' = 0.01, '*' = 0.05,
#install.packages("ggpubr")

cols <- c("black","#bd1e24")

box_plts <- info[order(info$sequencing_round),] 
p1 <- ggboxplot(box_plts, "kit", "total_norm_cts",color = cols, ylab = "normalized counts\n", xlab="library kit")+
  scale_y_continuous(labels = scales::scientific)+
  stat_compare_means(method = "t.test", label.y = max(info$total_norm_cts)+(max(info$total_norm_cts)*.1),label.x = 2.25)+      # Add global p-value
  stat_compare_means(label = "p.signif", method = "t.test",ref.group = ".all.", label.y = (max(info$total_norm_cts)+(max(info$total_norm_cts)*.05)))
fit <- t.test(total_norm_cts~kit, data=info)
summary(fit)
p2 <- ggboxplot(box_plts, "kit", "total_raw_cts",color = cols, ylab = "raw counts\n", xlab="kit")+
  scale_y_continuous(labels = scales::scientific)+
  stat_compare_means(method = "t.test", label.y = max(info$total_raw_cts)+(max(info$total_raw_cts)*.1),label.x = 2.05)+      # Add global p-value
  stat_compare_means(label = "p.signif", method = "t.test",ref.group = ".all.", label.y = (max(info$total_raw_cts)+(max(info$total_raw_cts)*.05)))
fit <- t.test(total_raw_cts~kit, data=info)
summary.lm(fit)
p3 <- ggboxplot(box_plts, "kit", "raw_fastq_reads",color = cols, ylab = "raw fastq reads\n", xlab="kit")+
  scale_y_continuous(labels = scales::scientific)+
  stat_compare_means(method = "t.test", label.y = max(info$raw_fastq_reads)+(max(info$raw_fastq_reads)*.1),label.x = 2.1)+      # Add global p-value
  stat_compare_means(label = "p.signif", method = "t.test",ref.group = ".all.", label.y = (max(info$raw_fastq_reads)+(max(info$raw_fastq_reads)*.05)))
fit <- t.test(raw_fastq_reads~kit, data=info)
summary.lm(fit)
p4 <- ggboxplot(box_plts, "kit", "tin_median",color = cols, ylab = "TIN\n", xlab="kit")+
  #scale_y_continuous(labels = scales::scientific)+
  stat_compare_means(method = "t.test", label.y = max(info$tin_median)+(max(info$tin_median)*.15),label.x = 2.1)+      # Add global p-value
  stat_compare_means(label = "p.signif", method = "t.test",ref.group = ".all.", label.y = (max(info$tin_median)+(max(info$tin_median)*.05)))
fit <- t.test(tin_median~kit, data=info)
summary.lm(fit)
p5 <- ggboxplot(box_plts, "kit", "median_gc_content",color = cols, ylab = "median GC content\n", xlab="kit")+
  #scale_y_continuous(labels = scales::scientific)+
  stat_compare_means(method = "t.test", label.y = max(info$median_gc_content)+(max(info$median_gc_content)*.02),label.x = 2.1)+      # Add global p-value
  stat_compare_means(label = "p.signif", method = "t.test",ref.group = ".all.", label.y = (max(info$median_gc_content)+(max(info$median_gc_content)*.01)))
fit <- t.test(median_gc_content~kit, data=info)
summary.lm(fit)
p6 <- ggboxplot(box_plts, "kit", "prop_duplicate_reads",color = cols, ylab = "proportion of duplicate reads\n", xlab="kit")+
  #scale_y_continuous(labels = scales::scientific)+
  stat_compare_means(method = "t.test", label.y = max(info$prop_duplicate_reads)+(max(info$prop_duplicate_reads)*.15),label.x = 2.1)+      # Add global p-value
  stat_compare_means(label = "p.signif", method = "t.test",ref.group = ".all.", label.y = (max(info$prop_duplicate_reads)+(max(info$prop_duplicate_reads)*.05)))
fit <- t.test(prop_duplicate_reads~kit, data=info)
summary.lm(fit)

tiff("D:/Martin Lab/RNA-seq/axm/post_reviews/qc_stats/plots/kit_counts.tiff", width = 11, height = 3, units = 'in', res = 1000)
grid.arrange(p1, p2,p3, ncol=3)
dev.off()
tiff("D:/Martin Lab/RNA-seq/axm/post_reviews/qc_stats/plots/kit_stats.tiff", width = 11, height = 3, units = 'in', res = 1000)
grid.arrange(p4, p5,p6, ncol=3)
dev.off()



h17 <- info[which(info$type_stage == "h 17dpf"),]
h8 <- info[which(info$type_stage == "h 8dpf"),]
p17 <- info[which(info$type_stage == "p 17dpf"),]
p8 <- info[which(info$type_stage == "p 8dpf"),]
box_names <- c("hybrids", "hybrids", "parents", "parents")
cols <- c("#0067a7","#bd1e24","#0067a7","#bd1e24")
reset.par() 
par(mfrow = c( 2, 3 ), mai= c(.5,.5,.5,.5))
boxplot(h17$total_norm_cts, h8$total_norm_cts, p17$total_norm_cts, p8$total_norm_cts,ylab = "", col = cols, names=box_names)
boxplot(h17$total_raw_cts, h8$total_raw_cts, p17$total_raw_cts, p8$total_raw_cts, col = cols, names=box_names)
boxplot(h17$raw_fastq_reads, h8$raw_fastq_reads, p17$raw_fastq_reads, p8$raw_fastq_reads, col = cols, names=box_names)

boxplot(h17$tin_median, h8$tin_median, p17$tin_median, p8$tin_median, col = cols, names=box_names)
boxplot(h17$median_gc_content, h8$median_gc_content, p17$median_gc_content, p8$median_gc_content, col = cols, names=box_names)
boxplot(h17$prop_duplicate_reads, h8$prop_duplicate_reads, p17$prop_duplicate_reads, p8$prop_duplicate_reads, col = cols, names=box_names)


fit <- aov(total_norm_cts~type_stage, data=info)
summary.lm(fit)
fit <- aov(total_raw_cts~type_stage, data=info)
summary.lm(fit)
fit <- aov(raw_fastq_reads~type_stage, data=info)
summary.lm(fit)


t.test(info$total_norm_cts~info$type_stage)

lm <- lm(info$total_norm_cts~info$tin_median+info$stage+info$sequencing_round+info$kit)
summary.lm(lm)
lm <- lm(info$total_norm_cts~info$tin_median+info$prop_duplicate_reads+info$median_gc_content+info$kit)
summary.lm(lm)


mapped <- read.csv("D:/Martin Lab/RNA-seq/axm/post_reviews/mapped_reads_to_features.csv")
mappedt <- read.csv("D:/Martin Lab/RNA-seq/axm/post_reviews/mapped_reads_to_features_transposed.csv")
cols <- 	c("#f6c700", "#FF00CC", "#bd1e24", "#0067a7","#0067a7")

head(mappedt)
mapped <- mapped[which(mapped$stage == 8),]
mapped <- mapped[which(mapped$species != "h"),]
mappedt <- t(mappedt)

counts <- table(mapped$sample, mapped$Assigned)
barplot(as.matrix(mappedt), col = cols)
axis(1, at=c(), cex.axis = 1.5)


mappedt$row <- seq_len(nrow(mappedt))
dat2 <- melt(mappedt, id.vars = "sample")

library(ggplot2)

ggplot(dat2, aes(x=variable, y=value, fill=sample)) + 
  geom_bar(stat="identity") +
  guides(fill=FALSE) +
  theme_bw()


t.test(prop_mapped_to_features ~ species, data=mapped,alternative = c("two.sided", "less", "greater"))
fit <- aov(prop_mapped_to_features ~ stage, data=mapped)
summary(fit)

cols <- 	c("#f6c700", "#E77200", "#FF00CC","#000000", "#bd1e24", "#0067a7")
barplot(mapped$prop_mapped_to_features, col = c("red3","purple3","blue3")[mapped$species],names.arg=info$sample,las=2, cex.names=.8)


collapser <- c(1)
for (collapse in collapser)
{
pdf("D:/Martin Lab/RNA-seq/axm/post_reviews/RSeq_QC.pdf",width=7,height=5)
  
par(xpd = T, mar = par()$mar + c(0,0,0,5))
lm <- lm(info$total_raw_cts~info$tin_median)
plot    (info$tin_median, info$total_raw_cts,pch = c(17, 17, 19, 19)[info$sequencing_round], col = cols_in,
         xlab = "TIN (transcript integrity number)", ylab = "raw read counts",
         main = "raw read counts")
legend("bottomright", title = c("4 sequencing rounds"),legend = c("A and M 17dpf", "MxA hybrid 17dpf", "round 3 & 4 8dpf"), inset=c(-0.3,0),bty = "n", xpd=TRUE, mar(c(7,7,7,10)), cex = 0.8, pch = c(19, 19, 17), col = c("#E77200","#0067a7","#000000"))
summary(lm)
par(xpd = F)
abline(lm)
textplot(capture.output(summary(lm)), mar = c(3,5,3,5))

par(xpd = T)
lm <- lm(info$total_norm_cts~info$tin_median)
plot    (info$tin_median, info$total_norm_cts,pch = c(17, 17, 19, 19)[info$sequencing_round], col = cols_in,
         xlab = "TIN (transcript integrity number)", ylab = "normalized read counts",
         main = "normalized read counts")
legend("bottomright", title = c("4 sequencing rounds"),legend = c("A and M 17dpf", "MxA hybrid 17dpf", "round 3 & 4 8dpf"), inset=c(-0.3,0),bty = "n", xpd=TRUE, mar(c(7,7,7,10)), cex = 0.8, pch = c(19, 19, 17), col = c("#E77200","#0067a7","#000000"))
summary(lm)
par(xpd = F)
abline(lm)
textplot(capture.output(summary(lm)), mar = c(3,5,3,5))

par(xpd = T)
lm <- lm(info$reads_mapped~info$tin_median)
plot    (info$tin_median, info$reads_mapped,pch = c(17, 17, 19, 19)[info$sequencing_round], col = cols_in,
         xlab = "TIN (transcript integrity number)", ylab = "mapped reads",
         main = "mapped reads")
legend("bottomright", title = c("4 sequencing rounds"),legend = c("A and M 17dpf", "MxA hybrid 17dpf", "round 3 & 4 8dpf"), inset=c(-0.3,0),bty = "n", xpd=TRUE, mar(c(7,7,7,10)), cex = 0.8, pch = c(19, 19, 17), col = c("#E77200","#0067a7","#000000"))
summary(lm)
par(xpd = F)
abline(lm)
textplot(capture.output(summary(lm)), mar = c(3,5,3,5))

par(xpd = T)
lm <- lm(info$median_gc_content~info$tin_median)
plot    (info$tin_median, info$median_gc_content,pch = c(17, 17, 19, 19)[info$sequencing_round], col = cols_in,
         xlab = "TIN (transcript integrity number)", ylab = "median GC content across all genes",
         main = "median GC content across all genes")
legend("bottomright", title = c("4 sequencing rounds"),legend = c("A and M 17dpf", "MxA hybrid 17dpf", "round 3 & 4 8dpf"), inset=c(-0.3,0),bty = "n", xpd=TRUE, mar(c(7,7,7,10)), cex = 0.8, pch = c(19, 19, 17), col = c("#E77200","#0067a7","#000000"))
summary(lm)
par(xpd = F)
abline(lm)
textplot(capture.output(summary(lm)), mar = c(3,5,3,5))

par(xpd = T)
lm <- lm(info$prop_duplicate_reads~info$tin_median)
plot    (info$tin_median, info$prop_duplicate_reads,pch = c(17, 17, 19, 19)[info$sequencing_round], col = cols_in,
         xlab = "TIN (transcript integrity number)", ylab = "proportion of duplicate reads",
                                                     main = "proportion of duplicate reads")
legend("bottomright", title = c("4 sequencing rounds"),legend = c("A and M 17dpf", "MxA hybrid 17dpf", "round 3 & 4 8dpf"), inset=c(-0.3,0),bty = "n", xpd=TRUE, mar(c(7,7,7,10)), cex = 0.8, pch = c(19, 19, 17), col = c("#E77200","#0067a7","#000000"))
par(xpd = F)
abline(lm)
summary(lm)
textplot(capture.output(summary(lm)), mar = c(3,5,3,5))

lm <- lm(info$total_norm_cts~info$tin_median+info$stage+info$sequencing_round+info$kit)
summary(lm)
textplot(capture.output(summary(lm)), mar = c(3,5,3,5))

dev.off()
}


plot(info$tin_median, info$total_norm_cts)
abline(lm, col = "red")

lm <- lm(info$total_norm_cts~info$stage)
summary (lm)
plot(info$stage, info$total_norm_cts)
abline(lm, col = "red")

lm <- lm(info$median_gc_content~info$tin_median)
summary(lm)
plot(info$tin_median, info$median_gc_content)
abline(lm, col = "red")

lm <- lm(info$total_raw_cts~info$tin_median+info$stage+info$sequencing_round+info$kit)
summary(lm)




#####
########################################
########################################
##### pleiotropy #######################

library(Rcpp)
library(reshape2)
library(plyr)
library(MASS)
zfin <- read.table("D:/Martin Lab/RNA-seq/pleiotropy/gene_association.jam.zfin", fill=TRUE,header = FALSE,quote = "", stringsAsFactors = FALSE, sep = "\t")
head(zfin)
zfin <- zfin[which(zfin$V9 == "P"),]
zfin <- zfin[which(zfin$V7 == "EXP" | zfin$V7 == "IDA"| zfin$V7 == "IPI" | zfin$V7 == "IMP" | zfin$V7 == "IGI" | zfin$V7 == "IEP"),]
head(zfin)

inheritance_8dpf  <- read.csv("D:/Martin Lab/RNA-seq/axm/post_reviews/inheritance_patterns/OAxOM_8dpf_inheritance.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
inheritance_17dpf <- read.csv("D:/Martin Lab/RNA-seq/axm/post_reviews/inheritance_patterns/LAxLM_17dpf_inheritance.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
inheritance <- inheritance_17dpf

over <- inheritance[which(inheritance$inheritance == 5),]
over <- unique(over$zeb_gene_symbol_one_way)
under <- inheritance[which(inheritance$inheritance == 6),]
under <- unique(under$zeb_gene_symbol_one_way)
con <- inheritance[which(inheritance$inheritance == 1),]
con <- unique(con$zeb_gene_symbol_one_way)


head(con)
nrow(over)
over_zeb <- zfin[zfin$V3 %in% over,]
over_zeb$go_index <- paste(over_zeb$V3, over_zeb$V5, sep = ";")
over_zeb <- subset(over_zeb, !duplicated(over_zeb[,18]))
over_zeb_count <- count(over_zeb, "V3")
over_zeb_count$type <- "over"
head(over_zeb_count)
nrow(over_zeb_count)
under_zeb <- zfin[zfin$V3 %in% under,]
under_zeb$go_index <- paste(under_zeb$V3, under_zeb$V5, sep = ";")
under_zeb <- subset(under_zeb, !duplicated(under_zeb[,18]))
under_zeb_count <- count(under_zeb, "V3")
under_zeb_count$type <- "under"
head(under_zeb_count)
nrow(under_zeb_count)
con_zeb <- zfin[zfin$V3 %in% con,]
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
#write.table(all_cranial_go,"D:/Martin Lab/RNA-seq/axm/post_reviews/pleiotropy/17dpf_one_way_hit_up_dn_con.txt", quote = FALSE, row.names = FALSE, sep = "\t")
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

## ppi

ppi <- read.table("D:/Martin Lab/RNA-seq/pleiotropy/danio_ppi_experimental", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
head(ppi)
range(ppi$experimental)
range(ppi$combined_score)
quantile(ppi$experimental,.95)
quantile(ppi$combined_score,.95)
ppi_sig <- ppi[which(ppi$combined_score >= 481 & ppi$experimental >= 265),]
#ppi_sig <- ppi
ppi_counts <-  count(ppi_sig, 'protein1')
head(ppi_counts)

ensemble <- read.table("D:/Martin Lab/RNA-seq/axm/post_reviews/pleiotropy/ensemble_ids_one_way_17dpf.txt",header = TRUE,quote = "", stringsAsFactors = FALSE, sep = "\t")
zeb <- read.table("D:/Martin Lab/RNA-seq/axm/post_reviews/pleiotropy/17dpf_one_way_hit_up_dn_con.txt",header = TRUE,quote = "", stringsAsFactors = FALSE, sep = "\t")
ensemble <- merge(ensemble, zeb, by = "V3")

head(ensemble)
ensemble$tag <- paste(ensemble$V3, ensemble$Ensembl.Protein.ID, sep= ":")
ensemble$protein1 <- ensemble$Ensembl.Protein.ID
all_ppi <- merge(ensemble, ppi_counts, by = c("protein1"))
head(all_ppi)
nrow(all_ppi)
unique(all_ppi$type)
#all_ppi$type[all_ppi$type == "con"]
#length(unique(ppi_genes$protein1))
#over$Gene <- over$V1
#under$Gene <- under$V1
#con$Gene <- con$V1
#over_ppi <- merge(over, ppi_genes, by = ("Gene"))
#over_ppi$type <- "over"
#under_ppi <- merge(under, ppi_genes, by = ("Gene"))
#under_ppi$type <- "under"
#con_ppi <- merge(con, ppi_genes, by = ("Gene"))
#con_ppi$type <- "con"
#con_ppi$type2 <- "con"
#under_ppi$type2 <- "mis"
#over_ppi$type2 <- "mis"
#nrow(con_ppi)
#nrow(over_ppi)
#nrow(under_ppi)
##under_ppi <- under_ppi[which(under_ppi$freq < 1000),]
#all_ppi <- rbind(con_ppi, under_ppi, over_ppi)


#tiff("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/axm/pleio_all_ppi.tif", width = 5.5, height = 6, units = 'in', res = 600)
p <- ggplot(all_ppi, aes(x=type, y=freq.y, fill = type)) + 
  #geom_boxplot()
  geom_violin(trim = TRUE, bw = 10)#+geom_dotplot(binaxis='y', stackdir='center',
#position=position_dodge(1),dotsize = .3)
p+ scale_fill_manual(values=c("#f6c700", "#bd1e24", "#0067a7")) + theme(axis.text.y = element_text(size=18, color = '#000000'),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                        panel.background = element_blank(),axis.line = element_line(colour = "black"),
                                                                        axis.text.x = element_text(size=18, color = '#000000'), legend.position="none",
                                                                        axis.title.y=element_blank(),axis.title.x=element_blank()) + stat_summary(fun.data=data_summary, col = "black")
#dev.off()
median(con_ppi$freq.y)
median(under_ppi$freq.y)
median(over_ppi$freq.y)
mean(con_ppi$freq.y)
mean(under_ppi$freq.y)
mean(over_ppi$freq.y)

ppi_lm <- glm(all_ppi$freq.y~all_ppi$type2, family = poisson)
ppi_lm <- glm.nb(all_ppi$freq.y~all_ppi$type)
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



#####
########################################
########################################
##### transcript integrity numbers #####
########################################

setwd("D:/Martin Lab/RNA-seq/axm/post_reviews/tins/")
samples <- read.table("D:/Martin Lab/RNA-seq/axm/post_reviews/table_maker.txt", header = TRUE, stringsAsFactors = FALSE)
tin_files <- samples$sample
#tin_file <- "axmJ1"

sample_names <- c()
tin_meds <- c()

for (tin_file in tin_files)
{
  
tin_in <- read.table(paste(tin_file, ".rna.sort.summary.txt", sep = ""), header = TRUE, stringsAsFactors = FALSE) 
sample_names <- c(tin_file, sample_names)
tin_meds <- c(tin_in$TIN.median., tin_meds)

}

tins <- data.frame(sample_names,tin_meds)
names(tins) <- c("sample","tin_median")
final_samples <- merge(samples, tins, by = ("sample"))
write.table(final_samples, "D:/Martin Lab/RNA-seq/axm/post_reviews/table_maker_tins.txt", sep = "\t", quote = FALSE, row.names = FALSE)



#####
########################################
########################################
##### danio genes overlap ##############
########################################
library(reshape2)
a <- read.delim("C:/Users/jmcgirr/Desktop/Danio_rerio.GRCz11.96.gtf", header = FALSE, sep = ";")
head(a)
b <- as.data.frame(a$V3)
b$biotype <- a$V5
b <- cbind(b, colsplit(b$biotype, "gene_biotype ", c("j", "type")))
b <- cbind(b, colsplit(b[,1], "gene_name ", c("j2", "gene")))
head(b)
b <- b[which(b$type == "protein_coding"),]
c <- unique(b$gene)
length(c)
zeb <- tolower(c)
pup <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/ase_data/genes.saf", header = TRUE, stringsAsFactors = FALSE)
head(pup)
nrow(pup)
pup<- cbind(pup, colsplit(pup$GeneID, ";", c("j2", "gene")))
pup <- tolower(pup$gene)
length(intersect(pup,zeb))



#####
########################################
########################################
##### MBASED ASE #######################
########################################

## all h no p ase
##h_jaw_ase 
#ase_files<- c("axmJ1_mbased_ase_unphased.txt","axmJ4_mbased_ase_unphased.txt","axmJ5_mbased_ase_unphased.txt","axmJ6_mbased_ase_unphased.txt")
#ase <- read.table(paste("D:/Martin Lab/RNA-seq/axm/plosOne/",ase_files[1], sep = ""), header = TRUE, stringsAsFactors = FALSE)
#sig <- ase[which(ase$pValueASE <= 0.05),]
#sig_ase <-sig$geneID
#for (i in c(2:length(ase_files)))
#{
#ase <- read.table(paste("D:/Martin Lab/RNA-seq/axm/plosOne/",ase_files[i], sep = ""), header = TRUE, stringsAsFactors = FALSE)
#sig <- ase[which(ase$pValueASE <= 0.05),]
#sig_ase <-intersect(sig_ase, sig$geneID) 
#
#}
##p_jaw_ase 
#ase_files<- c("CPAJ1_mbased_ase_unphased.txt","CPAJ2_mbased_ase_unphased.txt","CPAJ3_mbased_ase_unphased.txt","CPMJ1_mbased_ase_unphased.txt","CPMJ2_mbased_ase_unphased.txt","CPMJ3_mbased_ase_unphased.txt","LLAJ1_mbased_ase_unphased.txt","LLAJ2_mbased_ase_unphased.txt","LLAJ3_mbased_ase_unphased.txt","LLMJ1_mbased_ase_unphased.txt","LLMJ2_mbased_ase_unphased.txt","LLMJ3_mbased_ase_unphased.txt")
#sig_ase_p <-c()
#for (ase_file in ase_files)
#{
#  ase <- read.table(paste("D:/Martin Lab/RNA-seq/axm/plosOne/",ase_file, sep = ""), header = TRUE, stringsAsFactors = FALSE)
#  sig <- ase[which(ase$pValueASE <= 0.05),]
#  sig_ase_p <-c(sig_ase_p,ase$geneID) 
#  
#}
#
#length(setdiff(sig_ase, sig_ase_p))
#
##h_emb2_ase 
#ase_files<- c("OUE1_mbased_ase_unphased.txt","OUE3_mbased_ase_unphased.txt","OUE4_mbased_ase_unphased.txt")
#ase <- read.table(paste("D:/Martin Lab/RNA-seq/axm/plosOne/",ase_files[1], sep = ""), header = TRUE, stringsAsFactors = FALSE)
#sig <- ase[which(ase$pValueASE <= 0.05),]
#sig_ase <-sig$geneID
#for (i in c(2:length(ase_files)))
#{
#  ase <- read.table(paste("D:/Martin Lab/RNA-seq/axm/plosOne/",ase_files[i], sep = ""), header = TRUE, stringsAsFactors = FALSE)
#  sig <- ase[which(ase$pValueASE <= 0.05),]
#  sig_ase <-intersect(sig_ase, sig$geneID) 
#  
#}
##p_emb_ase  
##ase_files<- c("CPAE1_mbased_ase_unphased.txt","CPAE2_mbased_ase_unphased.txt","CPAE3_mbased_ase_unphased.txt","CPME1_mbased_ase_unphased.txt","CPME2_mbased_ase_unphased.txt","CPME3_mbased_ase_unphased.txt","LLAE1_mbased_ase_unphased.txt","LLAE2_mbased_ase_unphased.txt","LLAE3_mbased_ase_unphased.txt","LLME1_mbased_ase_unphased.txt","LLME2_mbased_ase_unphased.txt","LLME3_mbased_ase_unphased.txt")
##p_emb2_ase 
#ase_files<- c("CAE1_mbased_ase_unphased.txt","CAE2_mbased_ase_unphased.txt","CAE3_mbased_ase_unphased.txt","CME1_mbased_ase_unphased.txt","CME2_mbased_ase_unphased.txt","CME5_mbased_ase_unphased.txt","OAE2_mbased_ase_unphased.txt","OAE3_mbased_ase_unphased.txt","OAE4_mbased_ase_unphased.txt","OME1_mbased_ase_unphased.txt","OME2_mbased_ase_unphased.txt","OME3_mbased_ase_unphased.txt")
#sig_ase_p <-c()
#for (ase_file in ase_files)
#{
#  ase <- read.table(paste("D:/Martin Lab/RNA-seq/axm/plosOne/",ase_file, sep = ""), header = TRUE, stringsAsFactors = FALSE)
#  sig <- ase[which(ase$pValueASE <= 0.05),]
#  sig_ase_p <-c(sig_ase_p,ase$geneID) 
#  
#}
#
#length(setdiff(sig_ase, sig_ase_p))


# all h no all p

#h_jaw_ase 
ase_files<- c("axmJ1_mbased_ase_unphased.txt","axmJ4_mbased_ase_unphased.txt","axmJ5_mbased_ase_unphased.txt","axmJ6_mbased_ase_unphased.txt")
ase <- read.table(paste("D:/Martin Lab/RNA-seq/axm/plosOne/",ase_files[1], sep = ""), header = TRUE, stringsAsFactors = FALSE)
ase <- cbind(ase, colsplit(ase$geneID, ":", c("gene", "other")))
sig <- ase[which(ase$pValueASE <= 0.05),]
sig_ase <-sig$gene
for (i in c(2:length(ase_files)))
{
  ase <- read.table(paste("D:/Martin Lab/RNA-seq/axm/plosOne/",ase_files[i], sep = ""), header = TRUE, stringsAsFactors = FALSE)
  ase <- cbind(ase, colsplit(ase$geneID, ":", c("gene", "other")))
  sig <- ase[which(ase$pValueASE <= 0.05),]
  sig_ase <-intersect(sig_ase, sig$gene) 
  
}
#p_jaw_ase 
ase_files<- c("CPAJ1_mbased_ase_unphased.txt","CPAJ2_mbased_ase_unphased.txt","CPAJ3_mbased_ase_unphased.txt","CPMJ1_mbased_ase_unphased.txt","CPMJ2_mbased_ase_unphased.txt","CPMJ3_mbased_ase_unphased.txt","LLAJ1_mbased_ase_unphased.txt","LLAJ2_mbased_ase_unphased.txt","LLAJ3_mbased_ase_unphased.txt","LLMJ1_mbased_ase_unphased.txt","LLMJ2_mbased_ase_unphased.txt","LLMJ3_mbased_ase_unphased.txt")
ase <- read.table(paste("D:/Martin Lab/RNA-seq/axm/plosOne/",ase_files[1], sep = ""), header = TRUE, stringsAsFactors = FALSE)
ase <- cbind(ase, colsplit(ase$geneID, ":", c("gene", "other")))
sig <- ase[which(ase$pValueASE <= 0.05),]
sig_ase_p <-sig$gene
for (i in c(2:length(ase_files)))
{
  ase <- read.table(paste("D:/Martin Lab/RNA-seq/axm/plosOne/",ase_files[i], sep = ""), header = TRUE, stringsAsFactors = FALSE)
  ase <- cbind(ase, colsplit(ase$geneID, ":", c("gene", "other")))
  sig <- ase[which(ase$pValueASE <= 0.05),]
  sig_ase_p <-intersect(sig_ase_p, sig$gene) 
  
}

length(setdiff(sig_ase, sig_ase_p))

#h_emb2_ase 
ase_files<- c("OUE1_mbased_ase_unphased.txt","OUE3_mbased_ase_unphased.txt","OUE4_mbased_ase_unphased.txt")
ase <- read.table(paste("D:/Martin Lab/RNA-seq/axm/plosOne/",ase_files[1], sep = ""), header = TRUE, stringsAsFactors = FALSE)
ase <- cbind(ase, colsplit(ase$geneID, ":", c("gene", "other")))
sig <- ase[which(ase$pValueASE <= 0.05),]
sig_ase <-sig$gene
for (i in c(2:length(ase_files)))
{
  ase <- read.table(paste("D:/Martin Lab/RNA-seq/axm/plosOne/",ase_files[i], sep = ""), header = TRUE, stringsAsFactors = FALSE)
  ase <- cbind(ase, colsplit(ase$geneID, ":", c("gene", "other")))
  sig <- ase[which(ase$pValueASE <= 0.05),]
  sig_ase <-intersect(sig_ase, sig$gene) 
  
}
#p_emb_ase  
#ase_files<- c("CPAE1_mbased_ase_unphased.txt","CPAE2_mbased_ase_unphased.txt","CPAE3_mbased_ase_unphased.txt","CPME1_mbased_ase_unphased.txt","CPME2_mbased_ase_unphased.txt","CPME3_mbased_ase_unphased.txt","LLAE1_mbased_ase_unphased.txt","LLAE2_mbased_ase_unphased.txt","LLAE3_mbased_ase_unphased.txt","LLME1_mbased_ase_unphased.txt","LLME2_mbased_ase_unphased.txt","LLME3_mbased_ase_unphased.txt")
#p_emb2_ase 
ase_files<- c("CAE1_mbased_ase_unphased.txt","CAE2_mbased_ase_unphased.txt","CAE3_mbased_ase_unphased.txt","CME1_mbased_ase_unphased.txt","CME2_mbased_ase_unphased.txt","CME5_mbased_ase_unphased.txt","OAE2_mbased_ase_unphased.txt","OAE3_mbased_ase_unphased.txt","OAE4_mbased_ase_unphased.txt","OME1_mbased_ase_unphased.txt","OME2_mbased_ase_unphased.txt","OME3_mbased_ase_unphased.txt")
ase <- read.table(paste("D:/Martin Lab/RNA-seq/axm/plosOne/",ase_files[1], sep = ""), header = TRUE, stringsAsFactors = FALSE)
ase <- cbind(ase, colsplit(ase$geneID, ":", c("gene", "other")))
sig <- ase[which(ase$pValueASE <= 0.05),]
sig_ase_p <-sig$gene
for (i in c(2:length(ase_files)))
{
  ase <- read.table(paste("D:/Martin Lab/RNA-seq/axm/plosOne/",ase_files[i], sep = ""), header = TRUE, stringsAsFactors = FALSE)
  sig <- ase[which(ase$pValueASE <= 0.05),]
  ase <- cbind(ase, colsplit(ase$geneID, ":", c("gene", "other")))
  sig_ase_p <-intersect(sig_ase_p, sig$gene) 
  
}

length(setdiff(sig_ase, sig_ase_p))

mis <- read.csv("D:/Martin Lab/RNA-seq/axm/post_reviews/conditions/DE_(CRPA1_&_LILA_&_CRPM1_&_LILM)_vs_(LAxLM)_17dpf_genes.csv", header = TRUE, stringsAsFactors = FALSE)
am <- read.csv("D:/Martin Lab/RNA-seq/axm/post_reviews/conditions/DE_(CRPA1)_vs_(LILA)_17dpf_genes.csv", header = TRUE, stringsAsFactors = FALSE)
mis <- read.csv("D:/Martin Lab/RNA-seq/axm/post_reviews/conditions/DE_(OSPA_&_OSPM)_vs_(OAxOM)_8dpf_genes.csv", header = TRUE, stringsAsFactors = FALSE)
am <- read.csv("D:/Martin Lab/RNA-seq/axm/post_reviews/conditions/DE_(OSPA)_vs_(OSPM)_8dpf_genes.csv", header = TRUE, stringsAsFactors = FALSE)

am <- am[which(am$padj >0.05),]
mis <- mis[which(mis$padj <0.05),]
comp <- intersect(am$symbol,setdiff(sig_ase, sig_ase_p))
comp_mis <- intersect(comp, mis$symbol)
length(comp)
length(comp_mis)

inheritance_and_ase_8dpf  <- read.csv("D:/Martin Lab/RNA-seq/axm/post_reviews/inheritance_patterns/inheritance_and_ase_8.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
inheritance_and_ase_17dpf  <- read.csv("D:/Martin Lab/RNA-seq/axm/post_reviews/inheritance_patterns/inheritance_and_ase_17.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
comp <- inheritance_and_ase_17dpf[which(inheritance_and_ase_17dpf$ase_in_hybrids_not_parents == "yes" &
                                        inheritance_and_ase_17dpf$padj_p1p2 > 0.05),]
nrow(comp)
comp_mis <- intersect(comp$symbol, mis$symbol)
length(comp_mis)


#####













###################### old #####

#condition_species_ap_embryo.txt
#condition_species_ap_jaw.txt
#condition_species_apxm_embryo.txt
#condition_species_apxm_jaw.txt


col_data <-                                                          "condition_species_ham_jaw.txt"
cts_data <-                                                               "DESeq_counts_ham_jaw.txt"
cts_data_genes <-                                                                   "DE_ham_jaw.csv"
genes_out <-                                                                        "DE_ham_jaw_genes.csv"
comparison <-                                                                          "am x h"

sample_list <- read.table(col_data, header = TRUE, stringsAsFactors = FALSE)
head(sample_list)
keeps <- c("Geneid", sample_list$sample)
keeper <- sample_list$sample

cts <- cts[keeps]
head(cts)

setwd("D:/Martin Lab/RNA-seq/axm/revisions/")
#write.table(cts, cts_data, row.names = FALSE, quote= FALSE,sep="\t") 



#setwd("D:/Martin Lab/RNA-seq/axm/post_reviews/ASE")
#axmJ1 <- read.csv("axmJ1_gene_haplotype_counts.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
#axmJ4 <- read.csv("axmJ4_gene_haplotype_counts.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
#axmJ5 <- read.csv("axmJ5_gene_haplotype_counts.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
#axmJ6 <- read.csv("axmJ6_gene_haplotype_counts.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
#OUE1 <-  read.csv("OUE1_gene_haplotype_counts.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
#OUE3 <-  read.csv("OUE3_gene_haplotype_counts.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
#OUE4 <-  read.csv("OUE4_gene_haplotype_counts.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
#
#axmJ1 <- read.csv("D:/Martin Lab/RNA-seq/axm/post_reviews/ASE/min_cov_20/axmJ1_binomial.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
#axmJ4 <- read.csv("D:/Martin Lab/RNA-seq/axm/post_reviews/ASE/min_cov_20/axmJ4_binomial.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
#axmJ5 <- read.csv("D:/Martin Lab/RNA-seq/axm/post_reviews/ASE/min_cov_20/axmJ5_binomial.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
#axmJ6 <- read.csv("D:/Martin Lab/RNA-seq/axm/post_reviews/ASE/min_cov_20/axmJ6_binomial.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
#axmJ1$snp_tag <- paste(axmJ1$geneIndex, axmJ1$snpIndex, sep = ":")
#axmJ4$snp_tag <- paste(axmJ4$geneIndex, axmJ4$snpIndex, sep = ":")
#axmJ5$snp_tag <- paste(axmJ5$geneIndex, axmJ5$snpIndex, sep = ":")
#axmJ6$snp_tag <- paste(axmJ6$geneIndex, axmJ6$snpIndex, sep = ":")
#axmJ1$totalCount <- axmJ1$alleleOneCount + axmJ1$alleleTwoCount
#axmJ4$totalCount <- axmJ4$alleleOneCount + axmJ4$alleleTwoCount
#axmJ5$totalCount <- axmJ5$alleleOneCount + axmJ5$alleleTwoCount
#axmJ6$totalCount <- axmJ6$alleleOneCount + axmJ6$alleleTwoCount
#axmJ1 <- axmJ1[which(axmJ1$totalCount >= 20 & axmJ1$alleleOneCount > 0 & axmJ1$alleleTwoCount > 0),]
#axmJ4 <- axmJ4[which(axmJ4$totalCount >= 20 & axmJ4$alleleOneCount > 0 & axmJ4$alleleTwoCount > 0),]
#axmJ5 <- axmJ5[which(axmJ5$totalCount >= 20 & axmJ5$alleleOneCount > 0 & axmJ5$alleleTwoCount > 0),]
#axmJ6 <- axmJ6[which(axmJ6$totalCount >= 20 & axmJ6$alleleOneCount > 0 & axmJ6$alleleTwoCount > 0),]
#OUE1 <- read.csv("D:/Martin Lab/RNA-seq/axm/post_reviews/ASE/min_cov_20/OUE1_binomial.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
#OUE3 <- read.csv("D:/Martin Lab/RNA-seq/axm/post_reviews/ASE/min_cov_20/OUE3_binomial.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
#OUE4 <- read.csv("D:/Martin Lab/RNA-seq/axm/post_reviews/ASE/min_cov_20/OUE4_binomial.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
#OUE1$snp_tag <- paste(OUE1$geneIndex, OUE1$snpIndex, sep = ":")
#OUE3$snp_tag <- paste(OUE3$geneIndex, OUE3$snpIndex, sep = ":")
#OUE4$snp_tag <- paste(OUE4$geneIndex, OUE4$snpIndex, sep = ":")
#OUE1$totalCount <- OUE1$alleleOneCount + OUE1$alleleTwoCount
#OUE3$totalCount <- OUE3$alleleOneCount + OUE3$alleleTwoCount
#OUE4$totalCount <- OUE4$alleleOneCount + OUE4$alleleTwoCount
#OUE1 <- OUE1[which(OUE1$totalCount >= 20 & OUE1$alleleOneCount > 0 & OUE1$alleleTwoCount > 0),]
#OUE3 <- OUE3[which(OUE3$totalCount >= 20 & OUE3$alleleOneCount > 0 & OUE3$alleleTwoCount > 0),]
#OUE4 <- OUE4[which(OUE4$totalCount >= 20 & OUE4$alleleOneCount > 0 & OUE4$alleleTwoCount > 0),]
#
#
#range(axmJ6$alleleTwoCount)
#head(axmJ6)
#a <- merge(axmJ1, axmJ4, by = c("snp_tag"))
#a <- merge(a, axmJ5, by = c("snp_tag"))
#a <- merge(a, axmJ6, by = c("snp_tag"))
#head(a)
#nrow(a)
#a <- merge(axmJ1, axmJ4, by = c("geneIndex"))
#a <- merge(a, axmJ5, by = c("geneIndex"))
#a <- merge(a, axmJ6, by = c("geneIndex"))
#head(a)
#nrow(a)
#range(a$binom_test_p_value.x)
#
#a <- merge(OUE1, OUE3, by = c("snp_tag"))
#a <- merge(a, OUE4, by = c("snp_tag"))
#head(a)
#nrow(a)
#a <- merge(OUE1, OUE3, by = c("geneIndex"))
#a <- merge(a, OUE4, by = c("geneIndex"))
#head(a)
#nrow(a)
#range(a$binom_test_p_value.x)
#
#
#count_genes_with_het_snps <-  count(a, 'geneIndex')
#head(count_genes_with_het_snps)
#hist(count_genes_with_het_snps$freq)#, ylim = c(0,250))
## # of genes with 1 snp showing ase
#nrow(count_genes_with_het_snps[which(count_genes_with_het_snps$freq == 1),])
## # of genes with more than 1 snp showing ase
#nrow(count_genes_with_het_snps[which(count_genes_with_het_snps$freq > 1),])
#hybrid_ase_genes <- unique(a$geneIndex.x)
#length(unique(a$geneIndex.x))
## # of genes showing ase in hybrids but not parents
#length(setdiff(hybrid_ase_genes, parental_ase_genes_8))

#### DIFFERENTIAL EXPRESSION ANALYSIS ####
#vignette("DESeq2")

cts <- as.matrix(read.table(cts_data ,sep = "\t",header = TRUE,row.names=1))
cts <- na.omit(cts)
head(cts)
nrow(cts)

setwd("D:/Martin Lab/RNA-seq/axm/conditions")
colData <- as.matrix(read.table(col_data ,header = TRUE,row.names=1))
head(colData)


ncol(cts)
nrow(colData)
comp <- "species"
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = colData,
                              design= ~ species)

dds <- dds[ rowSums(counts(dds)) > 1, ]
dds <- DESeq(dds)

#res <- results(dds, independentFiltering=FALSE)
#res$pvalue[res$baseMean < 10] <- NA
#res$padj <- p.adjust(res$pvalue, method="BH")

res <- results(dds, alpha=0.05)
resOrdered <- res[order(res$padj),]
summary(res)

plotMA(res, ylim=c(-2,2), main = comparison)

resLFC <- lfcShrink(dds, coef=2, res=res)
plotMA(resLFC, ylim=c(-2,2), main = comparison)

plotCounts(dds, gene=which.min(res$padj), intgroup= comp)


res_ordered <- as.data.frame(resOrdered)
head(res_ordered)
#write.csv(res_ordered, file= cts_data_genes)

#### DE genes ####

cts <- read.table("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/post_qc/mapping_correction/counts_geneid", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
de_genes <- read.csv(cts_data_genes, header = TRUE, stringsAsFactors = FALSE)
head(cts)
length(unique(cts$Geneid))
head(de_genes)
de_genes$Geneid <- de_genes$X
cts$Chr <-  gsub(";.*", "", cts$Chr)
cts$Start <-  gsub(";.*", "", cts$Start)
cts$End <-  gsub(".*;", "", cts$End)
cts$Strand <-  gsub(";.*", "", cts$Strand)
keeps <- c("Geneid","Chr","Start","End","Strand")
genes <- cts[keeps]
head(genes)
de_genes["Geneid"] <- de_genes[,1]
head(de_genes)

de_genes <- merge(de_genes,genes, by = c("Geneid"))
head(de_genes)
features <- read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/feature_table.csv", header = TRUE, stringsAsFactors = FALSE)
features["Chr"] <- features$genomic_accession
features["Start"] <- features$start
features["End"] <- features$end
head(features)
unique(features$feature)
final <- merge(de_genes, features, by = c("Chr","Start", "End"))
final <- final[order(final$padj, decreasing = FALSE),]
head(final)
nrow(final)
#write.csv(final, genes_out)  


####  DE OUTLIER OVERLAP WITH DIVERGENT REGIONS AND INDELS ####
amxp_embryo <- read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/mapping_correction/DE_amxp_embryo_genes.csv", header = TRUE, stringsAsFactors = FALSE)
apxm_embryo <- read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/mapping_correction/DE_apxm_embryo_genes.csv", header = TRUE, stringsAsFactors = FALSE)
mp_embryo   <- read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/mapping_correction/DE_mp_embryo_genes.csv", header = TRUE, stringsAsFactors = FALSE)
am_embryo   <- read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/mapping_correction/DE_am_embryo_genes.csv", header = TRUE, stringsAsFactors = FALSE)
ap_embryo   <- read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/mapping_correction/DE_ap_embryo_genes.csv", header = TRUE, stringsAsFactors = FALSE)
amxp_jaw    <- read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/mapping_correction/DE_amxp_jaw_genes.csv", header = TRUE, stringsAsFactors = FALSE)
apxm_jaw    <- read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/mapping_correction/DE_apxm_jaw_genes.csv", header = TRUE, stringsAsFactors = FALSE)
mp_jaw      <- read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/mapping_correction/DE_mp_jaw_genes.csv", header = TRUE, stringsAsFactors = FALSE)
ap_jaw      <- read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/mapping_correction/DE_ap_jaw_genes.csv", header = TRUE, stringsAsFactors = FALSE)
am_jaw      <- read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/mapping_correction/DE_am_jaw_genes.csv", header = TRUE, stringsAsFactors = FALSE)



amxp_embryo$comp <- "amxp_embryo"
apxm_embryo$comp <- "apxm_embryo"
mp_embryo$comp <- "mp_embryo"
am_embryo$comp <- "am_embryo"
ap_embryo$comp <- "ap_embryo"
amxp_jaw$comp <- "amxp_jaw"
apxm_jaw$comp <- "apxm_jaw"
mp_jaw$comp <- "mp_jaw"
ap_jaw$comp <- "ap_jaw"
am_jaw$comp <- "am_jaw"

new <- rbind(amxp_embryo, 
             apxm_embryo, 
             mp_embryo, 
             am_embryo, 
             ap_embryo, 
             amxp_jaw, 
             apxm_jaw, 
             mp_jaw, 
             ap_jaw, 
             am_jaw) 

length(new)
tail(new)
head(new)

#write.csv(new, "C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/post_qc/mapping_correction/DESeq_genes_all.csv") 
sig_transcripts <- new[which(new$padj <= 0.05),]
#write.csv(sig_transcripts, "C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/post_qc/mapping_correction/DESeq_genes_sig.csv") 

amxp_embryo <- amxp_embryo[which(amxp_embryo$padj <= 0.05),]
apxm_embryo <- apxm_embryo[which(apxm_embryo$padj <= 0.05),]
mp_embryo   <- mp_embryo[which(mp_embryo$padj <= 0.05),]
am_embryo   <- am_embryo[which(am_embryo$padj <= 0.05),]
ap_embryo   <- ap_embryo[which(ap_embryo$padj <= 0.05),]
amxp_jaw    <- amxp_jaw[which(amxp_jaw$padj <= 0.05),]
apxm_jaw    <- apxm_jaw[which(apxm_jaw$padj <= 0.05),]
mp_jaw      <- mp_jaw[which(mp_jaw$padj <= 0.05),]
ap_jaw      <- ap_jaw[which(ap_jaw$padj <= 0.05),]
am_jaw      <- am_jaw[which(am_jaw$padj <= 0.05),]

head(amxp_embryo)
nrow(amxp_embryo)
length(unique(amxp_embryo$Geneid))
length(unique(apxm_embryo$Geneid))
length(unique(mp_embryo$Geneid))
length(unique(am_embryo$Geneid))
length(unique(ap_embryo$Geneid))
length(unique(amxp_jaw$Geneid))
length(unique(apxm_jaw$Geneid))
length(unique(mp_jaw$Geneid))
length(unique(ap_jaw$Geneid))
length(unique(am_jaw$Geneid))

amxp_embryo <- amxp_embryo[which(amxp_embryo$feature != "gene" & amxp_embryo$feature != "CDS"),]
apxm_embryo <- apxm_embryo[which(apxm_embryo$feature != "gene" & apxm_embryo$feature != "CDS"),]
mp_embryo   <- mp_embryo[which(mp_embryo$feature != "gene" &    mp_embryo$feature != "CDS"),]
am_embryo   <- am_embryo[which(am_embryo$feature != "gene" &    am_embryo$feature != "CDS"),]
ap_embryo   <- ap_embryo[which(ap_embryo$feature != "gene" &    ap_embryo$feature != "CDS"),]
amxp_jaw    <- amxp_jaw[which(amxp_jaw$feature != "gene" &      amxp_jaw$feature != "CDS"),]
apxm_jaw    <- apxm_jaw[which(apxm_jaw$feature != "gene" &      apxm_jaw$feature != "CDS"),]
mp_jaw      <- mp_jaw[which(mp_jaw$feature != "gene" &         mp_jaw$feature != "CDS"),]
ap_jaw      <- ap_jaw[which(ap_jaw$feature != "gene" &         ap_jaw$feature != "CDS"),]
am_jaw      <- am_jaw[which(am_jaw$feature != "gene" &         am_jaw$feature != "CDS"),]

nrow(amxp_embryo )
nrow(apxm_embryo )
nrow(mp_embryo   )
nrow(am_embryo   )
nrow(ap_embryo   )
nrow(amxp_jaw    )
nrow(apxm_jaw    )
nrow(mp_jaw      )
nrow(ap_jaw      )
nrow(am_jaw      )
length(unique(ap_jaw$related_accession))
length(unique(am_jaw$related_accession))
length(unique(am_embryo$related_accession))
length(unique(ap_embryo$related_accession))
length(unique(ap_jaw$Geneid))
length(unique(am_jaw$Geneid))
length(unique(am_embryo$Geneid))
length(unique(ap_embryo$Geneid))
head(ap_jaw)
#write.csv(am_embryo, "C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/blast/am_embryo_sig.csv") 
#write.csv(am_jaw, "C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/blast/am_jaw_sig.csv") 
#write.csv(ap_embryo, "C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/blast/ap_embryo_sig.csv") 
#write.csv(ap_jaw, "C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/blast/ap_jaw_sig.csv") 
#write.csv(mp_jaw, "C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/blast/mp_jaw_sig.csv") 
#write.csv(mp_embryo, "C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/blast/mp_embryo_sig.csv") 


m_and_pjaw <- merge(am_jaw, ap_jaw, by = c("related_accession"))
length(m_and_pjaw$related_accession)
m_and_pjaw$protein_accession <- m_and_pjaw$related_accession
keeps <- c("Geneid.x", "Chr.x", "log2FoldChange.x", "padj.x", "comp.x","protein_accession", "log2FoldChange.y", "padj.y","comp.y", "symbol.y", "name.y")
m_and_pjaw <- m_and_pjaw[keeps]
head(m_and_pjaw)
nrow(m_and_pjaw)
length(unique(m_and_pjaw$protein_accession))
length(m_and_pjaw$protein_accession)
#write.csv(m_and_pjaw, "C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/blast/shared_jaw.csv") 
#write.csv(m_and_pembryo, "C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/blast/shared_embryo.csv")
m_and_pembryo <- merge(am_embryo, ap_embryo, by = c("related_accession"))
m_and_pembryo$protein_accession <- m_and_pembryo$related_accession
keeps <- c("Geneid.x", "Chr.x", "log2FoldChange.x", "padj.x", "comp.x","protein_accession", "log2FoldChange.y", "padj.y","comp.y", "symbol.y", "name.y")
m_and_pembryo <- m_and_pembryo[keeps]
head(m_and_pembryo)
nrow(m_and_pembryo)
length(unique(m_and_pembryo$protein_accession))
#write.csv(m_and_pjaw, "C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/blast/shared_jaw.csv") 
#write.csv(m_and_pembryo, "C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/blast/shared_embryo.csv") 


m_and_pjaw1 <- m_and_pjaw[which(m_and_pjaw$log2FoldChange.x > 0 & m_and_pjaw$log2FoldChange.y > 0),] ## up in both
length(unique(m_and_pjaw1$Geneid.x))
m_and_pjaw1 <- m_and_pjaw[which(m_and_pjaw$log2FoldChange.x < 0 & m_and_pjaw$log2FoldChange.y < 0),] ## down in both
length(unique(m_and_pjaw1$Geneid.x))
m_and_pjaw1 <- m_and_pjaw[which(m_and_pjaw$log2FoldChange.x > 0 & m_and_pjaw$log2FoldChange.y < 0),] ## up in m down in p
length(unique(m_and_pjaw1$Geneid.x))
m_and_pjaw1 <- m_and_pjaw[which(m_and_pjaw$log2FoldChange.x < 0 & m_and_pjaw$log2FoldChange.y > 0),] ## up in p down in m
length(unique(m_and_pjaw1$Geneid.x))

m_and_pjaw1 <- m_and_pjaw[which(m_and_pjaw$log2FoldChange.x > 0 & m_and_pjaw$log2FoldChange.y > 0),] ## up in both
length(unique(m_and_pjaw1$protein_accession))
m_and_pjaw1 <- m_and_pjaw[which(m_and_pjaw$log2FoldChange.x < 0 & m_and_pjaw$log2FoldChange.y < 0),] ## down in both
length(unique(m_and_pjaw1$protein_accession))
m_and_pjaw1 <- m_and_pjaw[which(m_and_pjaw$log2FoldChange.x > 0 & m_and_pjaw$log2FoldChange.y < 0),] ## up in m down in p
length(unique(m_and_pjaw1$protein_accession))
m_and_pjaw1 <- m_and_pjaw[which(m_and_pjaw$log2FoldChange.x < 0 & m_and_pjaw$log2FoldChange.y > 0),] ## up in p down in m
length(unique(m_and_pjaw1$protein_accession))

m_and_pembryo1 <- m_and_pembryo[which(m_and_pembryo$log2FoldChange.x > 0 & m_and_pembryo$log2FoldChange.y > 0),] ## up in both
length(unique(m_and_pembryo1$protein_accession))
m_and_pembryo2 <- m_and_pembryo[which(m_and_pembryo$log2FoldChange.x < 0 & m_and_pembryo$log2FoldChange.y < 0),] ## down in both
length(unique(m_and_pembryo2$protein_accession))
m_and_pembryo3 <- m_and_pembryo[which(m_and_pembryo$log2FoldChange.x > 0 & m_and_pembryo$log2FoldChange.y < 0),] ## up in m down in p
length(unique(m_and_pembryo3$protein_accession))
m_and_pembryo4 <- m_and_pembryo[which(m_and_pembryo$log2FoldChange.x < 0 & m_and_pembryo$log2FoldChange.y > 0),] ## up in p down in m
length(unique(m_and_pembryo4$protein_accession))


m_and_pembryo1 <- m_and_pembryo[which(m_and_pembryo$log2FoldChange.x > 0 & m_and_pembryo$log2FoldChange.y > 0),] ## up in both
length(unique(m_and_pembryo1$Geneid.x))
m_and_pembryo2 <- m_and_pembryo[which(m_and_pembryo$log2FoldChange.x < 0 & m_and_pembryo$log2FoldChange.y < 0),] ## down in both
length(unique(m_and_pembryo1$Geneid.x))
m_and_pembryo3 <- m_and_pembryo[which(m_and_pembryo$log2FoldChange.x > 0 & m_and_pembryo$log2FoldChange.y < 0),] ## up in m down in p
length(unique(m_and_pembryo1$Geneid.x))
m_and_pembryo4 <- m_and_pembryo[which(m_and_pembryo$log2FoldChange.x < 0 & m_and_pembryo$log2FoldChange.y > 0),] ## up in p down in m
length(unique(m_and_pembryo1$Geneid.x))

same_direction <- rbind(m_and_pembryo1,m_and_pembryo2)
same_direction <- na.omit(same_direction)
length(same_direction)
#same_direction <- same_direction[order(same_direction$log2FoldChange.x, decreasing = TRUE),]
#same_direction["gene_order_tag"] <- paste(same_direction$log2FoldChange.x,same_direction$Geneid.x)
plot(same_direction$log2FoldChange.x, col = col2alpha("green"), ylim = c(-7, 7), cex.axis=1.5, xlab = "", xaxt = "n", ylab = "log 2 fold change", pch = 16)
par(new=TRUE) 
plot(same_direction$log2FoldChange.y, col = col2alpha("blue"), ylim = c(-7, 7), cex.axis=1.5, xlab = "", xaxt = "n", ylab = "log 2 fold change", pch = 16)

barplot(m_and_pembryo1$log2FoldChange.x, col = col2alpha("blue",0.2), ylim = c(0,7), cex.axis=4)
par(new=TRUE) 
barplot(m_and_pembryo1$log2FoldChange.y, col = col2alpha("green",0.2), ylim = c(0,7), cex.axis=4)

barplot(m_and_pembryo2$log2FoldChange.x, col = col2alpha("blue",0.2), ylim = c(-7,0), cex.axis=2)
par(new=TRUE) 
barplot(m_and_pembryo2$log2FoldChange.y, col = col2alpha("green",0.2), ylim = c(-7,0), cex.axis=2)



opp_direction <- rbind(m_and_pembryo3,m_and_pembryo4)
opp_direction <- na.omit(opp_direction)
length(opp_direction)
#opp_direction <- opp_direction[order(opp_direction$log2FoldChange.x, decreasing = TRUE),]
#opp_direction["gene_order_tag"] <- paste(opp_direction$log2FoldChange.x,opp_direction$Geneid.x)
barplot(opp_direction$log2FoldChange.x, col = col2alpha("blue",0.2), ylim = c(-7, 7))
par(new=TRUE) 
barplot(opp_direction$log2FoldChange.y, col = col2alpha("green",0.2), ylim = c(-7, 7), names.arg = opp_direction$symbol.y)

##opposites
m_and_pjaw1 <- m_and_pjaw[which(m_and_pjaw$log2FoldChange.x > 0 & m_and_pjaw$log2FoldChange.y < 0 & m_and_pjaw$protein_accession != "NA"),] ## up in m down in p
length(unique(m_and_pjaw1$protein_accession))
m_and_pjaw1["direction"] <- "up in m down in p"
m_and_pjaw2 <- m_and_pjaw[which(m_and_pjaw$log2FoldChange.x < 0 & m_and_pjaw$log2FoldChange.y > 0 & m_and_pjaw$protein_accession != "NA"),] ## up in p down in m
length(unique(m_and_pjaw1$protein_accession))
m_and_pjaw2["direction"] <- "up in p down in m"
m_and_pembryo1 <- m_and_pembryo[which(m_and_pembryo$log2FoldChange.x > 0 & m_and_pembryo$log2FoldChange.y < 0 & m_and_pembryo$protein_accession != "NA"),] ## up in m down in p
length(unique(m_and_pembryo1$protein_accession))
m_and_pembryo1["direction"] <- "up in m down in p"
m_and_pembryo2 <- m_and_pembryo[which(m_and_pembryo$log2FoldChange.x < 0 & m_and_pembryo$log2FoldChange.y > 0 & m_and_pembryo$protein_accession != "NA"),] ## up in p down in m
length(unique(m_and_pembryo1$protein_accession))
m_and_pembryo2["direction"] <- "up in p down in m"
#write.csv(opposites, "C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/mapping_correction/opposites.csv") 
head(opposites)
opposites <- rbind(m_and_pjaw1, m_and_pjaw2, m_and_pembryo1, m_and_pembryo2)


m_zeb <- read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/blast/am_jawTopHits.csv", header = TRUE)
p_zeb <- read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/blast/ap_jawTopHits.csv", header = TRUE)
m_zeb$related_accession.y <- m_zeb$V1
p_zeb$related_accession.y <- p_zeb$V1
m_zeb$comp <- "am_jaw"
p_zeb$comp <- "ap_jaw"
keeps <- c("GenBank.Protein.Accession", "related_accession.y", "V12", "Gene.Symbol", "comp")
p_zeb <- p_zeb[keeps]
m_zeb <- m_zeb[keeps]
head(p_zeb)
head(m_zeb)
zeb_key <- rbind(m_zeb, p_zeb)
nrow(zeb_key)
nrow(m_and_pjaw)
zeb_merge <- merge(m_and_pjaw, zeb_key, by = c("related_accession.y"))
nrow(zeb_merge)
head(zeb_merge)
#write.csv(zeb_merge, "C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/mapping_correction/shared_jaw_genes_zeb.csv") 


#write.csv(m_and_pjaw, "C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/post_qc/mapping_correction/m_and_p_DE_jaws_sig_test.csv") 
m_and_pembryo <- merge(am_embryo, ap_embryo, by = c("Geneid"))
keeps <- c("Geneid", "Chr.x", "log2FoldChange.x", "padj.x", "comp.x","related_accession.y", "log2FoldChange.y", "padj.y","comp.y", "symbol.y", "name.y", "related_accession.y")
m_and_pembryo <- m_and_pembryo[keeps]
length(unique(m_and_pembryo$related_accession.y))
#write.csv(m_and_pembryo, "C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/mapping_correction/m_and_p_DE_embryos_sig.csv") 

head(m_and_pembryo)

m_zeb <- read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/blast/am_embryoTopHits.csv", header = TRUE)
p_zeb <- read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/blast/ap_embryoTopHits.csv", header = TRUE)
m_zeb$related_accession.y <- m_zeb$V1
p_zeb$related_accession.y <- p_zeb$V1
m_zeb$comp <- "am_embryo"
p_zeb$comp <- "ap_embryo"
keeps <- c("GenBank.Protein.Accession", "related_accession.y", "V12", "Gene.Symbol", "comp")
p_zeb <- p_zeb[keeps]
m_zeb <- m_zeb[keeps]
head(p_zeb)
head(m_zeb)
zeb_key <- rbind(m_zeb, p_zeb)
nrow(zeb_key)
zeb_merge <- merge(m_and_pembryo, zeb_key, by = c("related_accession.y"))
nrow(zeb_merge)
head(zeb_merge)
#write.csv(m_and_pembryo, "C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/post_qc/mapping_correction/m_and_p_DE_embryos_sig.csv") 
#write.csv(zeb_merge, "C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/mapping_correction/shared_embryo_genes_zeb.csv") 

m_and_pjaw1 <- m_and_pjaw[which(m_and_pjaw$log2FoldChange.x > 0 & m_and_pjaw$log2FoldChange.y < 0),]
length(unique(m_and_pjaw1$Geneid))
m_and_pembryo1 <- m_and_pembryo[which(m_and_pembryo$log2FoldChange.x > 0 & m_and_pembryo$log2FoldChange.y > 0),]

pjaw <- c(amxp_jaw$Geneid, mp_jaw$Geneid, ap_jaw$Geneid)
length(unique(pjaw))
mjaw <- c(apxm_jaw$Geneid, mp_jaw$Geneid, am_jaw$Geneid)
length(unique(mjaw))
pembryo <- c(amxp_embryo$Geneid, mp_embryo$Geneid, ap_embryo$Geneid)
length(unique(pembryo))
membryo <- c(apxm_embryo$Geneid, mp_embryo$Geneid, am_embryo$Geneid)
length(unique(membryo))


unique(am_jaw$feature)

shared <- merge(am_jaw, ap_jaw, by = ("Geneid"))
keeps <- c("Geneid", "Chr.x", "log2FoldChange.x", "comp.x", "log2FoldChange.y", "comp.y", "symbol.y", "name.y")
shared <- shared[keeps]
shared_opposites_p_higher <- shared[which(shared$log2FoldChange.x < 0 & shared$log2FoldChange.y > 0),]
shared_opposites_m_higher <- shared[which(shared$log2FoldChange.x > 0 & shared$log2FoldChange.y < 0),]
head(shared_opposites)
nrow(shared_opposites)
#write.csv(shared, "C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/post_qc/mapping_correction/embryo_overlap_genes_sig.csv") 


species_overlap <- merge(am_jaw, ap_jaw, by = c("Geneid"))
keeps <- c("Geneid", "Chr.x", "log2FoldChange.x", "comp.x", "log2FoldChange.y", "comp.y", "symbol.y", "name.y")
species_overlap <- species_overlap[keeps]
#write.csv(species_overlap, "C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/post_qc/mapping_correction/species_jaw_overlap_genes_sig.csv") 

#head(m_stage_overlap)
#nrow(m_stage_overlap)
#p_stage_overlap <- merge(ap_jaw, mp_jaw, by = c("Geneid"))
#p_stage_overlap <- merge(p_stage_overlap, ap_embryo, by = c("Geneid"))
#p_stage_overlap <- merge(p_stage_overlap, amxp_embryo, by = c("Geneid"))
#p_stage_overlap <- merge(p_stage_overlap, amxp_jaw, by = c("Geneid"))
#p_stage_overlap <- merge(p_stage_overlap, mp_jaw, by = c("Geneid"))
#p_stage_overlap <- subset(p_stage_overlap, !duplicated(p_stage_overlap[,1])) 
#head(p_stage_overlap)
#nrow(p_stage_overlap)

#rna_seq_genes <- read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/mapping_correction/DESeq_genes_all.csv", header = TRUE, stringsAsFactors = FALSE)
rna_seq_genes <- read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/mapping_correction/DESeq_genes_sig.csv", header = TRUE, stringsAsFactors = FALSE)
indels <- read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/43_pupfish/vcfs/indels/fixed_indels_genes.csv", header = TRUE, stringsAsFactors = FALSE, fileEncoding="UTF-8-BOM")
snps <-read.table("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/43_pupfish/fst/43_nomaysim/fixed_snps_genes_nomaysim.txt", header = TRUE, stringsAsFactors = FALSE)
mbe_cans <- read.table("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/mbe_candidates.txt", header = TRUE, stringsAsFactors = FALSE)

embryos <- read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/pleiotropy/embryos.csv", header = TRUE, stringsAsFactors = FALSE, fileEncoding="UTF-8-BOM")
jaws <- read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/pleiotropy/jaws.csv", header = TRUE, stringsAsFactors = FALSE, fileEncoding="UTF-8-BOM")

head(jaws)
merge(mbe_cans, embryos, by = c("Geneid"))
merge(mbe_cans, jaws, by = c("Geneid"))

head(rna_seq_genes)
nrow(rna_seq_genes)
length(unique(rna_seq_genes$Geneid))
head(indels)
nrow(indels)
head(snps)
nrow(snps)
head(mbe_cans)

rna_seq_genes <- rna_seq_genes[which(rna_seq_genes$comp == "ap_jaw" | rna_seq_genes$comp == "am_jaw"),]
rna_seq_genes <- rna_seq_genes[which(rna_seq_genes$comp == "ap_embryo" | rna_seq_genes$comp == "am_embryo"),]
mbe_overlap <- merge(rna_seq_genes, mbe_cans, by = c("Geneid"))
length(unique(mbe_overlap$symbol))
#write.csv(mbe_overlap, "C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/post_qc/mapping_correction/mbe_overlap.csv") 


############## snps ####
snps <- read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/43_pupfish/fst/43_nomaysim/fixed_snp_genes_fixed_comparison_mp.csv", header = TRUE, stringsAsFactors = FALSE)
rna_seq_genes <- read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/mapping_correction/DESeq_genes_sig.csv", header = TRUE, stringsAsFactors = FALSE)
#write.csv(rna_seq_genes,"C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/mapping_correction/DESeq_genes_sig.csv")
rna_seq_genes <- merge(rna_seq_genes,scaffolds_key, by = c("Chr"))
rna_seq_genes <- rna_seq_genes[which(rna_seq_genes$comp == "ap_jaw"),]
#rna_seq_genes$comp == "ap_embryo"),]
#snps <- snps[which(snps$comp =="vxm"),]
head(rna_seq_genes,100) 
unique(rna_seq_genes$seq_type)
head(snps)
tail(snps)
head(scaffolds_key)
snps <- merge(snps, scaffolds_key, by = c("CHROM"))
test <- merge(rna_seq_genes, snps, by = c("CHROM"))
tail(test)


scaffs <- snps$Chr

scaffolds   <- c()
start_s     <- c()
stop_s      <- c()
gene_names <- c()
snp_pss     <- c()
features   <- c()

#test <- snps[which(snps$Chr =="NW_015150453.1"),]
#scaff <- "NW_015150453.1"

for (scaff in scaffs) 
{
  
  snp_table   <- snps[which(snps$Chr == scaff),]
  genes_table <- rna_seq_genes[which(rna_seq_genes$Chr == scaff),]
  gene_strt   <- genes_table$Start
  gene_stp    <- genes_table$End
  gene_name   <- genes_table$Geneid
  
  run <- c(1:length(gene_strt))
  if ((nrow(genes_table)) > 0 )
  {
    for (ps in snp_table$POS)
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
          
        }
        else if (ps >= (gene_strt[i] - 10000) & ps <= gene_stp[i] + 10000)
        {
          
          scaffolds <- c(scaffolds, scaff)
          start_s <- c(start_s, gene_strt[i])
          stop_s <- c(stop_s, gene_stp[i])
          gene_names <- c(gene_names, gene_name[i])
          snp_pss <- c(snp_pss,  ps)#paste(ps, "10kb", sep=':'))
          
        }
      }}}}

test5 <- data.frame(Chr = scaffolds, feature_start = start_s, feature_stop = stop_s, snp_ps = snp_pss, 
                    Geneid = gene_names, 
                    stringsAsFactors = FALSE)
head(test5)
nrow(test5)
#write.csv(test7, "C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/fst 0.8/fixed_m_snp_genes_am_jaw.csv") 
test6 <- merge(test5, rna_seq_genes, by = c("Geneid"))
head(test6)
test6$tag <- paste(test6$Chr.x, test6$snp_ps, test6$Geneid, sep=":")
head(test6)
length(test6)
length(unique(test6$tag))
length(unique(test6$Geneid))
test7 <- subset(test6, !duplicated(test6[,40]))
nrow(test7)

#all_de_genes <- read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/post_qc/mapping_correction/DESeq_genes_all.csv")
head(all_de_genes)
fixed_snp_genes <- merge(test5, all_de_genes, by = c("Chr", "Geneid"))
head(fixed_snp_genes)
nrow(fixed_snp_genes)
fixed_snp_genes$tag <- paste(fixed_snp_genes$Chr, fixed_snp_genes$Geneid, fixed_snp_genes$comp, sep=":")
length(fixed_snp_genes)
fixed_snp_genes <- subset(fixed_snp_genes, !duplicated(fixed_snp_genes[,39]))
fixed_snp_genes <- fixed_snp_genes[which(fixed_snp_genes$padj <= 0.05),]
#write.csv(fixed_snp_genes, "C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/post_qc/mapping_correction/fixed_snp_genes_overlap.sig.csv") 

fixed_p_de_ap_jaw    <-  read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/mapping_correction/fixed_p_snp_genes_ap_jaw.csv")
fixed_p_de_ap_embryo <-  read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/mapping_correction/fixed_p_snp_genes_ap_embryo.csv")
fixed_m_de_am_embryo <-  read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/mapping_correction/fixed_m_snp_genes_am_embryo.csv")
fixed_p_de_ap_jaw    <-  subset(fixed_p_de_ap_jaw     , !duplicated(fixed_p_de_ap_jaw[,41])) 
fixed_p_de_ap_embryo <-  subset(fixed_p_de_ap_embryo  , !duplicated(fixed_p_de_ap_embryo[,41])) 
fixed_m_de_am_embryo <-  subset(fixed_m_de_am_embryo  , !duplicated(fixed_m_de_am_embryo[,41])) 

m_and_p_DE_embryos <- read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/mapping_correction/m_and_p_DE_embryos_sig.csv")
m_and_p_DE_jaws <- read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/mapping_correction/m_and_p_DE_jaws_sig.csv")
head(m_and_p_DE_embryos)
nrow(m_and_p_DE_embryos)
length(unique(m_and_p_DE_embryos$Geneid))
shared <- merge(m_and_p_DE_embryos, m_and_p_DE_jaws, by = c("Geneid"))
length(unique(shared$Geneid))

m_snps <- merge(m_and_p_DE_embryos, fixed_m_de_am_embryo, by = c("Geneid"))
p_snps <- merge(m_and_p_DE_jaws, fixed_p_de_ap_jaw, by = c("Geneid"))
length(unique(p_snps$Geneid))
LOC107101197
LOC107103763


############## indels #####
scaffs <- indels$scaffold


scaffolds <- c()
start_s <- c()
stop_s <- c()
gene_names <- c()
indel_pss <- c()
features <- c()

for (scaff in scaffs) 
{
  
  indel_table <- indels[which(indels$scaffold == scaff),]
  genes_table <- rna_seq_genes[which(rna_seq_genes$Chr == scaff),]
  gene_strt <- genes_table$Start
  gene_stp <- genes_table$End
  gene_name <- genes_table$Geneid
  
  run <- c(1:length(gene_strt))
  if ((nrow(genes_table)) > 0 )
  {
    for (ps in indel_table$ind_ps)
    {
      for (i in run)
      {
        
        
        if (ps >= gene_strt[i] & ps <= gene_stp[i])
        {
          
          scaffolds <- c(scaffolds, scaff)
          start_s <- c(start_s, gene_strt[i])
          stop_s <- c(stop_s, gene_stp[i])
          gene_names <- c(gene_names, gene_name[i])
          indel_pss <- c(indel_pss, ps)
          
        }
        else if (ps >= (gene_strt[i] - 10000) & ps <= gene_stp[i] + 10000)
        {
          
          scaffolds <- c(scaffolds, scaff)
          start_s <- c(start_s, gene_strt[i])
          stop_s <- c(stop_s, gene_stp[i])
          gene_names <- c(gene_names, gene_name[i])
          indel_pss <- c(indel_pss,  ps)#paste(ps, "10kb", sep=':'))
          
        }
      }}}}

test6 <- data.frame(Chr = scaffolds, feature_start = start_s, feature_stop = stop_s, ind_ps = indel_pss, 
                    Geneid = gene_names, 
                    stringsAsFactors = FALSE)
head(test6)
nrow(test6)
#write.csv(test6, "C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/post_qc/mapping_correction/fixed_indels_overlap.csv") 

all_de_genes <- read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/post_qc/mapping_correction/DESeq_genes_all.csv")
head(all_de_genes)
fixed_indels_genes <- merge(test6, all_de_genes, by = c("Chr", "Geneid"))
head(fixed_indels_genes)
nrow(fixed_indels_genes)
fixed_indels_genes$tag <- paste(fixed_indels_genes$Chr, fixed_indels_genes$Geneid, fixed_indels_genes$comp, sep=":")
length(fixed_indels_genes)
fixed_indels_genes <- subset(fixed_indels_genes, !duplicated(fixed_indels_genes[,39]))
fixed_indels_genes <- fixed_indels_genes[which(fixed_indels_genes$padj <= 0.05),]
#write.csv(fixed_indels_genes, "C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/post_qc/mapping_correction/fixed_indels_genes_overlap.sig.csv") 


data1 <- subset(test6, !duplicated(test6[,5]))
head(data1)
nrow(data1)
data1 <- merge(data1, indels, by =c("scaffold", "ind_ps"))



test6$Geneid <- test6$gene_name
tester <- merge(test6, cts, by = c("Geneid"))
head(tester)
nrow(tester)

data <- merge(tester, rna_seq_genes, by = c("Geneid"))

data$tag <- paste(data$comp, data$Geneid, data$Chr.x,sep=":")
head(data)
length(data)
data1 <- subset(data, !duplicated(data[,75])) 
head(data1)
nrow(data1)


#write.csv(data1, "C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/post_qc/counts_fixed_indels_genes.csv")



############## any fixed snp ####
all_fst <- read.table("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/43_pupfish/fst/43_nomaysim/fst_43_no_maysim", header = TRUE, stringsAsFactor = FALSE)
head(all_fst)
head(de_genes)

snps <- all_fst[which(all_fst$vmixr_fst == 1| all_fst$vpxm_fst == 1),]
head(snps)

scaffs <- snps$CHROM

scaffolds <- c()
start_s <- c()
stop_s <- c()
gene_names <- c()
snp_pss <- c()
features <- c()

for (scaff in scaffs) 
{
  
  snp_table <- snps[which(snps$CHROM == scaff),]
  genes_table <- de_genes[which(de_genes$CHROM == scaff),]
  gene_strt <- genes_table$Start
  gene_stp <- genes_table$End
  gene_name <- genes_table$Geneid
  run <- c(1:length(gene_strt))
  
  if ((nrow(snp_table)) > 0 )
  {
    
    for (ps in snp_table$POS)
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
          
        }
        else if (ps >= (gene_strt[i] - 10000) & ps <= gene_stp[i] + 10000)
        {
          
          scaffolds <- c(scaffolds, scaff)
          start_s <- c(start_s, gene_strt[i])
          stop_s <- c(stop_s, gene_stp[i])
          gene_names <- c(gene_names, gene_name[i])
          snp_pss <- c(snp_pss,  paste(ps, "10kb", sep=':'))
          
        }
      }}}}

test7 <- data.frame(scaffold = scaffolds, feature_start = start_s, feature_stop = stop_s, snp_ps = snp_pss, 
                    gene_name = gene_names, 
                    stringsAsFactors = FALSE)
head(test7)
nrow(test7)







############## counts at fixed genic snps? ####
snps <-read.table("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/43_pupfish/fst/43_nomaysim/fixed_snps_genes_nomaysim.txt", header = TRUE, stringsAsFactors = FALSE)
cts <- read.table("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/counts_my_gtf_post_trim", header = TRUE, stringsAsFactors = FALSE)
snps$Start <- snps$feature_start
snps$End <- snps$feature_stop
snps$Chr <- snps$scaffold

cts$Chr <-  gsub(";.*", "", cts$Chr)
cts$Start <-  gsub(";.*", "", cts$Start)
cts$End <-  gsub(".*;", "", cts$End)
cts$Strand <-  gsub(";.*", "", cts$Strand)

head(cts)
head(snps)
nrow(cts)
nrow(snps)

scaffs <- snps$Chr

scaffolds <- c()
start_s <- c()
stop_s <- c()
gene_names <- c()
snp_pss <- c()
features <- c()

for (scaff in scaffs) 
{
  
  snp_table <- snps[which(snps$Chr == scaff),]
  genes_table <- cts[which(cts$Chr == scaff),]
  gene_strt <- genes_table$Start
  gene_stp <- genes_table$End
  gene_name <- genes_table$Geneid
  run <- c(1:length(gene_strt))
  
  if ((nrow(genes_table)) > 0 )
  {
    
    for (ps in snp_table$snp_ps)
    {
      for (i in run)
      {
        
        
        if (ps >= (as.numeric(gene_strt[i]) - 10000) & ps <= (as.numeric(gene_stp[i]) + 10000))
        {
          
          scaffolds <- c(scaffolds, scaff)
          start_s <- c(start_s, gene_strt[i])
          stop_s <- c(stop_s, gene_stp[i])
          gene_names <- c(gene_names, gene_name[i])
          snp_pss <- c(snp_pss, ps)
          
        }
        #else if (ps >= (gene_strt[i] - 10000) & ps <= gene_stp[i] + 10000)
        #{
        #  
        #  scaffolds <- c(scaffolds, scaff)
        #  start_s <- c(start_s, gene_strt[i])
        #  stop_s <- c(stop_s, gene_stp[i])
        #  gene_names <- c(gene_names, gene_name[i])
        #  snp_pss <- c(snp_pss,  paste(ps, "10kb", sep=':'))
        #  
        #}
      }}}}

test7 <- data.frame(scaffold = scaffolds, feature_start = start_s, feature_stop = stop_s, snp_ps = snp_pss, 
                    Geneid = gene_names, 
                    stringsAsFactors = FALSE)
head(test7)
nrow(test7)
tester <- merge(test7, cts, by = c("Geneid"))
head(tester)
nrow(tester)

data <- merge(tester, rna_seq_genes, by = c("Geneid"))

data$tag <- paste(data$comp, data$Geneid, data$Chr.x,sep=":")
head(data)
length(data)
data1 <- subset(data, !duplicated(data[,74])) 
head(data1)

#write.csv(data1, "C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/post_qc/counts_fixed_snps_genes.csv")




#### fold change conversion ####

install.packages("gtools")
library(gtools)

foldchange2logratio(15, base=2)

#### permutation ####

total_genes <- 23753
#total_genes <- 36511

#genes
#de_am_jaw <- 120
#de_ap_jaw <- 1903
#
#de_am_embryo <- 1014
#de_ap_embryo <- 5982

#isoforms
#de_am_jaw <- 123
#de_ap_jaw <- 1119
#
#de_am_embryo <- 2222 
#de_ap_embryo <- 6776

##cornell
#
actual_jaw <- 394
de_ap_jaw <- 1824
de_am_jaw <- 1008

dat <- read.table("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/cornell_data/12864_2017_3810_MOESM31_ESM", header = TRUE, stringsAsFactors = FALSE, quote = "",sep = "\t")
all_DE <- read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/blast/all_DE_master.csv", header = TRUE, stringsAsFactors = FALSE)
head(dat)
nrow(dat)
length(unique(dat$GeneGenbank))
#To identify genes in our dataset that might be contributing to skull morphological variation, 
#we found the intersection set of genes at each stage that were differentially expressed 
#(DE; FDR ??? 0.1 and log2 fold change ??? 0.2) 

#dat <- cbind(dat, colsplit(dat$GeneGenbank, ":", c("junk", "GeneID")))
#dat <- dat[which(dat$SvIO_8dpf_PValue <=0.1 | dat$DvIO_8dpf_PValue <=0.1),]
#dat_p <- dat[which(dat$SvIO_8dpf_PValue <=0.1),]
#dat_m <- dat[which(dat$DvIO_8dpf_PValue <=0.1),]
#keeps <- c(GeneGenbank, GeneName, SvIO_8dpf_PValue)
#dat$log2fcp <- abs(dat$SvIO_8dpf_logFC)
#dat$log2fcm <- abs(dat$DvIO_8dpf_logFC)
#dat_p <- dat[which(dat$SvIO_8dpf_PValue <=0.1 & dat$log2fcp  >= 0.2),]
#dat_m <- dat[which(dat$DvIO_8dpf_PValue <=0.1 & dat$log2fcm  >= 0.2),]
#dat <- cbind(dat, colsplit(dat$GeneGenbank, ":", c("junk", "GeneID")))
dat_p <- dat[which(dat$SvIO_15dpf_PValue <=0.1),]
dat_m <- dat[which(dat$DvIO_15dpf_PValue <=0.1),]

dat_p <- dat_p[which(dat_p$SvIO_15dpf_PValue >=0.2 | dat_p$SvIO_15dpf_logFC <= -0.2),]
dat_m <- dat_m[which(dat_m$DvIO_15dpf_PValue >=0.2 | dat_m$DvIO_15dpf_logFC <= -0.2),]
merge_mp <- merge(dat_m, dat_p, by = c("GeneGenbank"))
actual_jaw <- nrow(merge_mp)
de_ap_jaw  <- nrow(dat_p)
de_am_jaw  <- nrow(dat_m)
actual_jaw 
de_ap_jaw  
de_am_jaw  

dat_p <- dat[which(dat$SvIO_8dpf_PValue <=0.1),]
dat_m <- dat[which(dat$DvIO_8dpf_PValue <=0.1),]

dat_p <- dat_p[which(dat_p$SvIO_8dpf_PValue >=0.2 | dat_p$SvIO_8dpf_logFC <= -0.2),]
dat_m <- dat_m[which(dat_m$DvIO_8dpf_PValue >=0.2 | dat_m$DvIO_8dpf_logFC <= -0.2),]
merge_mp <- merge(dat_m, dat_p, by = c("GeneGenbank"))
actual_jaw<- nrow(merge_mp)
de_ap_jaw <- nrow(dat_p)
de_am_jaw <- nrow(dat_m)


dat_p <- dat[which(dat$SvIO_96hpf_PValue <=0.1 & dat$SvIO_96hpf_logFC >= 0.2),]
dat_m <- dat[which(dat$DvIO_96hpf_PValue <=0.1 & dat$DvIO_96hpf_logFC >= 0.2),]

dat_p <- dat[which(dat$SvIO_96hpf_PValue <=0.1 & dat$SvIO_96hpf_logFC <= -0.2),]
dat_m <- dat[which(dat$DvIO_96hpf_PValue <=0.1 & dat$DvIO_96hpf_logFC <= -0.2),]
merge_mp <- merge(dat_m, dat_p, by = c("GeneGenbank"))
actual_jaw<- nrow(merge_mp)
de_ap_jaw <- nrow(dat_p)
de_am_jaw <- nrow(dat_m)

overlaps <- c()
run <- (1:100000)

for (i in run)
{
  de_ap_sample <- sample(total_genes, de_am_jaw, replace = FALSE, prob = NULL)
  de_am_sample <- sample(total_genes, de_ap_jaw, replace = FALSE, prob = NULL)
  
  overlap <- length(intersect(de_ap_sample,de_am_sample))
  overlaps <- c(overlaps, overlap)
}

head(overlaps)
length(overlaps)

hist(overlaps)
quantile(overlaps, .999)
q <- quantile(overlaps, .999)
abline(v=q, col = "red", lty = 3, lwd = 1.2)
abline(v=actual_jaw , col = "green", lty = 3, lwd = 1.2)

quantile(overlaps, 1)
actual_jaw


overlaps_e <- c()
run <- (1:100000)
for (i in run)
{
  de_ap_sample <- sample(total_genes, de_am_embryo, replace = FALSE, prob = NULL)
  de_am_sample <- sample(total_genes, de_ap_embryo, replace = FALSE, prob = NULL)
  
  overlap <- length(intersect(de_ap_sample,de_am_sample))
  overlaps_e <- c(overlaps_e, overlap)
}

head(overlaps_e)
length(overlaps_e)

hist(overlaps_e)
quantile(overlaps_e, 1)
q <- quantile(overlaps, .999)
actual_embryo <- 833
abline(v=q, col = "red", lty = 3, lwd = 1.2)
abline(v=actual_embryo , col = "green", lty = 3, lwd = 1.2)
#write.csv(overlaps_e, "C:/Users/Joseph McGirr Lab/Documents/python_course/overlap_freqs_embryo.csv") 
#write.csv( overlaps,  "C:/Users/Joseph McGirr Lab/Documents/python_course/overlap_freqs_jaw.csv") 
#write.csv(overlaps_e, "C:/Users/Joseph McGirr Lab/Documents/python_course/overlap_freqs_gene_embryo.csv") 
#write.csv( overlaps,  "C:/Users/Joseph McGirr Lab/Documents/python_course/overlap_freqs_gene_jaw.csv") 


jaw_i <- read.csv("C:/Users/Joseph McGirr Lab/Documents/python_course/overlap_freqs_embryo.csv")
emb_i <- read.csv("C:/Users/Joseph McGirr Lab/Documents/python_course/overlap_freqs_embryo.csv")
emb_g <- as.data.frame(overlaps_e)
jaw_g <- as.data.frame(overlaps)
head(jaw_g)
head(jaw_i)
jaw_i$jaw_i <-  jaw_i$x           
emb_i$emb_i <-  emb_i$x             
emb_g$emb_g <-  emb_g$overlaps_e  
jaw_g$jaw_g <-  jaw_g$overlaps 
all <- data.frame()
all$jaw_i <- as.numeric(jaw_i$jaw_i)
all$emb_i <- as.numeric(emb_i$emb_i)
all$emb_g <- as.numeric(emb_g$emb_g)
all$jaw_g <- as.numeric(jaw_g$jaw_g)


##### fixed snps in coding regions? #####
##### by CDS #####
gene_locations <- read.table("C:/Users/Joseph McGirr Lab/Desktop/Cyprinodon/feature_table.txt", header = TRUE, stringsAsFactors = FALSE, quote = "", fill = TRUE, sep = "\t")
snps <- read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/43_pupfish/fst/43_nomaysim/fixed_snp_genes_fixed_comparison_mp.csv", header = TRUE, stringsAsFactors = FALSE)
head(gene_locations)
unique(gene_locations$feature)
head(snps)
gene_locations <- gene_locations[which(gene_locations$feature == "CDS"),]

snps$genomic_accession <- snps$scaffold

scaffs <- snps$genomic_accession

scaffolds <- c()
start_s <- c()
stop_s <- c()
gene_names <- c()
snp_pss <- c()
features <- c()

#test <- snps[which(snps$CHROM =="KL653420.1"),]

for (scaff in scaffs) 
{
  
  snp_table <- snps[which(snps$genomic_accession == scaff),]
  genes_table <- gene_locations[which(gene_locations$genomic_accession == scaff),]
  gene_strt <- genes_table$start
  gene_stp <- genes_table$end
  gene_name <- genes_table$symbol
  
  run <- c(1:length(gene_strt))
  if ((nrow(genes_table)) > 0 )
  {
    for (ps in snp_table$snp_ps)
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
          
        }
        
      }}}}

test5 <- data.frame(Chr = scaffolds, feature_start = start_s, feature_stop = stop_s, snp_ps = snp_pss, 
                    Geneid = gene_names, 
                    stringsAsFactors = FALSE)
head(test5)
nrow(test5)
#write.csv(test7 , "C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/mapping_correction/cds_snps.csv") 
test5$tag <- paste(test5$Chr, test5$snp_ps, sep=":")
test6   <- subset(test5    , !duplicated(test5[,6]))
head(test6)
nrow(test6)
head(snps)
snps$tag <- paste(snps$scaffold, snps$snp_ps, sep=":")
test7 <- merge(test6, snps, by = c("tag"))
head(test7)
nrow(test7)
length(unique(test6$tag))
length(unique(test6$Geneid))
tester <- merge(test6, m_and_pembryo, by =c("Geneid"))
tester2 <- tester[unique(tester$tag),]
length(unique(tester$Geneid))
head(tester)


##### by exons #####
gene_locations <- read.table("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/mapping_correction/my.gtf.geneid", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
snps <- read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/43_pupfish/fst/43_nomaysim/fixed_snp_genes_fixed_comparison_mp.csv", header = TRUE, stringsAsFactors = FALSE)
head(gene_locations)
gene_locations <- gene_locations[which(gene_locations$V3 == "exon"),]

head(snps)

scaffs <- snps$scaffold

scaffolds <- c()
start_s <- c()
stop_s <- c()
gene_names <- c()
snp_pss <- c()
features <- c()

#test <- snps[which(snps$CHROM =="KL653420.1"),]

for (scaff in scaffs) 
{
  
  snp_table <- snps[which(snps$scaffold == scaff),]
  genes_table <- gene_locations[which(gene_locations$V1 == scaff),]
  gene_strt <- genes_table$V4
  gene_stp <- genes_table$V5
  gene_name <- genes_table$V9
  
  run <- c(1:length(gene_strt))
  if ((nrow(genes_table)) > 0 )
  {
    for (ps in snp_table$snp_ps)
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
          
        }
        
      }}}}

test5 <- data.frame(Chr = scaffolds, feature_start = start_s, feature_stop = stop_s, snp_ps = snp_pss, 
                    Geneid = gene_names, 
                    stringsAsFactors = FALSE)
head(test5)
nrow(test5)
#write.csv(test8 , "C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/mapping_correction/exon_snps.csv") 
test5$tag <- paste(test5$Chr, test5$snp_ps, sep=":")
test6   <- subset(test5    , !duplicated(test5[,6]))
head(test6)
nrow(test6)
head(snps)
snps$tag <- paste(snps$scaffold, snps$snp_ps, sep=":")
test7 <- merge(test6, snps, by = c("tag"))
#fixed_genes_DE <- read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/mapping_correction/fixed_de_all.csv", header = TRUE, stringsAsFactors = FALSE)
test9 <- merge(test8, fixed_genes_DE, by = c("symbol"))
test8 <- merge(test7, all_DE, by = c("Chr", "symbol"))
head(all_DE)
head(test9)
nrow(test7)


#### fixed snps near genes DE for both ap and am ####
snps <- read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/43_pupfish/fst/43_nomaysim/fixed_snp_genes_fixed_comparison_mp.csv", header = TRUE, stringsAsFactors = FALSE)
#rna_seq_genes <- read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/mapping_correction/DESeq_genes_sig.csv", header = TRUE, stringsAsFactors = FALSE)
rna_seq_genes <- read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/blast/all_DE_master.csv", header = TRUE, stringsAsFactors = FALSE)

#rna_seq_genes <- rna_seq_genes[which(rna_seq_genes$comp == "ap_embryo"),]
#rna_seq_genes$comp == "am_embryo"),]
#snps <- snps[which(snps$comp =="vxr"),]
head(rna_seq_genes) 
head(snps)
snps$CHROM <- snps$scaffold
snps$Chr <- snps$CHROM

scaffs <- snps$CHROM

scaffolds <- c()
start_s <- c()
stop_s <- c()
gene_names <- c()
snp_pss <- c()
features <- c()

#test <- snps[which(snps$CHROM =="KL653420.1"),]

for (scaff in scaffs) 
{
  
  snp_table <- snps[which(snps$CHROM == scaff),]
  genes_table <- rna_seq_genes[which(rna_seq_genes$Chr == scaff),]
  gene_strt <- genes_table$Start
  gene_stp <- genes_table$End
  gene_name <- genes_table$Geneid
  
  run <- c(1:length(gene_strt))
  if ((nrow(genes_table)) > 0 )
  {
    for (ps in snp_table$snp_ps)
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
          
        }
        else if (ps >= (gene_strt[i] - 10000) & ps <= gene_stp[i] + 10000)
        {
          
          scaffolds <- c(scaffolds, scaff)
          start_s <- c(start_s, gene_strt[i])
          stop_s <- c(stop_s, gene_stp[i])
          gene_names <- c(gene_names, gene_name[i])
          snp_pss <- c(snp_pss,  ps)#paste(ps, "10kb", sep=':'))
          
        }
      }}}}

test5 <- data.frame(Chr = scaffolds, Start = start_s, End = stop_s, snp_ps = snp_pss, 
                    stringsAsFactors = FALSE)
head(test5)
nrow(test5)
#write.csv(test, "C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/blast/test.csv") 
test6 <- merge(test5, rna_seq_genes, by = c("Chr", "Start", "End"))
test6$tag <- paste(test6$Chr, test6$comp, test6$zeb_homolog, test6$symbol, sep=":")
head(test6)
length(unique(test6$tag))
length(unique(test6$Geneid))
length(test6)
test7   <- subset(test6    , !duplicated(test6[,28]))
head(test7)

test <- merge(test7, snps, by = c("Chr", "snp_ps"))

tester <- merge(test6, m_and_pembryo, by =c("Geneid"))
tester2 <- tester[unique(tester$tag),]
length(unique(tester$Geneid))
head(tester)

#all_de_genes <- read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/post_qc/mapping_correction/DESeq_genes_all.csv")
#head(all_de_genes)
#fixed_snp_genes <- merge(test5, all_de_genes, by = c("Chr", "Geneid"))
#head(fixed_snp_genes)
#nrow(fixed_snp_genes)
#fixed_snp_genes$tag <- paste(fixed_snp_genes$Chr, fixed_snp_genes$Geneid, fixed_snp_genes$comp, sep=":")
#length(fixed_snp_genes)
#fixed_snp_genes <- subset(fixed_snp_genes, !duplicated(fixed_snp_genes[,39]))
#fixed_snp_genes <- fixed_snp_genes[which(fixed_snp_genes$padj <= 0.05),]
#write.csv(fixed_snp_genes, "C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/post_qc/mapping_correction/fixed_snp_genes_overlap.sig.csv") 

#### blast zebrafish homologs ####
hits <- read.table("C:/Users/Joseph McGirr Lab/Desktop/ubuntushare/blast/ncbi-blast-2.6.0+/bin/fixed_snp_genes_blast_hits", header = FALSE, stringsAsFactors = FALSE, sep = "\t")
head(hits)
length(unique(hits$V1))
#hits[which(hits$V1 == "XP_015258556.1"),]

match_table <- data.frame()
head(match_table)
proteins <- unique(hits$V1) 

for (i in proteins)
{
  protein_table <- hits[which(hits$V1 == i),]
  match <- protein_table[order(protein_table$V12, decreasing = TRUE),]
  matches <- match[1,]
  match_table <- rbind(match_table, matches)
}


head(match_table)
nrow(match_table)
match_table <- cbind(match_table, colsplit(match_table$V2, "\\.", c("GenBank.Protein.Accession", "junk")))
head(match_table)
nrow(match_table)
#write.csv(match_table, "C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/blast/fixed_snp_genes_TopHits.csv") 
zeb <- read.table("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/blast/fixed_snp_genes_zebra_acession_symbol_key.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
head(zeb)
final <- merge(match_table, zeb, by = c("GenBank.Protein.Accession"))
head(final)
nrow(final)
final$protein_accession <- final$V1
final$zeb_homolog <- final$Gene.Symbol
final$zeb_protein_accession <- final$V2
keeps <- c("protein_accession", "zeb_protein_accession", "zeb_homolog")
final_zeb <- final[keeps]
head(final_zeb)
#write.csv(mp_DE_zeb , "C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/blast/mp_sig_zeb.csv") 
mp_DE <- read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/blast/mp_sig.csv", header = TRUE, stringsAsFactors = FALSE)
head(mp_DE)

mp_DE_zeb <- merge(mp_DE, final_zeb, all.x = TRUE, by = c("protein_accession"))
nrow(mp_DE)
nrow(mp_DE_zeb)
head(mp_DE_zeb)

?merge
pa <- read.table("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/sweeps/fixed_snp_genes_pa.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
fsg <- read.table("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/sweeps/fixed_snp_genes.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")

head(pa)
pa$protein_accession1 <- pa$Symbol
final_zeb <- cbind(final_zeb, colsplit(final_zeb$protein_accession, "\\.", c("protein_accession1", "junk")))
test <- merge(final_zeb, pa, by = c("protein_accession1"))
tail(fsg)
head(test)
test$symbol <- test$Gene
final <- merge(fsg, test, by = c("symbol"))
head(final)
unique(final$fixed_in)
unique(test$Gene)
#write.csv(final , "C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/sweeps/test.csv") 


#### pleiotropy ####
##### GO terms ####
#install.packages("Rcpp")
zfin <- read.table("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/pleiotropy/gene_association.jam.zfin", fill=TRUE,header = FALSE,quote = "", stringsAsFactors = FALSE, sep = "\t")
head(zfin)
zfin <- zfin[which(zfin$V9 == "P"),]
zfin <- zfin[which(zfin$V7 == "EXP" | zfin$V7 == "IDA"| zfin$V7 == "IPI" | zfin$V7 == "IMP" | zfin$V7 == "IGI" | zfin$V7 == "IEP"),]
all_DE <- read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/blast/all_DE_master.csv", header = TRUE, stringsAsFactors = FALSE)
head(all_DE)
ensembl <- read.table("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/pleiotropy/ensembl_protein_accessions.txt", header = TRUE, stringsAsFactors = FALSE)
head(ensembl)
ensembl$ones <- 1
ensembl$zeb_protein_accession <- paste(ensembl$GenBank_Protein_Accession, ensembl$ones, sep=".")
keeps <- c("zeb_protein_accession", "EnsemblProteinID")
ensembl <- ensembl[keeps]
#new_master <- merge(all_DE, ensembl, all.x = TRUE, by = c("zeb_protein_accession"))
#nrow(all_DE)
#nrow(new_master)
#write.csv(zeb_hom_fixed_genes_DE, "C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/mapping_correction/fixed_de_all_zeb_hom.csv")
head(fixed_genes_DE)
fixed_genes_DE$protein_accession <- fixed_genes_DE$related_accession
zeb_hom_fixed_genes_DE <- merge(fixed_genes_DE, all_DE, all.x = TRUE, by = c("protein_accession"))
head(zeb_hom_fixed_genes_DE)

zfin <- zfin[order(zfin$V3),]
zfin$zeb_homolog <- zfin$V3
zfin$go <- zfin$V5
keeps <- c("zeb_homolog", "go")
go <- zfin[keeps]
head(go)

go_merge <- merge(all_DE, go, by = c("zeb_homolog"))
head(go_merge)
nrow(go_merge)
length(unique(go_merge$zeb_homolog))
go_merge <- go_merge[order(go_merge$zeb_homolog),]
go_counts <-  count(go_merge, 'zeb_homolog')
head(go_counts)
nrow(go_counts)


#ap_jaw <- go_merge[which(go_merge$comp == "ap_jaw"),]
#am_jaw <- go_merge[which(go_merge$comp == "am_jaw"),]
#am_embryo <- go_merge[which(go_merge$comp == "ap_embryo"),]
#ap_embryo <- go_merge[which(go_merge$comp == "am_embryo"),]
#ap_jaw$tag <- paste(ap_jaw$zeb_homolog, ap_jaw$go, sep=":") 
#am_jaw$tag <- paste(am_jaw$zeb_homolog, am_jaw$go, sep=":") 
#am_embryo$tag <- paste(am_embryo$zeb_homolog, am_embryo$go, sep=":") 
#ap_embryo$tag <- paste(ap_embryo$zeb_homolog, ap_embryo$go, sep=":") 
#head(ap_jaw)
#length(ap_jaw)
#
#
#ap_j_go_counts <-  count(ap_jaw,    'zeb_homolog')
#am_j_go_counts <-  count(am_jaw,    'zeb_homolog')
#am_e_go_counts <-  count(am_embryo, 'zeb_homolog')
#ap_e_go_counts <-  count(ap_embryo, 'zeb_homolog')
#head(ap_j_go_counts)
#
#
#ap_j_go_freqs <- merge(ap_j_go_counts, ap_jaw   , by = c("zeb_homolog"))
#am_j_go_freqs <- merge(am_j_go_counts, am_jaw   , by = c("zeb_homolog"))
#am_e_go_freqs <- merge(am_e_go_counts, am_embryo, by = c("zeb_homolog"))
#ap_e_go_freqs <- merge(ap_e_go_counts, ap_embryo, by = c("zeb_homolog"))
#test <- ap_e_go_freqs[order(ap_e_go_freqs$freq, decreasing = TRUE),]
#head(test)
#
#ap_j_go_freqs <- subset(ap_j_go_freqs, !duplicated(ap_j_go_freqs[,1]))
#am_j_go_freqs <- subset(am_j_go_freqs, !duplicated(am_j_go_freqs[,1])) 
#am_e_go_freqs <- subset(am_e_go_freqs, !duplicated(am_e_go_freqs[,1])) 
#ap_e_go_freqs <- subset(ap_e_go_freqs, !duplicated(ap_e_go_freqs[,1])) 


ap_jaw <- go_merge[which(go_merge$comp == "ap_jaw"),]
am_jaw <- go_merge[which(go_merge$comp == "am_jaw"),]
am_embryo <- go_merge[which(go_merge$comp == "ap_embryo"),]
ap_embryo <- go_merge[which(go_merge$comp == "am_embryo"),]
ap_jaw$tag <- paste(ap_jaw$zeb_homolog, ap_jaw$go, sep=":") 
am_jaw$tag <- paste(am_jaw$zeb_homolog, am_jaw$go, sep=":") 
am_embryo$tag <- paste(am_embryo$zeb_homolog, am_embryo$go, sep=":") 
ap_embryo$tag <- paste(ap_embryo$zeb_homolog, ap_embryo$go, sep=":") 


jaws <- merge(ap_jaw, am_jaw, by = c("tag", "protein_accession"))
head(jaws)
nrow(jaws)
length(unique(jaws$tag))
length(unique(jaws$protein_accession))
head(ap_jaw)
shared <- unique(jaws$protein_accession)
length(shared)
in_p_not_m <- ap_jaw$protein_accession[!(ap_jaw$protein_accession %in% am_jaw$protein_accession)]
in_m_not_p <- am_jaw$protein_accession[!(am_jaw$protein_accession %in% ap_jaw$protein_accession)]
unshared <- c(in_p_not_m, in_m_not_p)
shared <- as.data.frame(shared)
shared$protein_accession <- shared$shared
keeps <- c("protein_accession")
shared <- shared[keeps]
unshared <- as.data.frame(unshared)
unshared$protein_accession <- unshared$unshared
unshared <- unshared[keeps]
head(unshared)

all_jaws <- rbind(ap_jaw, am_jaw)
head(all_jaws)

unshared_go <- merge(unshared, all_jaws, by = ("protein_accession"))
length(unshared_go)
unshared_go    <- subset(unshared_go    , !duplicated(unshared_go[,28]))
head(unshared_go)
length(unique(unshared_go$protein_accession))
unshared_counts <-  count(unshared_go,    'zeb_homolog')
head(unshared_counts)
unshared_freqs <- merge(unshared_counts, unshared_go, by = ("zeb_homolog"))
colnames(unshared_freqs)

head(shared)
shared_go <- merge(shared, all_jaws, by = ("protein_accession"))
length(shared_go)
shared_go    <- subset(shared_go    , !duplicated(shared_go[,28]))
shared_counts <-  count(shared_go,    'zeb_homolog')
head(shared_counts)
shared_freqs <- merge(shared_counts, shared_go, by = ("zeb_homolog"))
colnames(shared_freqs)
shared_freqs$overlap <- "shared"
unshared_freqs$overlap <- "unshared"
all_freqs <- rbind(shared_freqs, unshared_freqs)
all_freqs   <- subset(all_freqs    , !duplicated(all_freqs[,1]))
head(all_freqs)
nrow(all_freqs)

#install.packages("MASS")
#library(MASS)
#go_lm <- glm(all_freqs$freq~all_freqs$overlap, family = poisson)
go_lm <- glm.nb(all_freqs$freq~all_freqs$overlap)
summary (go_lm)
nrow(all_freqs[which(all_freqs$overlap == "shared"),])
nrow(all_freqs[which(all_freqs$overlap == "unshared"),])
boxplot(all_freqs$freq~all_freqs$overlap, main = "jaws")
abline(go_lm)

#write.csv(all_freqs, "C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/pleiotropy/jaws.csv") 
####
jaws <- merge(ap_embryo, am_embryo, by = c("tag", "protein_accession"))
head(jaws)
nrow(jaws)
length(unique(jaws$tag))
length(unique(jaws$protein_accession))
head(ap_jaw)
shared <- unique(jaws$protein_accession)
length(shared)
in_p_not_m <- ap_embryo$protein_accession[!(ap_embryo$protein_accession %in% am_embryo$protein_accession)]
in_m_not_p <- am_embryo$protein_accession[!(am_embryo$protein_accession %in% ap_embryo$protein_accession)]
unshared <- c(in_p_not_m, in_m_not_p)
length(unique(unshared))
shared <- as.data.frame(shared)
shared$protein_accession <- shared$shared
keeps <- c("protein_accession")
shared <- shared[keeps]
unshared <- as.data.frame(unshared)
unshared$protein_accession <- unshared$unshared
unshared <- unshared[keeps]
head(unshared)

all_jaws <- rbind(ap_embryo, am_embryo)
head(all_jaws)

unshared_go <- merge(unshared, all_jaws, by = ("protein_accession"))
length(unshared_go)
head(unshared_go)
unshared_go    <- subset(unshared_go    , !duplicated(unshared_go[,28]))
head(unshared_go)
length(unique(unshared_go$protein_accession))
unshared_counts <-  count(unshared_go,    'zeb_homolog')
head(unshared_counts)
unshared_freqs <- merge(unshared_counts, unshared_go, by = ("zeb_homolog"))
colnames(unshared_freqs)

head(shared)
shared_go <- merge(shared, all_jaws, by = ("protein_accession"))
length(shared_go)
shared_go    <- subset(shared_go    , !duplicated(shared_go[,28]))
shared_counts <-  count(shared_go,    'zeb_homolog')
head(shared_counts)
shared_freqs <- merge(shared_counts, shared_go, by = ("zeb_homolog"))
colnames(shared_freqs)
shared_freqs$overlap <- "shared"
unshared_freqs$overlap <- "unshared"
all_freqs <- rbind(shared_freqs, unshared_freqs)
all_freqs   <- subset(all_freqs    , !duplicated(all_freqs[,1]))
head(all_freqs)
nrow(all_freqs)


go_lm <- glm(all_freqs$freq~all_freqs$overlap, family = poisson)
go_lm <- glm(all_freqs$freq~all_freqs$overlap, family = quasipoisson)
go_lm <- glm.nb(all_freqs$freq~all_freqs$overlap)
dispersiontest(go_lm)
#go_lm <- lm(all_freqs$freq~all_freqs$overlap)
go_lm.resi <- resid(go_lm)
qqnorm(go_lm.resi)
qqline(go_lm.resi)
summary (go_lm)
abline(go_lm)
nrow(all_freqs[which(all_freqs$overlap == "shared"),])
nrow(all_freqs[which(all_freqs$overlap == "unshared"),])
boxplot(all_freqs$freq~all_freqs$overlap, main = "embryos", ylab = "GO Terms")



#write.csv(all_freqs, "C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/pleiotropy/embryos.csv") 


##### REVIGO #####
shared   <- read.table("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/blast/revigo_shared.txt",   fill=TRUE,header = TRUE,quote = "", stringsAsFactors = FALSE, sep = "\t")
unshared <- read.table("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/blast/revigo_unshared.txt", fill=TRUE,header = TRUE,quote = "", stringsAsFactors = FALSE, sep = "\t")
head(shared)
go_counts_shared <-  count(shared, 'representative')
head(go_counts,50)
go_counts_unshared <-  count(unshared, 'representative')
#write.csv(go_counts_unshared, "C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/blast/test2.csv")

shared   <- read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/blast/revigo_shared_pie.csv",   fill=TRUE,header = TRUE,quote = "", stringsAsFactors = FALSE)
unshared <- read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/blast/revigo_unshared_pie.csv", fill=TRUE,header = TRUE,quote = "", stringsAsFactors = FALSE)
head(unshared)
unique(shared$col_group_gen)
pie(shared$freq, col= c("firebrick1", "chartreuse4", "blue3")[as.numeric(shared$col_group_gen)])
pie(unshared$freq, col= c("firebrick1", "chartreuse4", "blue3")[as.numeric(unshared$col_group_gen)])

unshared_dev <- unshared[which(unshared$)]




##### GO term overlap between shared and unshared #####
#go_shared   <- read.table("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/pleiotropy/go_enrichment_shared.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
#go_unshared <- read.table("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/pleiotropy/go_enrichment_unshared.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")

go_m      <- read.table("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/pleiotropy/GO_enrich_unshared_am_embryo.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
go_p      <- read.table("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/pleiotropy/GO_enrich_unshared_ap_embryo.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
go_shared <- read.table("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/pleiotropy/GO_enrich_shared_embryo.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
go_mp     <- read.table("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/pleiotropy/GO_enrich_mp_embryo.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
nrow(go_m      )
nrow(go_p      )
nrow(go_shared )
nrow(go_mp     )
head(go_m      )


length(unique(go_m$GO.biological.process.complete))
length(unique(go_p$GO.biological.process.complete))
length(unique(go_shared$GO.biological.process.complete))
length(unique(go_mp$GO.biological.process.complete))


zfin <- read.table("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/pleiotropy/gene_association.jam.zfin", fill=TRUE,header = FALSE,quote = "", stringsAsFactors = FALSE, sep = "\t")

head(zfin)
nrow(zfin)
head(go_shared)
head(go_m)
nrow(go_p)
nrow(go_shared)
nrow(go_m)
colnames(go_shared)


same_go <- merge(go_m, go_p, by = c("GO.biological.process.complete"))
nrow(same_go)
head(same_go)

same_go <- cbind(same_go, colsplit(same_go$GO.biological.process.complete, "\\(", c("names", "go")))
head(same_go)
go <- gsub( ")", "", as.character(same_go$go))
same_go$V5 <- go
length(unique(same_go$V5))

go_types <- merge(same_go, zfin, by = c("V5"))
head(go_types)
nrow(go_types)
length(unique(go_types$V5))

same_go <- cbind(go_shared, colsplit(go_shared$GO.biological.process.complete, "\\(", c("names", "go")))
head(same_go)
go <- gsub( ")", "", as.character(same_go$go))
same_go$V5 <- go
go_ids <- unique(same_go$V5)
final <- as.data.frame(go_ids)
#write.table(final , "C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/blast/go_shared_ids.txt", sep = "\t", quote = FALSE) 



##### PPI ####

all_DE <- read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/pleiotropy/all_DE_master_with_ensembl.csv", header = TRUE, stringsAsFactors = FALSE)
head(all_DE)
ppi <- read.table("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/pleiotropy/danio_ppi_experimental", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
head(ppi)
range(ppi$experimental)
range(ppi$combined_score)
quantile(ppi$experimental,.99)
quantile(ppi$combined_score,.99)
ppi_sig <- ppi[which(ppi$combined_score >= 781 & ppi$experimental >= 427),]
ppi_sig <- ppi
ppi_counts <-  count(ppi_sig, 'protein1')
head(ppi_counts)

ensembl <- read.table("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/pleiotropy/ensembl_protein_accessions.txt", header = TRUE, stringsAsFactors = FALSE)
head(ensembl)
ensembl$ones <- 1
ensembl$zeb_protein_accession <- paste(ensembl$GenBank_Protein_Accession, ensembl$ones, sep=".")
keeps <- c("zeb_protein_accession", "EnsemblProteinID")
ensembl <- ensembl[keeps]
jaws <- read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/pleiotropy/jaws.csv", header = TRUE, stringsAsFactors = FALSE)
jaws <- merge(jaws, ensembl, all.x = TRUE, by = c("zeb_protein_accession"))
jaws$protein1 <- jaws$EnsemblProteinID
ppi_jaws <- merge(jaws, ppi_counts, by = c("protein1"))
ppi_jaws   <- subset(ppi_jaws    , !duplicated(ppi_jaws[,1]))
head(ppi_jaws)
nrow(ppi_jaws)
nrow(ppi_jaws[which(ppi_jaws$overlap == "shared"),])
nrow(ppi_jaws[which(ppi_jaws$overlap == "unshared"),])
boxplot(ppi_jaws$freq.y~ppi_jaws$overlap, main = "jaws")
ppi_lm <- glm(ppi_jaws$freq.y~ppi_jaws$overlap, family = poisson)
summary (ppi_lm)
#write.csv(ppi_jaws, "C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/pleiotropy/ppi_jaws.csv") 


#install.packages("AER")
#library(AER)
embryos <- read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/pleiotropy/embryos.csv", header = TRUE, stringsAsFactors = FALSE)
embryos <- merge(embryos, ensembl, all.x = TRUE, by = c("zeb_protein_accession"))
embryos$protein1 <- embryos$EnsemblProteinID
ppi_embryos <- merge(embryos, ppi_counts, by = c("protein1"))
ppi_embryos   <- subset(ppi_embryos    , !duplicated(ppi_embryos[,1]))
head(ppi_embryos)
nrow(ppi_embryos)
ppi_embryos <- ppi_embryos[order(ppi_embryos$freq.y, decreasing = TRUE),]
head(ppi_embryos)
nrow(ppi_embryos[which(ppi_embryos$overlap == "shared"),])
nrow(ppi_embryos[which(ppi_embryos$overlap == "unshared"),])
boxplot(ppi_embryos$freq.y~ppi_embryos$overlap, main = "embryos", ylab = "Protein Interactions")
#ppi_lm <- glm(ppi_embryos$freq.y~ppi_embryos$overlap, family = poisson)
#ppi_lm <- glm(ppi_embryos$freq.y~ppi_embryos$overlap, family = quasipoisson)
ppi_lm <- glm.nb(ppi_embryos$freq.y~ppi_embryos$overlap)
dispersiontest(ppi_lm)
ppi_lm.resi <- resid(ppi_lm)
qqnorm(ppi_lm.resi)
qqline(ppi_lm.resi)
summary (ppi_lm)
1 - pchisq(deviance(ppi_lm),df.residual(ppi_lm))

#write.csv(ppi_embryos, "C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/pleiotropy/ppi_embryos.csv") 

hist(residuals(ppi_lm,type = c("deviance")))
?glm
?log


#source("https://bioconductor.org/biocLite.R")
#biocLite("BgeeDB")
library(BgeeDB)
myTopAnatData <- loadTopAnatData(species=7955, datatype="in_situ")
head(myTopAnatData$gene2anatomy)


##### tissue expression #####
tis <- read.table("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/pleiotropy/Danio_rerio_expr_simple_development.tsv", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
#head(tis)
#tis <- tis[which(tis$Call.quality == "gold quality"),]
#unique(tis$Anatomical.entity.name)
#tis$tag <- paste(tis$Gene.name, tis$Anatomical.entity.name, sep = ":")
#length(tis)
#tis <- subset(tis    , !duplicated(tis_jaws[,10]))


#embryos <- read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/pleiotropy/embryos.csv", header = TRUE, stringsAsFactors = FALSE)
#jaws <- read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/pleiotropy/jaws.csv", header = TRUE, stringsAsFactors = FALSE)
#tis_counts <- count(tis, 'Gene.name')
#nrow(tis_counts)
#head(tis_counts)
#jaws$Gene.name <- jaws$zeb_homolog
#head(jaws)
#tis_jaws <- merge(jaws, tis_counts, all.x = TRUE, by = c("Gene.name"))
#head(tis_jaws)
#tis_jaws   <- subset(tis_jaws    , !duplicated(tis_jaws[,1]))
#nrow(tis_jaws[which(tis_jaws$overlap == "shared"),])
#nrow(tis_jaws[which(tis_jaws$overlap == "unshared"),])
#boxplot(tis_jaws$freq.y~tis_jaws$overlap, main = "jaws", ylab = "tissues expressed")
#tis_lm <- glm(tis_jaws$freq.y~tis_jaws$overlap, family = poisson)
#tis_lm.resi <- resid(tis_lm)
#qqnorm(tis_lm.resi)
#qqline(tis_lm.resi)
#summary (tis_lm)
##write.csv(tis_jaws, "C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/pleiotropy/tis_jaws.csv") 
#
#embryos$Gene.name <- embryos$zeb_homolog
#head(embryos)
#tis_embryos <- merge(embryos, tis_counts, all.x = TRUE, by = c("Gene.name"))
#head(tis_embryos)
#tis_embryos   <- subset(tis_embryos    , !duplicated(tis_embryos[,1]))
#tis_embryos <- tis_embryos[order(tis_embryos$freq.y, decreasing = TRUE),]
#head(tis_embryos)
#nrow(tis_embryos[which(tis_embryos$overlap == "shared"),])
#nrow(tis_embryos[which(tis_embryos$overlap == "unshared"),])
#boxplot(tis_embryos$freq.y~tis_embryos$overlap, main = "embryos", ylab = "tissues expressed")
#stripchart(tis_embryos$freq.y~tis_embryos$overlap, 
#           vertical = TRUE, method = "jitter", 
#           pch = 21, 
#           add = TRUE) 
#tis_lm <- glm(tis_embryos$freq.y~tis_embryos$overlap, family = poisson)
#tis_lm.resi <- resid(tis_lm)
#qqnorm(tis_lm.resi)
#qqline(tis_lm.resi)
#summary (tis_lm)
#abline(tis_lm)
##write.csv(tis_embryos, "C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/pleiotropy/tis_embryos.csv") 

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
#unique(tis_test$Developmental.stage.name)
tis <- tis2
head(tis)
colnames(tis2)


embryos <- read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/pleiotropy/embryos.csv", header = TRUE, stringsAsFactors = FALSE)
jaws <- read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/pleiotropy/jaws.csv", header = TRUE, stringsAsFactors = FALSE)
tis_counts <- count(tis, 'Gene.name')
nrow(tis_counts)
head(tis_counts)
jaws$Gene.name <- jaws$zeb_homolog
head(jaws)
tis_jaws <- merge(jaws, tis_counts, all.x = TRUE, by = c("Gene.name"))
head(tis_jaws)
tis_jaws   <- subset(tis_jaws    , !duplicated(tis_jaws[,1]))
nrow(tis_jaws[which(tis_jaws$overlap == "shared"),])
nrow(tis_jaws[which(tis_jaws$overlap == "unshared"),])
boxplot(tis_jaws$freq.y~tis_jaws$overlap, main = "jaws", ylab = "tissues expressed")
tis_lm <- glm(tis_jaws$freq.y~tis_jaws$overlap, family = poisson)
stripchart(tis_jaws$freq.y~tis_jaws$overlap, 
           vertical = TRUE, method = "jitter", 
           pch = 21, 
           add = TRUE)
tis_lm.resi <- resid(tis_lm)
qqnorm(tis_lm.resi)
qqline(tis_lm.resi)
summary (tis_lm)
#write.csv(tis_jaws, "C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/pleiotropy/tis2_jaws.csv") 


embryos$Gene.name <- embryos$zeb_homolog
head(embryos)
tis_embryos <- merge(embryos, tis_counts, all.x = TRUE, by = c("Gene.name"))
head(tis_embryos)
tis_embryos   <- subset(tis_embryos    , !duplicated(tis_embryos[,1]))
tis_embryos <- tis_embryos[order(tis_embryos$freq.y, decreasing = TRUE),]
head(tis_embryos)
nrow(tis_embryos[which(tis_embryos$overlap == "shared"),])
nrow(tis_embryos[which(tis_embryos$overlap == "unshared"),])
boxplot(tis_embryos$freq.y~tis_embryos$overlap, main = "embryos", ylab = "tissues expressed")
stripchart(tis_embryos$freq.y~tis_embryos$overlap, 
           vertical = TRUE, method = "jitter", 
           pch = 21, 
           add = TRUE) 
#tis_lm <- glm(tis_embryos$freq.y~tis_embryos$overlap, family = poisson)
#tis_lm <- glm(tis_embryos$freq.y~tis_embryos$overlap, family = quasipoisson)
tis_lm <- glm.nb(tis_embryos$freq.y~tis_embryos$overlap)
dispersiontest(tis_lm)
tis_lm.resi <- resid(tis_lm)
qqnorm(tis_lm.resi)
qqline(tis_lm.resi)
summary (tis_lm)
abline(tis_lm)
#write.csv(tis_embryos, "C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/pleiotropy/tis2_embryos.csv") 



##### duplicated genes #####
dup <- read.table("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/pleiotropy/dgd_Danio_all_v71.tsv", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
embryos <- read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/pleiotropy/embryos.csv", header = TRUE, stringsAsFactors = FALSE)
jaws <- read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/pleiotropy/jaws.csv", header = TRUE, stringsAsFactors = FALSE)

head(dup)
head(embryos)

dup$zeb_homolog <- dup$Name


dup_e <- merge(dup, embryos, by = c("zeb_homolog"))
dup_j <- merge(dup, jaws, by = c("zeb_homolog"))
dup_e$duplicated <- "dup"
dup_j$duplicated <- "dup"
head(dup_e)
keeps <- c("zeb_homolog", "duplicated")
dup_e <- dup_e[keeps]
dup_j <- dup_j[keeps]

embryos <- merge(embryos, dup_e, all.x = TRUE, by = c("zeb_homolog"))
jaws <- merge(jaws, dup_e, all.x = TRUE, by = c("zeb_homolog"))
head(embryos)
head(jaws)

num_shared_dup_e   <- nrow(embryos[which(embryos$overlap == "shared" & embryos$duplicated == "dup"),])
num_shared_dup_j   <- nrow(jaws[which(jaws$overlap == "shared" & jaws$duplicated == "dup"),])
num_unshared_dup_e <- nrow(embryos[which(embryos$overlap == "unshared" & embryos$duplicated == "dup"),])
num_unshared_dup_j <- nrow(jaws[which(jaws$overlap == "unshared" & jaws$duplicated == "dup"),])
total_shared_e   <- nrow(embryos[which(embryos$overlap == "shared"),])
total_shared_j   <- nrow(jaws[which(jaws$overlap == "shared"),])
total_unshared_e <- nrow(embryos[which(embryos$overlap == "unshared"),])
total_unshared_j <- nrow(jaws[which(jaws$overlap == "unshared"),])
#write.csv(jaws, "C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/pleiotropy/jaws_dups.csv") 
## change duplicated column to 1 for unduplicated and 2 for duplicated in excel
embryos <- read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/pleiotropy/embryos_dups.csv", header = TRUE, stringsAsFactors = FALSE)
jaws <- read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/pleiotropy/jaws_dups.csv", header = TRUE, stringsAsFactors = FALSE)
table(embryos$overlap, embryos$duplicated)
table(jaws$overlap, jaws$duplicated)
chisq.test(embryos$overlap, embryos$duplicated, correct=FALSE)
chisq.test(jaws$overlap, jaws$duplicated, correct=FALSE)

num_shared_dup_e   / total_unshared_e   
num_unshared_dup_e / total_unshared_e 
num_shared_dup_j   / total_shared_j 
num_unshared_dup_j / total_unshared_j 



##### skeletal ontology #####
co <- read.table("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/candidates/craniofacial_skeletal_system_development_genes.txt", header = FALSE, stringsAsFactors = FALSE, sep = "\t")
cans <- read.table("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/candidates/zfin_gene_ids_de_genes_with_fixed_snps.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
sk <- read.table("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/candidates/skeletal_system_development_genes.txt", header = FALSE, stringsAsFactors = FALSE, sep = "\t")
mu <- read.table("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/candidates/muscle_organ_development_genes.txt", header = FALSE, stringsAsFactors = FALSE, sep = "\t")


co <- cbind(co, colsplit(co$V1, ":", c("junk", "ZFIN.ID")))
sk <- cbind(sk, colsplit(sk$V1, ":", c("junk", "ZFIN.ID")))
mu <- cbind(mu, colsplit(mu$V1, ":", c("junk", "ZFIN.ID")))
head(sk,20)
head(co)
head(cans)

merge(cans, co, by = c("ZFIN.ID"))
merge(cans, sk, by = c("ZFIN.ID"))
merge(cans, mu, by = c("ZFIN.ID"))



##### fst 0.8 gene overlap ####
#head(all_fst)
##snps <- read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/fst 0.8/axm_embryo_.8_overlap.csv", header = TRUE, stringsAsFactors = FALSE)
#head(snps)
#snps["CHROM"] <- snps$CHROM_x
#snps["POS"] <- snps["snp_ps"]
#test <- merge(snps, all_fst, by = c("CHROM", "POS"))
#test["snp_tag"] <- paste(test$CHROM, test$POS, sep = ":")
#length(unique(test$snp_tag))
#length(unique(test$Geneid))
#head(test)
##range(test$vxm_fst)
##nrow(all_fst[which(all_fst$vxm_fst >= 0.8),])
##write.csv(test, "C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/fst 0.8/axm_embryo_.8_overlap_dedup.csv") 
all_fst <- read.table("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/43_pupfish/fst/43_nomaysim/fst_43_fixed_comparisons", header = TRUE, stringsAsFactor = FALSE)
head(all_fst)
nrow(all_fst[which(all_fst$vxr >= 0.8),])

ap_embryo <- read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/fst 0.8/axp_embryo_.8_overlap_dedup.csv", header = TRUE, stringsAsFactors = FALSE, fileEncoding="UTF-8-BOM")
ap_jaw    <-   read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/fst 0.8/axp_jaw_.8_overlap_dedup.csv", header = TRUE, stringsAsFactors = FALSE, fileEncoding="UTF-8-BOM")
am_embryo <- read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/fst 0.8/axm_embryo_.8_overlap_dedup.csv", header = TRUE, stringsAsFactors = FALSE, fileEncoding="UTF-8-BOM")
am_jaw    <-   read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/fst 0.8/axm_jaw_.8_overlap_dedup.csv", header = TRUE, stringsAsFactors = FALSE, fileEncoding="UTF-8-BOM")
length(unique(ap_embryo$Geneid)) 
length(unique(ap_jaw$Geneid   )) 
length(unique(am_embryo$Geneid)) 
length(unique(am_jaw$Geneid))
length(unique(ap_embryo$snp_tag)) 
length(unique(ap_jaw$snp_tag)) 
length(unique(am_embryo$snp_tag)) 
length(unique(am_jaw$snp_tag))

p_snps <- c(ap_embryo$snp_tag,ap_jaw$snp_tag)
length(unique(p_snps))
m_snps <- c(am_embryo$snp_tag,am_jaw$snp_tag)
length(unique(m_snps))

p_genes <- c(ap_embryo$Geneid,ap_jaw$Geneid)
length(unique(p_genes))
m_genes <- c(am_embryo$Geneid,am_jaw$Geneid)
length(unique(m_genes))

ap_embryo <-ap_embryo[which(ap_embryo$vxr_fst == 1),] 
ap_jaw    <-   ap_jaw[which(ap_jaw$vxr_fst == 1),]    
am_embryo <-am_embryo[which(am_embryo$vxm_fst == 1),] 
am_jaw    <-   am_jaw[which(am_jaw$vxm_fst == 1),]    

p_snps <- c(ap_embryo$snp_tag,ap_jaw$snp_tag)
length(unique(p_snps))
m_snps <- c(am_embryo$snp_tag,am_jaw$snp_tag)
length(unique(m_snps))

p_genes <- c(ap_embryo$Geneid,ap_jaw$Geneid)
length(unique(p_genes))
m_genes <- c(am_embryo$Geneid,am_jaw$Geneid)
length(unique(m_genes))




shared_em <- merge(ap_embryo, am_embryo, by = c("Geneid"))
shared_jaw <- merge(ap_jaw, am_jaw, by = c("Geneid"))
length(unique(shared_jaw$Geneid   ))
length(unique(shared_em$Geneid   ))
Geneid <- unique(shared_em$Geneid)
shared_genes_em <- as.data.frame(Geneid, col.names= "Geneid")
Geneid <- unique(shared_jaw$Geneid)
shared_genes_jaw <- as.data.frame(Geneid, col.names="Geneid")

rna_seq_genes <- read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/blast/all_DE_master.csv", header = TRUE, stringsAsFactors = FALSE)
shared_em <- merge(shared_genes_em, rna_seq_genes, by = c("Geneid"))
shared_jaw <- merge(shared_genes_jaw, rna_seq_genes, by = c("Geneid"))
length(unique(shared_jaw$Geneid   ))
length(unique(shared_em$Geneid   ))

fixed_ap_embryo <- ap_embryo[which(ap_embryo$vxr_fst == 1),]
head(fixed_ap_embryo)

#write.csv(shared_em,"C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/fst 0.8/shared_embryo_genes.csv")


snps <- all_fst[which(all_fst$vxm_fst >= 0.8 & all_fst$vxm_fst <= 1),]
#snps <- read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/43_pupfish/fst/43_nomaysim/fixed_snp_genes_fixed_comparison_mp.csv", header = TRUE, stringsAsFactors = FALSE)
rna_seq_genes <- read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/mapping_correction/DESeq_genes_sig.csv", header = TRUE, stringsAsFactors = FALSE)
#write.csv(rna_seq_genes,"C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/mapping_correction/DESeq_genes_sig.csv")
rna_seq_genes <- merge(rna_seq_genes,scaffolds_key, by = c("Chr"))
rna_seq_genes <- rna_seq_genes[which(rna_seq_genes$comp == "ap_jaw"),]
#rna_seq_genes$comp == "ap_embryo"),]
#snps <- snps[which(snps$comp =="vxm"),]
head(rna_seq_genes,100) 
unique(rna_seq_genes$seq_type)
head(snps)
tail(snps)
head(scaffolds_key)
snps <- merge(snps, scaffolds_key, by = c("CHROM"))
test <- merge(rna_seq_genes, snps, by = c("CHROM"))
tail(test)


scaffs <- snps$Chr

scaffolds   <- c()
start_s     <- c()
stop_s      <- c()
gene_names <- c()
snp_pss     <- c()
features   <- c()

#test <- snps[which(snps$Chr =="NW_015150453.1"),]
#scaff <- "NW_015150453.1"

for (scaff in scaffs) 
{
  
  snp_table   <- snps[which(snps$Chr == scaff),]
  genes_table <- rna_seq_genes[which(rna_seq_genes$Chr == scaff),]
  gene_strt   <- genes_table$Start
  gene_stp    <- genes_table$End
  gene_name   <- genes_table$Geneid
  
  run <- c(1:length(gene_strt))
  if ((nrow(genes_table)) > 0 )
  {
    for (ps in snp_table$POS)
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
          
        }
        else if (ps >= (gene_strt[i] - 10000) & ps <= gene_stp[i] + 10000)
        {
          
          scaffolds <- c(scaffolds, scaff)
          start_s <- c(start_s, gene_strt[i])
          stop_s <- c(stop_s, gene_stp[i])
          gene_names <- c(gene_names, gene_name[i])
          snp_pss <- c(snp_pss,  ps)#paste(ps, "10kb", sep=':'))
          
        }
      }}}}

test5 <- data.frame(Chr = scaffolds, feature_start = start_s, feature_stop = stop_s, snp_ps = snp_pss, 
                    Geneid = gene_names, 
                    stringsAsFactors = FALSE)
head(test5)
nrow(test5)
#write.csv(test7, "C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/fst 0.8/fixed_m_snp_genes_am_jaw.csv") 
test6 <- merge(test5, rna_seq_genes, by = c("Geneid"))
head(test6)
test6$tag <- paste(test6$Chr.x, test6$snp_ps, test6$Geneid, sep=":")
head(test6)
length(test6)
length(unique(test6$tag))
length(unique(test6$Geneid))
test7 <- subset(test6, !duplicated(test6[,40]))
nrow(test7)





##### sweep stats #####
#fixed snps come from fst >= 0.8 above
ap_embryo <-ap_embryo[which(ap_embryo$vxr_fst == 1),] 
ap_jaw    <-   ap_jaw[which(ap_jaw$vxr_fst == 1),]    
am_embryo <-am_embryo[which(am_embryo$vxm_fst == 1),] 
am_jaw    <-   am_jaw[which(am_jaw$vxm_fst == 1),]
snp_de_genes <- rbind(ap_embryo, 
                      ap_jaw,    
                      am_embryo, 
                      am_jaw)
#snp_de_genes <- read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/mapping_correction/fixed_snps_near_de_genes_all_genes.csv")
#r_sweeps <- read.table("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/Jaw_length_candidate_project/recombined/Sweed/r_sweeps.txt", header = TRUE, stringsAsFactors = FALSE)
#r_sweeps <- cbind(r_sweeps, colsplit(r_sweeps$Position, "\\.", c("ps", "junk")))
#m_sweeps <- read.table("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/43_pupfish/SweeD/m_sweeps_43.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
#m_sweeps <- cbind(m_sweeps, colsplit(m_sweeps$Position, "\\.", c("ps", "junk")))
r_sweeps <- read.table("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/43_pupfish/SweeD/r_sweeps_bottle_43.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
r_sweeps <- cbind(r_sweeps, colsplit(r_sweeps$Position, "\\.", c("ps", "junk")))
m_sweeps <- read.table("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/43_pupfish/SweeD/m_sweeps_bottle_43.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
m_sweeps <- cbind(m_sweeps, colsplit(m_sweeps$Position, "\\.", c("ps", "junk")))

#snp_de_genes <- merge(snp_de_genes, scaffolds_key, by =c("Chr"))
length(snp_de_genes)
snp_de_genes    <- subset(snp_de_genes    , !duplicated(snp_de_genes[,45]))
head(snp_de_genes)
nrow(snp_de_genes)
head(m_sweeps)
length(unique(snp_de_genes$Geneid))

scaffs <- sort(unique(snp_de_genes$CHROM))
chrs <- c()
positions <- c()
quantiles_r <- c()
sweed_sigs <- c()
comparisons <- c()

#snp_de_genes <- test7

for (scaff in scaffs) 
{
  #par(mfrow = c(2,1))
  snps <- snp_de_genes[which(snp_de_genes$CHROM == scaff),]
  snp_list <- snps$snp_ps
  for (snp in snp_list)
  {
    
    sweep <- r_sweeps[which(r_sweeps$Scaffold == scaff),]
    max <- quantile(sweep$Likelihood, 1)
    min <- quantile(sweep$Likelihood, 0)
    sweep$norm <- ((sweep$Likelihood) - min) / (max - min) 
    sweed_sig <- sweep[which(sweep$ps >= snp),]
    sweed_sig <- sweed_sig$norm
    sweed_sig <- sweed_sig[1]
    sweed_sigs <- c(sweed_sigs, sweed_sig)
    quantile_r <- quantile(sweep$norm, .95)
    quantiles_r <- c(quantiles_r, quantile_r)
    positions <- c(positions, snp)
    chrs <- c(chrs, scaff)
    
  }
}

test5 <- data.frame(CHROM = chrs, snp_ps = positions, CLR_95_r = sweed_sigs, 
                    CLR_snp_window_r = quantiles_r, stringsAsFactors = FALSE)
test6 <- merge(snp_de_genes, test5, all.x = TRUE, by = c("CHROM", "snp_ps"))
head(test5)
head(test6)
nrow(test6)
length(test6)
count_genes <- count(test6, "Geneid")
freqs <- merge(test6, count_genes, all.x = TRUE, by=c("Geneid"))
nrow(freqs)
length(freqs)
head(freqs)
length(unique(freqs$Geneid))

#write.csv(freqs, "C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/sweeps/sweeps_bottle_r_m_de_genes_snps.csv") 

genes <- read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/sweeps/sweeps_bottle_slim_genes_snps.csv", header = TRUE, stringsAsFactors = FALSE)
all_DE <- read.csv("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/pleiotropy/all_DE_master_with_ensembl.csv", header = TRUE, stringsAsFactors = FALSE)
mbe_cans <- read.table("C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/mbe_candidates.txt", header = TRUE, stringsAsFactors = FALSE)

head(mbe_cans)
head(all_DE)
test <- merge(genes, all_DE, by = c("Geneid"))
head(test)
nrow(test)
length(test)
test    <- subset(test    , !duplicated(test[,1]))
test2 <- merge(test, mbe_cans, by = c("Geneid"))
unique(test2$symbol.x)
#write.csv(test , "C:/Users/Joseph McGirr Lab/Documents/Martin Lab/RNA-seq/sweeps/sweeps_bottle_slim__zeb_genes.csv") 
range(all_DE$baseMean)
median(all_DE$baseMean)
range((abs(all_DE$log2FoldChange)))
quantile(abs(all_DE$log2FoldChange))
median(all_DE$log2FoldChange)
quantile(abs(all_DE$log2FoldChange), .4)
median(abs(all_DE$log2FoldChange))



