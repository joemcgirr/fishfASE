######## ase jam way ########## 
###############################

######## install packages #####
###############################
#install.packages('VennDiagram')
#install.packages("ggplot2")
#install.packages("reshape2")
#install.packages("seqinr")
#install.packages("plyr")
#install.packages("MASS")
#install.packages("AER")
#install.packages("pheatmap")
#install.packages("vegan")
#install.packages("ape")
#install.packages("rgl")
source("http://bioconductor.org/biocLite.R")
biocLite("DESeq2")
biocLite("IRanges")
biocLite("stringi")
biocLite("XML")
biocLite("survival")
biocLite("foreign")
##install.packages("rlang")
#source("https://bioconductor.org/biocLite.R")
#pkgs <- rownames(installed.packages())
#biocLite(pkgs, type="source")
#install.packages("DESeq2")
#####
##### load libraries ##########
###############################

library(DESeq2)
library(reshape2)
library(seqinr)
library(plyr)
library(MASS)
library(AER)
library(rlang)

red <- "#E8000B" 
blu <- "#023EFF" 
yel <- "#FFC400"
bla <- "#000000" 
pur <- "#8B2BE2"
grb <- "#00D7FF" 
gre <- "#1AC938" 
dkr <- "#9F4800" 
whi <- "#ffffff"
lir <- "#F14CC1"


#####
###############################################################################################
###############################################################################################
##### perform DE analyses for comparisons. Output MA plots, gene tables, and summary stats ####
###############################################################################################
###############################################################################################

# compare hybrids
all_comps <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/ase_data/ase_comparisons.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")

head(all_comps)
final_features <- read.table("D:/Martin Lab/rna_2018/all_2018_samples/features_gff.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t", quote = "")
# mrna counts
all_cts <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/ase_data/all_samples_2018_counts_mrna_saf.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")

blast_key <- read.table("D:/Martin Lab/rna_2018/all_2018_samples/blast/cyprinodon_to_danio_one_way_best_hit_symbols.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
head(blast_key)
#i <- 84
#subset_comps <- c(31:34)

sum_stats <- data.frame(pop1=character(), 
                        pop2=character(),
                        stage=character(),
                        n_pop1=numeric(),
                        n_pop2=numeric(),
                        informative_transcripts=numeric(),
                        DE=numeric(),
                        up=numeric(),
                        down=numeric(),
                        #bs_1000_10q=numeric(),
                        #bs_1000_90q=numeric(),
                        #ds_3v3_10q=numeric(),
                        #ds_3v3_90q=numeric(),
                        stringsAsFactors=FALSE)

#for (i in (1:nrow(all_comps)))
for (i in subset_comps)
{
  #compare parent vs parent
  pops1_name <- all_comps$parent_pops1[i]
  pops2_name <- all_comps$parent_pops2[i]
  
  #compare parents vs hybrids
  #pops1_a <- all_comps$parent_pops1[i]
  #pops1_b <- all_comps$parent_pops2[i]
  #pops1_name <- paste(pops1_a, pops1_b, sep = "_and_")
  #pops2_name <- all_comps$hybrids[i]

  #compare parent1 vs hybrids
  #pops1_name <- all_comps$parent_pops1[i]
  #pops2_name <- all_comps$hybrids[i]
  
  #compare parent2 vs hybrids
  #pops1_name <- all_comps$parent_pops2[i]
  #pops2_name <- all_comps$hybrids[i]
  
  pops1_t <- pops1_name
  pops2_t <- pops2_name
  pops1 <- strsplit(pops1_name,"_and_")[[1]]
  pops2 <- strsplit(pops2_name,"_and_")[[1]]
  stage <- all_comps$stage[i]
  comp_name <- paste("condition_species_",pops1_name, "_vs_", pops2_name, "_" ,stage,".txt", sep = "")
  single_run <- paste(pops1_name, "_vs_", pops2_name, "_" ,stage, sep = "")
  
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
  write.table(comp_table, comp_name, row.names = FALSE, quote= FALSE,sep="\t")
  
  comp_file <- single_run
  
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
  #length(idx)
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
  res$padj <- p.adjust(res$pvalue, method="BH")
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
  
  #### bootstraps and downsampling ###
  #boot_10q_t <- 999
  #boot_90q_t <- 999
  #dns_10q_t <-  999
  #dns_90q_t <-  999
  #if (sample_size_a >2 & sample_size_b >2)
  #{
  #  
  #boots_in <- paste("D:/Martin Lab/rna_2018/all_2018_samples/bootstrapping_and_downsampling/", comp_file, "_boots.txt", sep = "")
  #dns_in    <- paste("D:/Martin Lab/rna_2018/all_2018_samples/bootstrapping_and_downsampling/", comp_file, "_downsamples.txt", sep = "")
  #boots <- read.table(boots_in, header = TRUE, stringsAsFactors = FALSE)
  #dns   <- read.table(dns_in, header = TRUE, stringsAsFactors = FALSE)
  #
  #boot_ci <- paste("BS ",round(as.numeric(quantile(boots$prop_DE,.1)[1]),digits = 4)*100,"-",round(as.numeric(quantile(boots$prop_DE,.9)[1]),digits = 4)*100,"%" , sep = "")
  #dns_ci  <- paste("DS ",round(as.numeric(quantile(dns$prop_DE,.1)[1]),digits = 4)*100,"-",round(as.numeric(quantile(dns$prop_DE,.9)[1]),digits = 4)*100,"%" , sep = "")
  #boot_10q_t <- round(as.numeric(quantile(boots$prop_DE,.1)[1]),digits = 4)*100
  #boot_90q_t <- round(as.numeric(quantile(boots$prop_DE,.9)[1]),digits = 4)*100
  #dns_10q_t <- round(as.numeric(quantile(dns$prop_DE,.1)[1]),digits = 4)*100
  #dns_90q_t <- round(as.numeric(quantile(dns$prop_DE,.9)[1]),digits = 4)*100
  #
  #}
  #
  i_stats <- data.frame(pop1=pops1_t, 
                        pop2=pops2_t,
                        stage=stage,
                        n_pop1=sample_size_a,
                        n_pop2=sample_size_b,
                        informative_transcripts=total_genes,
                        DE=de_total,
                        up=nrow(res_ordered[which(res_ordered$log2FoldChange < 0 & res_ordered$padj < 0.05),]),
                        down=nrow(res_ordered[which(res_ordered$log2FoldChange > 0 & res_ordered$padj < 0.05),]),
                        #bs_1000_10q=boot_10q_t,
                        #bs_1000_90q=boot_90q_t,
                        #ds_3v3_10q=dns_10q_t,
                        #ds_3v3_90q=dns_90q_t,
                        stringsAsFactors=FALSE)
  sum_stats <- rbind(sum_stats, i_stats)
  total_genes_plot <- paste(total_genes, "transcripts", sep = " ")
  de_total_plot <- paste(de_total, "DE", sep = " ")
  prop_de <- paste((100*(round(prop_de, digits = 3))), "% DE", sep = "")
  plotMA(resLFC, ylim=c(-6,5), main = comp_file)
  legend("topright", legend=c(sample_sizes_plot,total_genes_plot, de_total_plot, prop_de),cex=1.0, bty = 'n')
  legend("bottomleft", legend=c(paste((100*(round(de_dn, digits = 3))), "% DE down", sep = "")),cex=0.8, bty = 'n')
  legend("topleft", legend=c(paste((100*(round(de_up, digits = 3))), "% DE up", sep = "")),cex=0.8, bty = 'n')
  #legend('bottomright', legend=c(boot_ci, dns_ci),cex=0.8, bty = 'n')
  de_plot <- paste("C:/Users/jmcgirr/Documents/all_2018_samples/de_plots/",comp_file, "_de_plot.tiff", sep = "")
  tiff(de_plot, width = 7, height = 6, units = 'in', res = 600)
  plotMA(resLFC, ylim=c(-6,5), main = comp_file)
  legend("topright", legend=c(sample_sizes_plot,total_genes_plot, de_total_plot, prop_de),cex=1.0, bty = 'n')
  legend("bottomleft", legend=c(paste((100*(round(de_dn, digits = 3))), "% DE down", sep = "")),cex=0.8, bty = 'n')
  legend("topleft", legend=c(paste((100*(round(de_up, digits = 3))), "% DE up", sep = "")),cex=0.8, bty = 'n')
  #legend('bottomright', legend=c(boot_ci, dns_ci),cex=0.8, bty = 'n')
  dev.off()
  
  
  #### overlap with genes ###
  
  write.csv(res_ordered, file= cts_data_genes)
  genes <- read.csv(cts_data_genes, header = TRUE, stringsAsFactors = FALSE)
  genes <- cbind(genes, colsplit(genes$Geneid, ";", c("related_accession", "gene_name")))
  final <- merge(genes, final_features, by = c("related_accession"))
  length(unique(final$related_accession))
  length(unique(final$related_accession)) / length(unique(res_ordered$Geneid))
  head(final)
  nrow(final)
  head(blast_key)
  final <- merge(final, blast_key,all.x = TRUE, by = c("product_accession"))
  final <- final[order(final$padj, decreasing = FALSE),]
  write.csv(final, genes_out, row.names = FALSE)
  
}
#write.table(sum_stats, "C:/Users/jmcgirr/Documents/all_2018_samples/summary_stats_mrna_ph.txt", row.names = FALSE, quote= FALSE,sep="\t")


#####
########################################
########################################
##### inheritance patterns #############
########################################

#### make plot like McManus 2010 4B

inheritance_comps <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/ase_data/ase_comparisons.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
DE_genes_dir <- "C:/Users/jmcgirr/Documents/all_2018_samples/conditions/"
inheritance_plots_dir <- "C:/Users/jmcgirr/Documents/all_2018_samples/inheritance_patterns/plots/"
#i <- 28
subset_comps <- c(9,13,24,28)

#for (i in subset_comps)
for (i in (1:nrow(inheritance_comps)))
{
  parent_pops1_name <- inheritance_comps$parent_pops1[i]
  parent_pops2_name <- inheritance_comps$parent_pops2[i]
  hy_name <- inheritance_comps$hybrids[i]
  stage <- inheritance_comps$stage[i]
  parent_pops_v_hy  <- paste(DE_genes_dir, "DE_", parent_pops1_name,"_and_",parent_pops2_name, "_vs_",hy_name, "_",stage,"_genes",".csv",sep = "")
  parent_pops1_v_hy <- paste(DE_genes_dir, "DE_", parent_pops1_name, "_vs_",hy_name, "_",stage,"_genes",".csv",sep = "")
  parent_pops2_v_hy <- paste(DE_genes_dir, "DE_", parent_pops2_name, "_vs_",hy_name, "_",stage,"_genes",".csv",sep = "")
  parent_pops1_v_parent_pops2 <- paste(DE_genes_dir, "DE_", parent_pops1_name, "_vs_",parent_pops2_name, "_",stage,"_genes",".csv",sep = "")
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
  
  #Hybrid inheritance was considered additive if gene expression was intermediate between
  #generalists and molluscivores with significant differential expression between generalists and
  #molluscivores. Inheritance was dominant if expression was intermediate between parental
  #species and hybrid expression was significantly different from one parent but not the other.
  #Genes showing misexpression in hybrids showed transgressive inheritance, where hybrid gene
  #expression was significantly higher (overdominant) or lower (underdominant) than parental
  #populations.
  
  total_genes <- nrow(all_genes)
  con <- all_genes[which(all_genes$padj_ph > 0.05 & all_genes$padj_p1h > 0.05 & all_genes$padj_p2h > 0.05 & all_genes$padj_p1p2 > 0.05),]
  con$inheritance <- "conserved"
  con$inheritance <- 1
  nrow(con)
  add <- all_genes[which(all_genes$padj_ph > 0.05 & all_genes$padj_p1p2 < 0.05),]
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
  
  #cols <- 	c("#f6c700", "#E77200", "#FF00CC","#000000", "#bd1e24", "#0067a7")
  cols <- 	c(yel, "black", lir,gre, red, blu)
  
  #con = yellow ="#f6c700"
  #add = orange = "#E77200"
  #p1_dom = pink = "#FF00CC"
  #p2_dom = black = "#000000"
  #over_dom = red = #bd1e24
  #under_dom = blue = #0067a7
  
  cols_in <- cols[inheritance$inheritance]
  x_name <- paste(parent_pops1_name, hy_name, sep = " vs. ")
  y_name <- paste(parent_pops2_name, hy_name, sep = " vs. ")
  
  tiff(outfile_plot, width = 5, height = 5, units = 'in', res = 1000)
  plot(inheritance$lfc_p1h, inheritance$lfc_p2h, pch =16, col = cols_in, 
       ylab = paste("log2 fold change", y_name),
       xlab = paste("log2 fold change", x_name),cex.axis=1.5,
       main = plot_title )
  abline(v=0, col = "black", lty = 3, lwd = 1.8)
  abline(h=0, col = "black", lty = 3, lwd = 1.8)
  legend("topleft", inset=0,
         c(prop_con,
           prop_add,
           prop_p1_dom,
           prop_p2_dom, 
           prop_overdom, 
          prop_underdom),
         cex = 0.9, fill=cols, horiz=FALSE, bty="n")
  dev.off()
  
  #parallel mp inheritance overlap 
  #length(pmp_o_8)
  #length(intersect(pmp_o_8,con$related_accession))
  #length(intersect(pmp_o_8,con$related_accession))/length(pmp_o_8)
  #length(intersect(pmp_o_8,add$related_accession))
  #length(intersect(pmp_o_8,add$related_accession))/length(pmp_o_8)
  #length(intersect(pmp_o_8,p1_dom$related_accession))
  #length(intersect(pmp_o_8,p1_dom$related_accession))/length(pmp_o_8)
  #length(intersect(pmp_o_8,p2_dom$related_accession))
  #length(intersect(pmp_o_8,p2_dom$related_accession))/length(pmp_o_8)
  #length(intersect(pmp_o_8,over_dom$related_accession))
  #length(intersect(pmp_o_8,over_dom$related_accession))/length(pmp_o_8)
  #length(intersect(pmp_o_8,under_dom$related_accession))
  #length(intersect(pmp_o_8,under_dom$related_accession))/length(pmp_o_8)
  
  
}  



#####
###########################################
###########################################
#### Clustering ###########################
###########################################
###########################################

library(DESeq2)
library(pheatmap)
library(vegan)
library(rgl)
library(ape)

#b <- read.table("D:/Martin Lab/rna_2018/all_2018_samples/all_samples_48hpf.txt", header = TRUE, stringsAsFactors = FALSE)
#b <- b[b$sample %in% a$sample, ]
#write.table(b,"C:/Users/jmcgirr/Documents/all_2018_samples/48hpf_outlier_rm.txt", quote = FALSE, row.names = FALSE)

setwd("C:/Users/jmcgirr/Documents/all_2018_samples/")
comp <- "8dpf_outlier_rm"
comp <- "48hpf_outlier_rm"
comp <- "all_outlier_rm"
comp_file <- paste(comp , ".txt", sep = "")
#comp_file <- "table_maker_master_outlier_rm.txt"
all_cts <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/ase_data/all_samples_2018_counts_mrna_saf.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
sample_list <- read.table(comp_file, header = TRUE, stringsAsFactors = FALSE)
keeps <- c("Geneid", sample_list$sample)
keeper <- sample_list$sample
cts <- all_cts[keeps]
cts_out <- paste(comp, "_delete.txt", sep = "")
write.table(cts, cts_out, row.names = FALSE, quote= FALSE,sep="\t")
cts <- as.matrix(read.table(cts_out ,sep = "\t",header = TRUE,row.names=1)) 
head(cts)
nrow(cts)


colData <- as.matrix(read.table(comp_file ,header = TRUE,row.names=1))
head(colData)
ncol(cts)
nrow(colData)

#if (length(unique(sample_list$sequencing_round)) > 1)
#{
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = colData,
                              design= ~sequencing_round+f1+stage)
#}else 
#{
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = colData,
                              design= ~f1)
#}

dds <- estimateSizeFactors(dds)
idx <- rowSums(counts(dds, normalized=TRUE) >= 10 ) >= 10
length(rowSums(counts(dds, normalized=TRUE) >= 10 )) >= 10
length(idx)
dds <- dds[idx,]
dds <- DESeq(dds)
norm_cts <- data.frame(counts(dds, normalized=TRUE))
head(norm_cts)
nrow(norm_cts)
norm_cts$Geneid <- rownames(norm_cts)

#write.table(norm_cts_48, "C:/Users/jmcgirr/Documents/all_2018_samples/norm_cts_48hpf_no_seq_round.txt", row.names = FALSE, quote= FALSE,sep="\t")
#norm_cts_48 <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/norm_cts_48hpf_no_seq_round.txt", header=TRUE, stringsAsFactors = FALSE)
#norm_cts_8 <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/norm_cts_8dpf.txt", header=TRUE, stringsAsFactors = FALSE)

head(norm_cts_48)


# PCA
library(ggplot2)

#rld <- rlog(dds)
rld <- vst(dds)
dists <- dist(t(assay(rld)))
plot(hclust(dists))

# contrast groups
results(dds, contrast=c("f1","CRPA","NCA"))


#plotPCA(rld, intgroup=c("stage"))

rld.sub <- rld[ , rld$stage %in% c("8dpf") ]
plotPCA(rld.sub, intgroup=c("f1"))+theme_bw()
#+geom_point(aes(shape='cross_type'))+scale_shape_manual(values=c(25))

pcaData <- plotPCA(rld.sub, intgroup=c("f1",'cross_type'), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
p1 <- ggplot(pcaData, aes(PC1, PC2, color=f1, shape=cross_type)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()+theme_bw()+scale_shape_manual(values=c(17, 16))
#tiff("C:/Users/jmcgirr/Documents/all_2018_samples/manuscript_figs/supp/pca_8dpf.tiff", width = 6, height = 6, units = 'in', res = 1000)
p1 
#dev.off()
#+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#panel.background = element_blank(), axis.line = element_line(colour = "black"))

# hierarchical clustering of samples and heatmap of sample similarities

vst=varianceStabilizingTransformation(dds)  
vsd=assay(vst)  # vsd is now the normalized log2-transformed data

#tiff("D:/Martin Lab/rna_2018/de_plots/pheatmap_round_1_48hpf.tiff", width = 12, height = 12, units = 'in', res = 800)
pheatmap(cor(vsd))  # yes it is that simple
#dev.off()
out <- pheatmap(cor(vsd)) 
#tiff("D:/Martin Lab/rna_2018/de_plots/dendro_round_1_48hpf.tiff", width = 10, height = 8, units = 'in', res = 800)
plot(out$tree_row)
#dev.off()



# Principle coordiante analysis

# assembling table of conditions to lable PCoA plot:
# (in the chunk below, replace factor1 and factor2 with your actual factor names from myConditions table)
factor1=as.character(colData(dds)$species)
factor2=as.character(colData(dds)$stage)
oneByTwo=paste(factor1,factor2,sep=".")
conditions=data.frame(cbind(factor1,factor2,oneByTwo))

# actual PCoA analysis
dds.pcoa=pcoa(vegdist(t(vsd),method="manhattan")/1000)
scores=dds.pcoa$vectors

# plotting
plot(scores[,1], scores[,2],col=as.numeric(as.factor(factor1)))
ordispider(scores,factor2,label=T)
ordiellipse(scores,factor2)

# interactive 3d plot - can rotate it by dragging mouse
radiusScale=2 # change it if the spheres come out too small or too big in the next one
plot3d(scores[,1], scores[,2],scores[,3],col=as.numeric(as.factor(factor1)),type="s",radius=radiusScale*as.numeric(as.factor(factor2)))

# formal permutation-based analysis of variance 
adonis(t(vsd)~factor1*factor2,data=conditions,method="manhattan") 

#####
######################################################
######################################################
##### uniquely DE genes ##############################
######################################################
####### 8dpf ####
DE_genes_dir <- "C:/Users/jmcgirr/Documents/all_2018_samples/conditions/"
venn_dir <- "C:/Users/jmcgirr/Documents/all_2018_samples/de_venns/"
stage <- "8dpf"

#load DE genes #
{
  # am crosses
  caxcm <- read.csv(paste(DE_genes_dir,"DE_CRPA_vs_CRPM_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  oaxom <- read.csv(paste(DE_genes_dir,"DE_OSPA_vs_OSPM_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  # ap crosses
  caxcp <- read.csv(paste(DE_genes_dir,"DE_CRPA_vs_CRPP_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  oaxop <- read.csv(paste(DE_genes_dir,"DE_OSPA_vs_OSPP_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  # specialist crosses
  cmxcp <- read.csv(paste(DE_genes_dir,"DE_CRPM_vs_CRPP_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  omxop <- read.csv(paste(DE_genes_dir,"DE_OSPM_vs_OSPP_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  # outgroup crosses
  caxna <-  read.csv(paste(DE_genes_dir,"DE_CRPA_vs_NCA_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  oaxna <-  read.csv(paste(DE_genes_dir,"DE_OSPA_vs_NCA_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  caxup <-  read.csv(paste(DE_genes_dir,"DE_CRPA_vs_UPxUA_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  oaxup <-  read.csv(paste(DE_genes_dir,"DE_OSPA_vs_UPxUA_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  
  # lake crosses
  caxoa <- read.csv(paste(DE_genes_dir,"DE_CRPA_vs_OSPA_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  cmxom <- read.csv(paste(DE_genes_dir,"DE_CRPM_vs_OSPM_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  cpxop <- read.csv(paste(DE_genes_dir,"DE_CRPP_vs_OSPP_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  
  # sig DE
  
  caxcm <- caxcm[which(caxcm$padj <= 0.05),] 
  oaxom <- oaxom[which(oaxom$padj <= 0.05),]
  caxcp <- caxcp[which(caxcp$padj <= 0.05),]
  oaxop <- oaxop[which(oaxop$padj <= 0.05),]
  cmxcp <- cmxcp[which(cmxcp$padj <= 0.05),]
  omxop <- omxop[which(omxop$padj <= 0.05),]
  caxna <- caxna[which(caxna$padj <= 0.05),]
  oaxna <- oaxna[which(oaxna$padj <= 0.05),]
  caxup <- caxup[which(caxup$padj <= 0.05),]
  oaxup <- oaxup[which(oaxup$padj <= 0.05),]
  caxoa <- caxoa[which(caxoa$padj <= 0.05),]
  cmxom <- cmxom[which(cmxom$padj <= 0.05),]
  cpxop <- cpxop[which(cpxop$padj <= 0.05),]
  
  caxcm <- caxcm$related_accession
  oaxom <- oaxom$related_accession
  caxcp <- caxcp$related_accession
  oaxop <- oaxop$related_accession
  cmxcp <- cmxcp$related_accession
  omxop <- omxop$related_accession
  
  caxna <- caxna$related_accession
  oaxna <- oaxna$related_accession
  caxup <- caxup$related_accession
  oaxup <- oaxup$related_accession
  caxoa <- caxoa$related_accession
  cmxom <- cmxom$related_accession
  cpxop <- cpxop$related_accession
}

tiff(paste(venn_dir, "mp_only_8dpf.tiff"), width = 5, height = 5, units = 'in', res = 300)
dev.off()
# all p
venn(list(caxcp=caxcp,oaxop=oaxop), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
all_ap <- intersect(caxcp,oaxop)
venn(list(all_ap=all_ap,caxna=caxna,oaxna=oaxna,caxup=caxup,oaxup=oaxup,caxoa=caxoa,cpxop=cpxop), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
ap_only <- setdiff(all_ap, c(caxna,oaxna,caxup,oaxup,caxoa,cpxop))
length(all_ap)
length(ap_only)
length(ap_only)-length(all_ap)
length(ap_only)/length(all_ap)

# all m
venn(list(caxcm=caxcm,oaxom=oaxom), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
all_am <- intersect(caxcm,oaxom)
venn(list(all_am=all_ap,caxna=caxna,oaxna=oaxna,caxup=caxup,oaxup=oaxup,caxoa=caxoa,cmxom=cmxom), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
am_only <- setdiff(all_ap, c(caxna,oaxna,caxup,oaxup,caxoa,cmxom))
length(all_am)
length(am_only)
length(am_only)-length(all_am)
length(am_only)/length(all_am)

# all mp
venn(list(cmxcp=cmxcp,omxop=omxop), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
all_mp <- intersect(cmxcp,omxop)
venn(list(all_mp=all_mp,caxna=caxna,oaxna=oaxna,caxup=caxup,oaxup=oaxup,caxoa=caxoa,cmxom=cmxom,cpxop=cpxop), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
mp_only <- setdiff(all_mp, c(caxna,oaxna,caxup,oaxup,caxoa,cmxom,cpxop))
length(all_mp)
length(mp_only)
length(mp_only)-length(all_mp)
length(mp_only)/length(all_mp)

cis_genes_specialist_comps <- c("XM_015371309.1","XM_015372239.1","XM_015375549.1","XM_015383575.1","XM_015396813.1","XM_015398422.1","XM_015399203.1","XM_015399204.1","XM_015373937.1","XM_015375704.1","XM_015379828.1","XM_015381009.1","XM_015386924.1","XM_015386952.1","XM_015392732.1","XM_015395494.1","XM_015402407.1","XM_015404979.1","XM_015371309.1","XM_015372239.1","XM_015390653.1","XM_015392922.1","XM_015396691.1","XM_015371217.1","XM_015371319.1","XM_015386908.1","XM_015371364.1","XM_015375868.1","XM_015380501.1")
intersect(ap_only,cis_genes_specialist_comps)
intersect(am_only,cis_genes_specialist_comps)
intersect(mp_only,cis_genes_specialist_comps)
st <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/ase_data/ase_and_mse/strict_ase_unphased_mbased_transtest_by_individual.txt", stringsAsFactors = FALSE, header = TRUE, sep = "\t")
ap_48 <- st[c(8,12,15),]
am_48 <- st[c(7,10,11,14),]
mp_48 <- st[c(9,13),]
ap_8  <- st[c(23,27,29),]
am_8  <- st[c(22,25,26),]
mp_8  <- st[c(24,28),]
cool_gene <- "XM_015399204.1"

for(cool_gene in mp_only)
{
  hits <- mp_8[,c(1,2,3,(grep(cool_gene, mp_8)))]
  if (length(hits) >3)
  {
    print(cool_gene)
    print(hits)
  }
}
for(cool_gene in ap_only)
{
  hits <- ap_8[,c(1,2,3,(grep(cool_gene, ap_8)))]
  if (length(hits) >3)
  {
    print(cool_gene)
    print(hits)
  }
}
for(cool_gene in am_only)
{
  hits <- am_8[,c(1,2,3,(grep(cool_gene, am_8)))]
  if (length(hits) >3)
  {
    print(cool_gene)
    print(hits)
  }
}

####### 48hpf ####
DE_genes_dir <- "C:/Users/jmcgirr/Documents/all_2018_samples/conditions/"
venn_dir <- "C:/Users/jmcgirr/Documents/all_2018_samples/de_venns/"
stage <- "48hpf"

#load DE genes #
{
# am crosses
caxcm <- read.csv(paste(DE_genes_dir,"DE_CRPA_vs_CRPM_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
oaxom <- read.csv(paste(DE_genes_dir,"DE_OSPA_vs_OSPM_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
# ap crosses
caxcp <- read.csv(paste(DE_genes_dir,"DE_CRPA_vs_CRPP_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
oaxop <- read.csv(paste(DE_genes_dir,"DE_OSPA_vs_OSPP_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
# specialist crosses
cmxcp <- read.csv(paste(DE_genes_dir,"DE_CRPM_vs_CRPP_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
omxop <- read.csv(paste(DE_genes_dir,"DE_OSPM_vs_OSPP_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
# outgroup crosses
caxna <-  read.csv(paste(DE_genes_dir,"DE_CRPA_vs_NCA_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
oaxna <-  read.csv(paste(DE_genes_dir,"DE_OSPA_vs_NCA_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
caxup <-  read.csv(paste(DE_genes_dir,"DE_CRPA_vs_UPxUA_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
oaxup <-  read.csv(paste(DE_genes_dir,"DE_OSPA_vs_UPxUA_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)

# lake crosses
caxoa <- read.csv(paste(DE_genes_dir,"DE_CRPA_vs_OSPA_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
cmxom <- read.csv(paste(DE_genes_dir,"DE_CRPM_vs_OSPM_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
cpxop <- read.csv(paste(DE_genes_dir,"DE_CRPP_vs_OSPP_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)

# sig DE

caxcm <- caxcm[which(caxcm$padj <= 0.05),] 
oaxom <- oaxom[which(oaxom$padj <= 0.05),]
caxcp <- caxcp[which(caxcp$padj <= 0.05),]
oaxop <- oaxop[which(oaxop$padj <= 0.05),]
cmxcp <- cmxcp[which(cmxcp$padj <= 0.05),]
omxop <- omxop[which(omxop$padj <= 0.05),]
caxna <- caxna[which(caxna$padj <= 0.05),]
oaxna <- oaxna[which(oaxna$padj <= 0.05),]
caxup <- caxup[which(caxup$padj <= 0.05),]
oaxup <- oaxup[which(oaxup$padj <= 0.05),]
caxoa <- caxoa[which(caxoa$padj <= 0.05),]
cmxom <- cmxom[which(cmxom$padj <= 0.05),]
cpxop <- cpxop[which(cpxop$padj <= 0.05),]

caxcm <- caxcm$related_accession
oaxom <- oaxom$related_accession
caxcp <- caxcp$related_accession
oaxop <- oaxop$related_accession
cmxcp <- cmxcp$related_accession
omxop <- omxop$related_accession

caxna <- caxna$related_accession
oaxna <- oaxna$related_accession
caxup <- caxup$related_accession
oaxup <- oaxup$related_accession
caxoa <- caxoa$related_accession
cmxom <- cmxom$related_accession
cpxop <- cpxop$related_accession
}

tiff(paste(venn_dir, "mp_only_48hpf.tiff"), width = 5, height = 5, units = 'in', res = 300)
dev.off()
# all p
venn(list(caxcp=caxcp,oaxop=oaxop), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
all_ap <- intersect(caxcp,oaxop)
venn(list(all_ap=all_ap,caxna=caxna,oaxna=oaxna,caxup=caxup,oaxup=oaxup,caxoa=caxoa,cpxop=cpxop), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
ap_only <- setdiff(all_ap, c(caxna,oaxna,caxup,oaxup,caxoa,cpxop))
length(all_ap)
length(ap_only)
length(ap_only)-length(all_ap)
length(ap_only)/length(all_ap)

# all m
venn(list(caxcm=caxcm,oaxom=oaxom), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
all_am <- intersect(caxcm,oaxom)
venn(list(all_am=all_ap,caxna=caxna,oaxna=oaxna,caxup=caxup,oaxup=oaxup,caxoa=caxoa,cmxom=cmxom), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
am_only <- setdiff(all_ap, c(caxna,oaxna,caxup,oaxup,caxoa,cmxom))
length(all_am)
length(am_only)
length(am_only)-length(all_am)
length(am_only)/length(all_am)

# all mp
venn(list(cmxcp=cmxcp,omxop=omxop), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
all_mp <- intersect(cmxcp,omxop)
venn(list(all_mp=all_mp,caxna=caxna,oaxna=oaxna,caxup=caxup,oaxup=oaxup,caxoa=caxoa,cmxom=cmxom,cpxop=cpxop), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
mp_only <- setdiff(all_mp, c(caxna,oaxna,caxup,oaxup,caxoa,cmxom,cpxop))
length(all_mp)
length(mp_only)
length(mp_only)-length(all_mp)
length(mp_only)/length(all_mp)

for(cool_gene in mp_only)
{
  hits <- mp_48[,c(1,2,3,(grep(cool_gene, mp_48)))]
  if (length(hits) >3)
  {
    print(cool_gene)
    print(hits)
  }
}
for(cool_gene in ap_only)
{
  hits <- ap_48[,c(1,2,3,(grep(cool_gene, ap_48)))]
  if (length(hits) >3)
  {
    print(cool_gene)
    print(hits)
  }
}
for(cool_gene in am_only)
{
  hits <- am_48[,c(1,2,3,(grep(cool_gene, am_48)))]
  if (length(hits) >3)
  {
    print(cool_gene)
    print(hits)
  }
}


#####
######################################################
######################################################
##### uniquely misexpressed genes ####################
######################################################


library(venn)
DE_genes_dir <- "C:/Users/jmcgirr/Documents/all_2018_samples/conditions/"
venn_dir <- "C:/Users/jmcgirr/Documents/all_2018_samples/mse_venns/"
stage <- "8dpf"


# am crosses
caxcm <- read.csv(paste(DE_genes_dir,"DE_CRPA_and_CRPM_vs_CAxCM_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
cmxca <- read.csv(paste(DE_genes_dir,"DE_CRPM_and_CRPA_vs_CMxCA_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
oaxom <- read.csv(paste(DE_genes_dir,"DE_OSPA_and_OSPM_vs_OAxOM_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
omxoa <- read.csv(paste(DE_genes_dir,"DE_OSPM_and_OSPA_vs_OMxOA_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
# ap crosses
caxcp <- read.csv(paste(DE_genes_dir,"DE_CRPA_and_CRPP_vs_CAxCP_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
oaxop <- read.csv(paste(DE_genes_dir,"DE_OSPA_and_OSPP_vs_OAxOP_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
opxoa <- read.csv(paste(DE_genes_dir,"DE_OSPP_and_OSPA_vs_OPxOA_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
# specialist crosses
cmxcp <- read.csv(paste(DE_genes_dir,"DE_CRPM_and_CRPP_vs_CMxCP_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
opxom <- read.csv(paste(DE_genes_dir,"DE_OSPM_and_OSPP_vs_OPxOM_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
# outgroup crosses
naxca <-  read.csv(paste(DE_genes_dir,"DE_CRPA_and_NCA_vs_NAxCA_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
upxca <-  read.csv(paste(DE_genes_dir,"DE_CRPA_and_UPxUA_vs_UPxCA_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
# lake crosses
oaxca <- read.csv(paste(DE_genes_dir,"DE_CRPA_and_OSPA_vs_OAxCA_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
caxoa <- read.csv(paste(DE_genes_dir,"DE_CRPA_and_OSPA_vs_CAxOA_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
cmxom <- read.csv(paste(DE_genes_dir,"DE_CRPM_and_OSPM_vs_CMxOM_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
opxcp <- read.csv(paste(DE_genes_dir,"DE_CRPP_and_OSPP_vs_OPxCP_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)

# sig misexpressed
caxcm <- caxcm[which(caxcm$padj <= 0.05),] 
cmxca <- cmxca[which(cmxca$padj <= 0.05),] 
oaxom <- oaxom[which(oaxom$padj <= 0.05),] 
omxoa <- omxoa[which(omxoa$padj <= 0.05),] 
caxcp <- caxcp[which(caxcp$padj <= 0.05),] 
oaxop <- oaxop[which(oaxop$padj <= 0.05),] 
opxoa <- opxoa[which(opxoa$padj <= 0.05),] 
cmxcp <- cmxcp[which(cmxcp$padj <= 0.05),] 
opxom <- opxom[which(opxom$padj <= 0.05),] 
naxca <- naxca[which(naxca$padj <= 0.05),] 
upxca <- upxca[which(upxca$padj <= 0.05),] 
oaxca <- oaxca[which(oaxca$padj <= 0.05),] 
caxoa <- caxoa[which(caxoa$padj <= 0.05),] 
cmxom <- cmxom[which(cmxom$padj <= 0.05),] 
opxcp <- opxcp[which(opxcp$padj <= 0.05),] 

caxcm <- caxcm$related_accession
cmxca <- cmxca$related_accession
oaxom <- oaxom$related_accession
omxoa <- omxoa$related_accession

caxcp <- caxcp$related_accession
oaxop <- oaxop$related_accession
opxoa <- opxoa$related_accession

cmxcp <- cmxcp$related_accession
opxom <- opxom$related_accession

naxca <- naxca$related_accession
upxca <- upxca$related_accession

oaxca <- oaxca$related_accession
caxoa <- caxoa$related_accession
cmxom <- cmxom$related_accession
opxcp <- opxcp$related_accession


### by species
# m only
tiff("C:/Users/jmcgirr/Documents/all_2018_samples/mse_venns/mp_only_8dpf.tiff", width = 5, height = 5, units = 'in', res = 300)
dev.off()
venn(list(caxcm=caxcm,cmxca=cmxca,oaxom=oaxom,omxoa=omxoa), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
venn(list(caxcm=caxcm,cmxca=cmxca,oaxom=oaxom), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
m_mis <- intersect(intersect(caxcm,cmxca),oaxom)
# 48hpf
m_mis <- intersect(intersect(caxcm,cmxca),intersect(oaxom,omxoa))
venn(list(m_mis=m_mis,cmxom=cmxom,naxca=naxca,upxca=upxca), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)

m_only <- setdiff(intersect(intersect(caxcm,cmxca),oaxom),c(cmxom,naxca,upxca))
# 48hpf
m_only <- setdiff(intersect(intersect(caxcm,cmxca),intersect(oaxom,omxoa)),c(cmxom,naxca,upxca))

# p only
venn(list(caxcp=caxcp,oaxop=oaxop,opxoa=opxoa), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
p_mis <- intersect(intersect(caxcp,oaxop),opxoa)
venn(list(p_mis=p_mis,opxcp=opxcp,naxca=naxca,upxca=upxca), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
p_only <- setdiff(intersect(caxcp,intersect(opxoa,oaxop)),c(opxcp,naxca,upxca))

# mp only
venn(list(cmxcp =cmxcp,opxom=opxom), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
mp_mis <- intersect(cmxcp,opxom)
venn(list(mp_mis=mp_mis,opxcp=opxcp,naxca=naxca,upxca=upxca,cmxom=cmxom), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
mp_only <- setdiff(intersect(cmxcp,opxom),c(opxcp,naxca,upxca,cmxom))



### by lake
# m only
venn(list(caxcm=caxcm,cmxca=cmxca), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
m_mis <- intersect(caxcm,cmxca)
venn(list(m_mis=m_mis,cmxom=cmxom,naxca=naxca,upxca=upxca), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
m_only <- setdiff(intersect(caxcm,cmxca),c(cmxom,naxca,upxca))

venn(list(oaxom=oaxom,omxoa=omxoa), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
m_mis <- intersect(oaxom,omxoa)
venn(list(m_mis=m_mis,cmxom=cmxom,naxca=naxca,upxca=upxca), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
m_only <- setdiff(oaxom,c(cmxom,naxca,upxca))
#48hpf
m_only <- setdiff(intersect(oaxom,omxoa),c(cmxom,naxca,upxca))
# p_only
venn(list(oaxop=oaxop,opxoa=opxoa), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
p_mis <- intersect(oaxop,opxoa)
venn(list(p_mis=p_mis,opxcp=opxcp,naxca=naxca,upxca=upxca), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
p_only <- setdiff(p_mis,c(opxcp,naxca,upxca))
p_only <- setdiff(caxcp,c(opxcp,naxca,upxca))


##load txt files from GO analyses section
p_only_o_48  <- data.frame(related_accession =p_only  )
p_only_c_48  <- data.frame(related_accession =p_only  )
m_only_o_48  <- data.frame(related_accession =m_only  )
m_only_c_48  <- data.frame(related_accession =m_only  )

p_only_o_8   <- data.frame(related_accession =p_only  )
p_only_c_8   <- data.frame(related_accession =p_only  )
m_only_o_8   <- data.frame(related_accession =m_only  )
m_only_c_8   <- data.frame(related_accession =m_only  )

p_only_o_48<- merge(p_only_o_48,final_features, by = c("related_accession"))
p_only_c_48<- merge(p_only_c_48,final_features, by = c("related_accession"))
m_only_o_48<- merge(m_only_o_48,final_features, by = c("related_accession"))
m_only_c_48<- merge(m_only_c_48,final_features, by = c("related_accession"))
p_only_o_8 <- merge(p_only_o_8 ,final_features, by = c("related_accession"))
p_only_c_8 <- merge(p_only_c_8 ,final_features, by = c("related_accession"))
m_only_o_8 <- merge(m_only_o_8 ,final_features, by = c("related_accession"))
m_only_c_8 <- merge(m_only_c_8 ,final_features, by = c("related_accession"))

p_only_o_48<- merge(p_only_o_48, blast_key, by = c("product_accession"))
p_only_c_48<- merge(p_only_c_48, blast_key, by = c("product_accession"))
m_only_o_48<- merge(m_only_o_48, blast_key, by = c("product_accession"))
m_only_c_48<- merge(m_only_c_48, blast_key, by = c("product_accession"))
p_only_o_8 <- merge(p_only_o_8 , blast_key, by = c("product_accession"))
p_only_c_8 <- merge(p_only_c_8 , blast_key, by = c("product_accession"))
m_only_o_8 <- merge(m_only_o_8 , blast_key, by = c("product_accession"))
m_only_c_8 <- merge(m_only_c_8 , blast_key, by = c("product_accession"))

m_only_o_8$zeb_gene_symbol_one_way







# overlap with regulatory mechanisms
st <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/ase_data/ase_and_mse/strict_ase_unphased_mbased_transtest_by_individual.txt", stringsAsFactors = FALSE, header = TRUE, sep = "\t")
st2 <- st[which(st$stage_s == "8dpf"),]
st8 <- st[which(st$stage_s == "48hpf"),]

caxcm_2 <-  strsplit((st2[which(st2$hy_name_s == "CAxCM"),]$cis_genes),";")[[1]]  
cmxca_2 <-  strsplit((st2[which(st2$hy_name_s == "CMxCA"),]$cis_genes),";")[[1]] 
oaxom_2 <-  strsplit((st2[which(st2$hy_name_s == "OAxOM"),]$cis_genes),";")[[1]] 
omxoa_2 <-  strsplit((st2[which(st2$hy_name_s == "OMxOA"),]$cis_genes),";")[[1]] 
caxcp_2 <-  strsplit((st2[which(st2$hy_name_s == "CAxCP"),]$cis_genes),";")[[1]]
oaxop_2 <-  strsplit((st2[which(st2$hy_name_s == "OAxOP"),]$cis_genes),";")[[1]]
opxoa_2 <-  strsplit((st2[which(st2$hy_name_s == "OPxOA"),]$cis_genes),";")[[1]]
cmxcp_2 <-  strsplit((st2[which(st2$hy_name_s == "CMxCP"),]$cis_genes),";")[[1]]
opxom_2 <-  strsplit((st2[which(st2$hy_name_s == "OPxOM"),]$cis_genes),";")[[1]]
naxca_2 <-  strsplit((st2[which(st2$hy_name_s == "NAxCA"),]$cis_genes),";")[[1]]
upxca_2 <-  strsplit((st2[which(st2$hy_name_s == "UPxCA"),]$cis_genes),";")[[1]]
oaxca_2 <-  strsplit((st2[which(st2$hy_name_s == "OAxCA"),]$cis_genes),";")[[1]]
caxoa_2 <-  strsplit((st2[which(st2$hy_name_s == "CAxOA"),]$cis_genes),";")[[1]]
cmxom_2 <-  strsplit((st2[which(st2$hy_name_s == "CMxOM"),]$cis_genes),";")[[1]]
opxcp_2 <-  strsplit((st2[which(st2$hy_name_s == "OPxCP"),]$cis_genes),";")[[1]]

caxcm_8 <-  strsplit((st8[which(st8$hy_name_s == "CAxCM"),]$cis_genes),";")[[1]]  
cmxca_8 <-  strsplit((st8[which(st8$hy_name_s == "CMxCA"),]$cis_genes),";")[[1]] 
oaxom_8 <-  strsplit((st8[which(st8$hy_name_s == "OAxOM"),]$cis_genes),";")[[1]] 
omxoa_8 <-  strsplit((st8[which(st8$hy_name_s == "OMxOA"),]$cis_genes),";")[[1]] 
caxcp_8 <-  strsplit((st8[which(st8$hy_name_s == "CAxCP"),]$cis_genes),";")[[1]]
oaxop_8 <-  strsplit((st8[which(st8$hy_name_s == "OAxOP"),]$cis_genes),";")[[1]]
opxoa_8 <-  strsplit((st8[which(st8$hy_name_s == "OPxOA"),]$cis_genes),";")[[1]]
cmxcp_8 <-  strsplit((st8[which(st8$hy_name_s == "CMxCP"),]$cis_genes),";")[[1]]
opxom_8 <-  strsplit((st8[which(st8$hy_name_s == "OPxOM"),]$cis_genes),";")[[1]]
naxca_8 <-  strsplit((st8[which(st8$hy_name_s == "NAxCA"),]$cis_genes),";")[[1]]
upxca_8 <-  strsplit((st8[which(st8$hy_name_s == "UPxCA"),]$cis_genes),";")[[1]]
oaxca_8 <-  strsplit((st8[which(st8$hy_name_s == "OAxCA"),]$cis_genes),";")[[1]]
caxoa_8 <-  strsplit((st8[which(st8$hy_name_s == "CAxOA"),]$cis_genes),";")[[1]]
cmxom_8 <-  strsplit((st8[which(st8$hy_name_s == "CMxOM"),]$cis_genes),";")[[1]]
opxcp_8 <-  strsplit((st8[which(st8$hy_name_s == "OPxCP"),]$cis_genes),";")[[1]]

# by lake
intersect(caxcm_2,cmxca_2)
intersect(oaxom_2,omxoa_2)
intersect(opxoa_2,oaxop_2)
intersect(cmxcp_2,opxom_2)
intersect(caxcm_8,cmxca_8)
intersect(oaxom_8,omxoa_8)
intersect(opxoa_8,oaxop_8)
intersect(cmxcp_8,opxom_8)
intersect(caxcm_8,caxcm_2)
intersect(oaxom_8,oaxom_2)
intersect(opxoa_8,opxoa_2)
intersect(cmxcp_8,cmxcp_2)

st <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/ase_data/ase_and_mse/super_lenient_ase_unphased_mbased.txt", stringsAsFactors = FALSE, header = TRUE, sep = "\t")
stage <- "48hpf"
st <- st[which(st$stage_s == stage),]

#compensatory mse overlap
caxcm_g <-  strsplit((st[which(st$hy_name_s == "CAxCM"),]$comp_genes),";")[[1]]  
cmxca_g <-  strsplit((st[which(st$hy_name_s == "CMxCA"),]$comp_genes),";")[[1]] 
oaxom_g <-  strsplit((st[which(st$hy_name_s == "OAxOM"),]$comp_genes),";")[[1]] 
omxoa_g <-  strsplit((st[which(st$hy_name_s == "OMxOA"),]$comp_genes),";")[[1]] 
caxcp_g <-  strsplit((st[which(st$hy_name_s == "CAxCP"),]$comp_genes),";")[[1]]
oaxop_g <-  strsplit((st[which(st$hy_name_s == "OAxOP"),]$comp_genes),";")[[1]]
opxoa_g <-  strsplit((st[which(st$hy_name_s == "OPxOA"),]$comp_genes),";")[[1]]
cmxcp_g <-  strsplit((st[which(st$hy_name_s == "CMxCP"),]$comp_genes),";")[[1]]
opxom_g <-  strsplit((st[which(st$hy_name_s == "OPxOM"),]$comp_genes),";")[[1]]
naxca_g <-  strsplit((st[which(st$hy_name_s == "NAxCA"),]$comp_genes),";")[[1]]
upxca_g <-  strsplit((st[which(st$hy_name_s == "UPxCA"),]$comp_genes),";")[[1]]
oaxca_g <-  strsplit((st[which(st$hy_name_s == "OAxCA"),]$comp_genes),";")[[1]]
caxoa_g <-  strsplit((st[which(st$hy_name_s == "CAxOA"),]$comp_genes),";")[[1]]
cmxom_g <-  strsplit((st[which(st$hy_name_s == "CMxOM"),]$comp_genes),";")[[1]]
opxcp_g <-  strsplit((st[which(st$hy_name_s == "OPxCP"),]$comp_genes),";")[[1]]


venn(list(compensatory=caxcm_g,misexpressed=caxcm), ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
venn(list(compensatory=cmxca_g,misexpressed=cmxca), ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
venn(list(compensatory=oaxom_g,misexpressed=oaxom), ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
venn(list(compensatory=omxoa_g,misexpressed=omxoa), ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
venn(list(compensatory=caxcp_g,misexpressed=caxcp), ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
venn(list(compensatory=oaxop_g,misexpressed=oaxop), ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
venn(list(compensatory=opxoa_g,misexpressed=opxoa), ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
venn(list(compensatory=cmxcp_g,misexpressed=cmxcp), ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
venn(list(compensatory=opxom_g,misexpressed=opxom), ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
venn(list(compensatory=naxca_g,misexpressed=naxca), ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
venn(list(compensatory=upxca_g,misexpressed=upxca), ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
venn(list(compensatory=oaxca_g,misexpressed=oaxca), ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
venn(list(compensatory=caxoa_g,misexpressed=caxoa), ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
venn(list(compensatory=cmxom_g,misexpressed=cmxom), ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
venn(list(compensatory=opxcp_g,misexpressed=opxcp), ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)


cis         <- c()  
comp        <- c()  
trans       <- c()  
mis_ase_comp<- c()  
mis_ase_cis <- c()  

for (row_num in c(1:nrow(st)))
{
  cis                      <- c(cis, strsplit(st$cis_genes[row_num],";")[[1]])
  comp                    <- c(comp, strsplit(st$comp_genes[row_num],";")[[1]])
  trans                  <- c(trans, strsplit(st$trans_genes[row_num],";")[[1]])
  mis_ase_comp    <- c(mis_ase_comp, strsplit(st$mis_ase_comp_genes[row_num],";")[[1]])
  mis_ase_cis     <- c(mis_ase_cis,  strsplit(st$mis_ase_cis_genes[row_num],";")[[1]])
}




#####
######################################################
######################################################
##### DE/ME overlap ##################################
######################################################
library(venn)
DE_genes_dir <- "C:/Users/jmcgirr/Documents/all_2018_samples/conditions/"
venn_dir <- "C:/Users/jmcgirr/Documents/all_2018_samples/mse_venns/"
stage <- "8dpf"
# load misexpressed
{
# am crosses
caxcm <- read.csv(paste(DE_genes_dir,"DE_CRPA_and_CRPM_vs_CAxCM_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
cmxca <- read.csv(paste(DE_genes_dir,"DE_CRPM_and_CRPA_vs_CMxCA_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
oaxom <- read.csv(paste(DE_genes_dir,"DE_OSPA_and_OSPM_vs_OAxOM_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
omxoa <- read.csv(paste(DE_genes_dir,"DE_OSPM_and_OSPA_vs_OMxOA_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
# ap crosses
caxcp <- read.csv(paste(DE_genes_dir,"DE_CRPA_and_CRPP_vs_CAxCP_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
oaxop <- read.csv(paste(DE_genes_dir,"DE_OSPA_and_OSPP_vs_OAxOP_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
opxoa <- read.csv(paste(DE_genes_dir,"DE_OSPP_and_OSPA_vs_OPxOA_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
# specialist crosses
cmxcp <- read.csv(paste(DE_genes_dir,"DE_CRPM_and_CRPP_vs_CMxCP_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
opxom <- read.csv(paste(DE_genes_dir,"DE_OSPM_and_OSPP_vs_OPxOM_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
# outgroup crosses
naxca <-  read.csv(paste(DE_genes_dir,"DE_CRPA_and_NCA_vs_NAxCA_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
upxca <-  read.csv(paste(DE_genes_dir,"DE_CRPA_and_UPxUA_vs_UPxCA_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
# lake crosses
oaxca <- read.csv(paste(DE_genes_dir,"DE_CRPA_and_OSPA_vs_OAxCA_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
caxoa <- read.csv(paste(DE_genes_dir,"DE_CRPA_and_OSPA_vs_CAxOA_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
cmxom <- read.csv(paste(DE_genes_dir,"DE_CRPM_and_OSPM_vs_CMxOM_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
opxcp <- read.csv(paste(DE_genes_dir,"DE_CRPP_and_OSPP_vs_OPxCP_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)

# sig misexpressed
caxcm <- caxcm[which(caxcm$padj <= 0.05),] 
cmxca <- cmxca[which(cmxca$padj <= 0.05),] 
oaxom <- oaxom[which(oaxom$padj <= 0.05),] 
omxoa <- omxoa[which(omxoa$padj <= 0.05),] 
caxcp <- caxcp[which(caxcp$padj <= 0.05),] 
oaxop <- oaxop[which(oaxop$padj <= 0.05),] 
opxoa <- opxoa[which(opxoa$padj <= 0.05),] 
cmxcp <- cmxcp[which(cmxcp$padj <= 0.05),] 
opxom <- opxom[which(opxom$padj <= 0.05),] 
naxca <- naxca[which(naxca$padj <= 0.05),] 
upxca <- upxca[which(upxca$padj <= 0.05),] 
oaxca <- oaxca[which(oaxca$padj <= 0.05),] 
caxoa <- caxoa[which(caxoa$padj <= 0.05),] 
cmxom <- cmxom[which(cmxom$padj <= 0.05),] 
opxcp <- opxcp[which(opxcp$padj <= 0.05),] 

caxcm_8_ME <- caxcm$related_accession
cmxca_8_ME <- cmxca$related_accession
oaxom_8_ME <- oaxom$related_accession
omxoa_8_ME <- omxoa$related_accession
caxcp_8_ME <- caxcp$related_accession
oaxop_8_ME <- oaxop$related_accession
opxoa_8_ME <- opxoa$related_accession
cmxcp_8_ME <- cmxcp$related_accession
opxom_8_ME <- opxom$related_accession
naxca_8_ME <- naxca$related_accession
upxca_8_ME <- upxca$related_accession
oaxca_8_ME <- oaxca$related_accession
caxoa_8_ME <- caxoa$related_accession
cmxom_8_ME <- cmxom$related_accession
opxcp_8_ME <- opxcp$related_accession
}
# load DE
{
  # am crosses
  caxcm <- read.csv(paste(DE_genes_dir,"DE_CRPA_vs_CRPM_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  oaxom <- read.csv(paste(DE_genes_dir,"DE_OSPA_vs_OSPM_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  # ap crosses
  caxcp <- read.csv(paste(DE_genes_dir,"DE_CRPA_vs_CRPP_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  oaxop <- read.csv(paste(DE_genes_dir,"DE_OSPA_vs_OSPP_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  # specialist crosses
  cmxcp <- read.csv(paste(DE_genes_dir,"DE_CRPM_vs_CRPP_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  omxop <- read.csv(paste(DE_genes_dir,"DE_OSPM_vs_OSPP_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  # outgroup crosses
  caxna <-  read.csv(paste(DE_genes_dir,"DE_CRPA_vs_NCA_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  oaxna <-  read.csv(paste(DE_genes_dir,"DE_OSPA_vs_NCA_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  caxup <-  read.csv(paste(DE_genes_dir,"DE_CRPA_vs_UPxUA_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  oaxup <-  read.csv(paste(DE_genes_dir,"DE_OSPA_vs_UPxUA_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  
  # lake crosses
  caxoa <- read.csv(paste(DE_genes_dir,"DE_CRPA_vs_OSPA_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  cmxom <- read.csv(paste(DE_genes_dir,"DE_CRPM_vs_OSPM_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  cpxop <- read.csv(paste(DE_genes_dir,"DE_CRPP_vs_OSPP_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  
  # number of informitive transcripts
  caxcm_n_8_DE <- nrow(caxcm)
  oaxom_n_8_DE <- nrow(oaxom)
  caxcp_n_8_DE <- nrow(caxcp)
  oaxop_n_8_DE <- nrow(oaxop)
  cmxcp_n_8_DE <- nrow(cmxcp)
  omxop_n_8_DE <- nrow(omxop)
  caxna_n_8_DE <- nrow(caxna)
  oaxna_n_8_DE <- nrow(oaxna)
  caxup_n_8_DE <- nrow(caxup)
  oaxup_n_8_DE <- nrow(oaxup)
  caxoa_n_8_DE <- nrow(caxoa)
  cmxom_n_8_DE <- nrow(cmxom)
  cpxop_n_8_DE <- nrow(cpxop)
  
  # sig DE
  caxcm <- caxcm[which(caxcm$padj <= 0.05),] 
  oaxom <- oaxom[which(oaxom$padj <= 0.05),]
  caxcp <- caxcp[which(caxcp$padj <= 0.05),]
  oaxop <- oaxop[which(oaxop$padj <= 0.05),]
  cmxcp <- cmxcp[which(cmxcp$padj <= 0.05),]
  omxop <- omxop[which(omxop$padj <= 0.05),]
  caxna <- caxna[which(caxna$padj <= 0.05),]
  oaxna <- oaxna[which(oaxna$padj <= 0.05),]
  caxup <- caxup[which(caxup$padj <= 0.05),]
  oaxup <- oaxup[which(oaxup$padj <= 0.05),]
  caxoa <- caxoa[which(caxoa$padj <= 0.05),]
  cmxom <- cmxom[which(cmxom$padj <= 0.05),]
  cpxop <- cpxop[which(cpxop$padj <= 0.05),]
  
  caxcm_8_DE <- caxcm$related_accession
  oaxom_8_DE <- oaxom$related_accession
  caxcp_8_DE <- caxcp$related_accession
  oaxop_8_DE <- oaxop$related_accession
  cmxcp_8_DE <- cmxcp$related_accession
  omxop_8_DE <- omxop$related_accession
  caxna_8_DE <- caxna$related_accession
  oaxna_8_DE <- oaxna$related_accession
  caxup_8_DE <- caxup$related_accession
  oaxup_8_DE <- oaxup$related_accession
  caxoa_8_DE <- caxoa$related_accession
  cmxom_8_DE <- cmxom$related_accession
  cpxop_8_DE <- cpxop$related_accession
}
stage <- "48hpf"
# load misexpressed
{
  # am crosses
  caxcm <- read.csv(paste(DE_genes_dir,"DE_CRPA_and_CRPM_vs_CAxCM_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  cmxca <- read.csv(paste(DE_genes_dir,"DE_CRPM_and_CRPA_vs_CMxCA_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  oaxom <- read.csv(paste(DE_genes_dir,"DE_OSPA_and_OSPM_vs_OAxOM_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  omxoa <- read.csv(paste(DE_genes_dir,"DE_OSPM_and_OSPA_vs_OMxOA_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  # ap crosses
  caxcp <- read.csv(paste(DE_genes_dir,"DE_CRPA_and_CRPP_vs_CAxCP_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  oaxop <- read.csv(paste(DE_genes_dir,"DE_OSPA_and_OSPP_vs_OAxOP_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  opxoa <- read.csv(paste(DE_genes_dir,"DE_OSPP_and_OSPA_vs_OPxOA_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  # specialist crosses
  cmxcp <- read.csv(paste(DE_genes_dir,"DE_CRPM_and_CRPP_vs_CMxCP_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  opxom <- read.csv(paste(DE_genes_dir,"DE_OSPM_and_OSPP_vs_OPxOM_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  # outgroup crosses
  naxca <-  read.csv(paste(DE_genes_dir,"DE_CRPA_and_NCA_vs_NAxCA_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  upxca <-  read.csv(paste(DE_genes_dir,"DE_CRPA_and_UPxUA_vs_UPxCA_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  # lake crosses
  oaxca <- read.csv(paste(DE_genes_dir,"DE_CRPA_and_OSPA_vs_OAxCA_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  caxoa <- read.csv(paste(DE_genes_dir,"DE_CRPA_and_OSPA_vs_CAxOA_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  cmxom <- read.csv(paste(DE_genes_dir,"DE_CRPM_and_OSPM_vs_CMxOM_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  opxcp <- read.csv(paste(DE_genes_dir,"DE_CRPP_and_OSPP_vs_OPxCP_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  
  # sig misexpressed
  caxcm <- caxcm[which(caxcm$padj <= 0.05),] 
  cmxca <- cmxca[which(cmxca$padj <= 0.05),] 
  oaxom <- oaxom[which(oaxom$padj <= 0.05),] 
  omxoa <- omxoa[which(omxoa$padj <= 0.05),] 
  caxcp <- caxcp[which(caxcp$padj <= 0.05),] 
  oaxop <- oaxop[which(oaxop$padj <= 0.05),] 
  opxoa <- opxoa[which(opxoa$padj <= 0.05),] 
  cmxcp <- cmxcp[which(cmxcp$padj <= 0.05),] 
  opxom <- opxom[which(opxom$padj <= 0.05),] 
  naxca <- naxca[which(naxca$padj <= 0.05),] 
  upxca <- upxca[which(upxca$padj <= 0.05),] 
  oaxca <- oaxca[which(oaxca$padj <= 0.05),] 
  caxoa <- caxoa[which(caxoa$padj <= 0.05),] 
  cmxom <- cmxom[which(cmxom$padj <= 0.05),] 
  opxcp <- opxcp[which(opxcp$padj <= 0.05),] 
  
  caxcm_48_ME <- caxcm$related_accession
  cmxca_48_ME <- cmxca$related_accession
  oaxom_48_ME <- oaxom$related_accession
  omxoa_48_ME <- omxoa$related_accession
  caxcp_48_ME <- caxcp$related_accession
  oaxop_48_ME <- oaxop$related_accession
  opxoa_48_ME <- opxoa$related_accession
  cmxcp_48_ME <- cmxcp$related_accession
  opxom_48_ME <- opxom$related_accession
  naxca_48_ME <- naxca$related_accession
  upxca_48_ME <- upxca$related_accession
  oaxca_48_ME <- oaxca$related_accession
  caxoa_48_ME <- caxoa$related_accession
  cmxom_48_ME <- cmxom$related_accession
  opxcp_48_ME <- opxcp$related_accession
}
# load DE
{
  # am crosses
  caxcm <- read.csv(paste(DE_genes_dir,"DE_CRPA_vs_CRPM_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  oaxom <- read.csv(paste(DE_genes_dir,"DE_OSPA_vs_OSPM_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  # ap crosses
  caxcp <- read.csv(paste(DE_genes_dir,"DE_CRPA_vs_CRPP_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  oaxop <- read.csv(paste(DE_genes_dir,"DE_OSPA_vs_OSPP_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  # specialist crosses
  cmxcp <- read.csv(paste(DE_genes_dir,"DE_CRPM_vs_CRPP_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  omxop <- read.csv(paste(DE_genes_dir,"DE_OSPM_vs_OSPP_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  # outgroup crosses
  caxna <-  read.csv(paste(DE_genes_dir,"DE_CRPA_vs_NCA_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  oaxna <-  read.csv(paste(DE_genes_dir,"DE_OSPA_vs_NCA_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  caxup <-  read.csv(paste(DE_genes_dir,"DE_CRPA_vs_UPxUA_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  oaxup <-  read.csv(paste(DE_genes_dir,"DE_OSPA_vs_UPxUA_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  
  # lake crosses
  caxoa <- read.csv(paste(DE_genes_dir,"DE_CRPA_vs_OSPA_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  cmxom <- read.csv(paste(DE_genes_dir,"DE_CRPM_vs_OSPM_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  cpxop <- read.csv(paste(DE_genes_dir,"DE_CRPP_vs_OSPP_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  
  # number of informitive transcripts
  caxcm_n_48_DE <- nrow(caxcm)
  oaxom_n_48_DE <- nrow(oaxom)
  caxcp_n_48_DE <- nrow(caxcp)
  oaxop_n_48_DE <- nrow(oaxop)
  cmxcp_n_48_DE <- nrow(cmxcp)
  omxop_n_48_DE <- nrow(omxop)
  caxna_n_48_DE <- nrow(caxna)
  oaxna_n_48_DE <- nrow(oaxna)
  caxup_n_48_DE <- nrow(caxup)
  oaxup_n_48_DE <- nrow(oaxup)
  caxoa_n_48_DE <- nrow(caxoa)
  cmxom_n_48_DE <- nrow(cmxom)
  cpxop_n_48_DE <- nrow(cpxop)
  
  # sig DE
  caxcm <- caxcm[which(caxcm$padj <= 0.05),] 
  oaxom <- oaxom[which(oaxom$padj <= 0.05),]
  caxcp <- caxcp[which(caxcp$padj <= 0.05),]
  oaxop <- oaxop[which(oaxop$padj <= 0.05),]
  cmxcp <- cmxcp[which(cmxcp$padj <= 0.05),]
  omxop <- omxop[which(omxop$padj <= 0.05),]
  caxna <- caxna[which(caxna$padj <= 0.05),]
  oaxna <- oaxna[which(oaxna$padj <= 0.05),]
  caxup <- caxup[which(caxup$padj <= 0.05),]
  oaxup <- oaxup[which(oaxup$padj <= 0.05),]
  caxoa <- caxoa[which(caxoa$padj <= 0.05),]
  cmxom <- cmxom[which(cmxom$padj <= 0.05),]
  cpxop <- cpxop[which(cpxop$padj <= 0.05),]
  
  caxcm_48_DE <- caxcm$related_accession
  oaxom_48_DE <- oaxom$related_accession
  caxcp_48_DE <- caxcp$related_accession
  oaxop_48_DE <- oaxop$related_accession
  cmxcp_48_DE <- cmxcp$related_accession
  omxop_48_DE <- omxop$related_accession
  caxna_48_DE <- caxna$related_accession
  oaxna_48_DE <- oaxna$related_accession
  caxup_48_DE <- caxup$related_accession
  oaxup_48_DE <- oaxup$related_accession
  caxoa_48_DE <- caxoa$related_accession
  cmxom_48_DE <- cmxom$related_accession
  cpxop_48_DE <- cpxop$related_accession
}

venn(list(cmxcp_48_DE=cmxcp_48_DE,cmxcp_48_ME=cmxcp_48_ME), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
venn(list(omxop_48_DE=omxop_48_DE,opxom_48_ME=opxom_48_ME), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
venn(list(cmxcp_8_DE=cmxcp_8_DE ,cmxcp_8_ME=cmxcp_8_ME), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
venn(list(omxop_8_DE=omxop_8_DE,opxom_8_ME=opxom_8_ME), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)

#other venn overlaps
{
venn(list(cmxcp_48_DE=cmxcp_48_DE,cmxcp_48_ME=cmxcp_48_ME,
          cmxcp_8_DE=cmxcp_8_DE ,cmxcp_8_ME=cmxcp_8_ME), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
venn(list(omxop_48_DE=omxop_48_DE,opxom_48_ME=opxom_48_ME,
          omxop_8_DE=omxop_8_DE,opxom_8_ME=opxom_8_ME), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)

venn(list(cmxcp_48_DE=cmxcp_48_DE,cmxcp_48_ME=cmxcp_48_ME,
          omxop_48_DE=omxop_48_DE,opxom_48_ME=opxom_48_ME), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
venn(list(cmxcp_8_DE=cmxcp_8_DE,cmxcp_8_ME=cmxcp_8_ME,
          omxop_8_DE=omxop_8_DE,opxom_8_ME=opxom_8_ME), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
}

# parallel incompatibilities
venn(list(caxcm_48_ME=caxcm_48_ME,cmxca_48_ME=cmxca_48_ME,caxcp_48_ME=caxcp_48_ME), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
intersect(caxcp_48_ME,intersect(caxcm_48_ME,cmxca_48_ME))
venn(list(caxcm_8_ME=caxcm_8_ME,cmxca_8_ME=cmxca_8_ME,caxcp_8_ME=caxcp_8_ME), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
intersect(caxcp_8_ME,intersect(caxcm_8_ME,cmxca_8_ME))

venn(list(oaxom_48_ME=oaxom_48_ME,omxoa_48_ME=omxoa_48_ME,oaxop_48_ME=oaxop_48_ME,opxoa_48_ME=opxoa_48_ME), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
venn(list(oaxom_8_ME=oaxom_8_ME,omxoa_8_ME=omxoa_8_ME,oaxop_8_ME=oaxop_8_ME,opxoa_8_ME=opxoa_8_ME), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
#intersect(caxcp_48_ME,intersect(caxcm_48_ME,cmxca_48_ME))
venn(list(caxcm_8_ME=caxcm_8_ME,cmxca_8_ME=cmxca_8_ME,caxcp_8_ME=caxcp_8_ME), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
#intersect(caxcp_8_ME,intersect(caxcm_8_ME,cmxca_8_ME))

# parallel DE http://nemates.org/MA/progs/overlap_stats.html
# crescent pond
venn(list(caxcm_48_DE=caxcm_48_DE,caxcp_48_DE=caxcp_48_DE), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
#intersect(caxcp_48_DE,caxcm_48_DE)
test <- matrix(c(length(intersect(caxcp_48_DE,caxcm_48_DE)),
                 length(setdiff(caxcp_48_DE,caxcm_48_DE)),
                 length(setdiff(caxcm_48_DE,caxcp_48_DE)),
                 (caxcp_n_48_DE - length(unique(c(caxcm_48_DE,caxcp_48_DE))))),ncol=2)
fisher.test(test, alternative = 'g')
venn(list(caxcm_8_DE=caxcm_8_DE,caxcp_8_DE=caxcp_8_DE), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
#intersect(caxcp_8_DE,intersect(caxcm_8_DE,cmxca_8_DE))
#intersect(intersect(caxcm_48_DE,caxcp_48_DE),intersect(caxcm_8_DE,caxcp_8_DE))
test <- matrix(c(length(intersect(caxcp_8_DE,caxcm_8_DE)),
                 length(setdiff(caxcp_8_DE,caxcm_8_DE)),
                 length(setdiff(caxcm_8_DE,caxcp_8_DE)),
                 (caxcp_n_8_DE - length(unique(c(caxcm_8_DE,caxcp_8_DE))))),ncol=2)
fisher.test(test, alternative = 'g')

# osprey
venn(list(oaxom_48_DE=oaxom_48_DE,oaxop_48_DE=oaxop_48_DE), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
#intersect(oaxom_48_DE,oaxop_48_DE)
test <- matrix(c(length(intersect(oaxop_48_DE,oaxom_48_DE)),
                 length(setdiff(oaxop_48_DE,oaxom_48_DE)),
                 length(setdiff(oaxom_48_DE,oaxop_48_DE)),
                 (oaxop_n_48_DE - length(unique(c(oaxom_48_DE,oaxop_48_DE))))),ncol=2)
fisher.test(test, alternative = 'g')
venn(list(oaxom_8_DE=oaxom_8_DE,oaxop_8_DE=oaxop_8_DE), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
#intersect(oaxom_8_DE,oaxop_8_DE)
#intersect(intersect(oaxom_48_DE,oaxop_48_DE),intersect(oaxom_8_DE,oaxop_8_DE))
test <- matrix(c(length(intersect(oaxop_8_DE,oaxom_8_DE)),
                 length(setdiff(oaxop_8_DE,oaxom_8_DE)),
                 length(setdiff(oaxom_8_DE,oaxop_8_DE)),
                 (oaxop_n_8_DE - length(unique(c(oaxom_8_DE,oaxop_8_DE))))),ncol=2)
fisher.test(test, alternative = 'g')


adaptive_incompatibilities <- c(intersect(cmxcp_48_DE,cmxcp_48_ME),intersect(omxop_48_DE,opxom_48_ME),
                                intersect(cmxcp_8_DE,cmxcp_8_ME),intersect(omxop_8_DE,opxom_8_ME))
intersect(adaptive_incompatibilities,cis)

length(adaptive_incompatibilities)
intersect(intersect(cmxcp_48_DE,cmxcp_48_ME),cis)
intersect(intersect(omxop_48_DE,opxom_48_ME),cis)
intersect(intersect(cmxcp_8_DE,cmxcp_8_ME),cis)
intersect(intersect(omxop_8_DE,opxom_8_ME),cis)



xm_to_gene <- function(xms){

xm_transcripts <- data.frame(related_accession = xms)
xm_transcripts <- merge(xm_transcripts,final_features, by = c("related_accession"))
xm_transcripts<- merge(xm_transcripts, blast_key, by = c("product_accession"))
print("danio paralogs")
print(xm_transcripts$zeb_gene_symbol_one_way)
print("danio paralogs annotated for skeletal effects")
xm_cranial <- merge(xm_transcripts, cranial_genes, by = c("zeb_gene_symbol_one_way"))
print(xm_cranial$zeb_gene_symbol_one_way)
print("cyprinodon annotated for skeletal effects")
xm_cranial_cyp <- merge(xm_transcripts, cranial_genes, by = c("symbol"))
print(xm_cranial_cyp$symbol)
}

xm_to_gene(adaptive_incompatibilities)
xm_to_gene(intersect(cmxcp_48_DE,cmxcp_48_ME))
xm_to_gene(intersect(omxop_48_DE,omxop_48_ME))
xm_to_gene(intersect(cmxcp_8_DE,cmxcp_8_ME))
xm_to_gene(intersect(omxop_8_DE,omxop_8_ME))
xms <- adaptive_incompatibilities
xms <- test

test <- c("XM_015404656.1", "XM_015394487.1", "XM_015395742.1")
test <- c("XM_015404656.1","XM_015389364.1","XM_015402498.1","XM_015398548.1","XM_015400588.1","XM_015402320.1","XM_015405217.1","XM_015381015.1","XM_015394487.1","XM_015387235.1","XM_015399782.1","XM_015394033.1","XM_015387156.1","XM_015400083.1","XM_015391078.1","XM_015391944.1","XM_015379642.1","XM_015396313.1","XM_015394620.1","XM_015388212.1","XM_015390249.1","XM_015384304.1","XM_015395742.1","XM_015381604.1","XM_015373232.1","XM_015371807.1","XM_015377296.1","XM_015397675.1","XM_015395749.1")
xm_to_gene(test)

#####
######################################################
######################################################
##### ASE ############################################
######################################################


library(reshape2)
inheritance_comps <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/ase_data/ase_comparisons.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
DE_genes_dir <- "C:/Users/jmcgirr/Documents/all_2018_samples/conditions/"
ase_plots_dir <- "D:/Martin Lab/rna_2018/ASE/plots/"
geneiase_dir <- "C:/Users/jmcgirr/Documents/all_2018_samples/ase_data/geneiase/"
mbased_dir <- "C:/Users/jmcgirr/Documents/all_2018_samples/ase_data/mbased/"
maternal_counts_dir <- "C:/Users/jmcgirr/Documents/all_2018_samples/ase_data/maternal_counts/"
master <- read.table("D:/Martin Lab/rna_2018/all_2018_samples/table_maker_master_outlier_rm.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
norm_cts <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/norm_cts_all_2018_no_seq_round_control.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
all_snps <- read.csv("C:/Users/jmcgirr/Documents/all_2018_samples/ase_data/rna_2018_snps.txt",na.strings=c("","NA"), header = TRUE, stringsAsFactors = FALSE, sep = "\t")
snp_features <- read.csv("C:/Users/jmcgirr/Documents/all_2018_samples/ase_data/rna_2018_features.txt",na.strings=c("","NA"), header = TRUE, stringsAsFactors = FALSE, sep = "\t")
#all_snps <- all_snps[, -grep("HP", colnames(all_snps))]
#write.table(snp_features,"C:/Users/jmcgirr/Documents/all_2018_samples/ase_data/rna_2018_features.txt", row.names = FALSE, quote= FALSE,sep="\t")
#i <- 27
#subset_comps <- c(1:3)

parent_pops1_name_s<- c()
parent_pops2_name_s<- c()
hy_name_s<- c()
stage_s<- c()
st_informative_genes_s<- c()
st_cis_s<- c()
st_trans_s<- c()
st_cis_plus_trans_s<- c()
st_cis_by_trans_s<- c()
st_comp_s<- c()
st_mis_s<- c()
st_mis_ase_s<- c()
st_mis_ase_comp_s<- c()
st_mis_ase_cis_s<- c()
st_cis_genes<- c()
st_trans_genes<- c()
st_cis_plus_trans_genes<- c()
st_cis_by_trans_genes<- c()
st_comp_genes<- c()
st_mis_ase_comp_genes<- c()
st_mis_ase_cis_genes<- c()
informative_genes_s<- c()
cis_s<- c()
trans_s<- c()
cis_plus_trans_s<- c()
cis_by_trans_s<- c()
comp_s<- c()
mis_s<- c()
mis_ase_s<- c()
mis_ase_comp_s<- c()
mis_ase_cis_s<- c()
cis_genes<- c()
trans_genes<- c()
cis_plus_trans_genes<- c()
cis_by_trans_genes<- c()
comp_genes<- c()
mis_ase_comp_genes<- c()
mis_ase_cis_genes<- c()
sl_informative_genes_s<- c()
sl_cis_s<- c()
sl_trans_s<- c()
sl_comp_s<- c()
sl_mis_s<- c()
sl_mis_ase_s<- c()
sl_mis_ase_comp_s<- c()
sl_mis_ase_cis_s<- c()
sl_cis_genes<- c()
sl_trans_genes<- c()
sl_comp_genes<- c()
sl_mis_ase_comp_genes<- c()
sl_mis_ase_cis_genes<- c()

st_mis_trans_genes<- c()
mis_trans_genes<- c()
sl_mis_trans_genes<- c()

#subset_comps <- c(1,19)

#for (i in subset_comps)
for (i in (1:nrow(inheritance_comps)))
{
  parent_pops1_name <- inheritance_comps$parent_pops1[i]
  parent_pops2_name <- inheritance_comps$parent_pops2[i]
  hy_name <- inheritance_comps$hybrids[i]
  stage <- inheritance_comps$stage[i]
  parent_pops_v_hy  <- paste(DE_genes_dir, "DE_", parent_pops1_name,"_and_",parent_pops2_name, "_vs_",hy_name, "_",stage,"_genes",".csv",sep = "")
  parent_pops1_v_hy <- paste(DE_genes_dir, "DE_", parent_pops1_name, "_vs_",hy_name, "_",stage,"_genes",".csv",sep = "")
  parent_pops2_v_hy <- paste(DE_genes_dir, "DE_", parent_pops2_name, "_vs_",hy_name, "_",stage,"_genes",".csv",sep = "")
  parent_pops1_v_parent_pops2 <- paste(DE_genes_dir, "DE_", parent_pops1_name, "_vs_",parent_pops2_name, "_",stage,"_genes",".csv",sep = "")
  outfile_plot  <- paste(ase_plots_dir, hy_name, "_",stage,"_misexpression",".tiff",sep = "")
  plot_title  <- paste(hy_name,stage,"misexpression",sep = " ")
  
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
  
  #Hybrid inheritance was considered additive if gene expression was intermediate between
  #generalists and molluscivores with significant differential expression between generalists and
  #molluscivores. Inheritance was dominant if expression was intermediate between parental
  #species and hybrid expression was significantly different from one parent but not the other.
  #Genes showing misexpression in hybrids showed transgressive inheritance, where hybrid gene
  #expression was significantly higher (overdominant) or lower (underdominant) than parental
  #populations.
  
  p1_inds <- master[which(master$f1 == parent_pops1_name & master$stage == stage),]
  p2_inds <- master[which(master$f1 == parent_pops2_name & master$stage == stage),]
  hy_names <- strsplit(hy_name,"_&_")[[1]]
  hy_inds <- master[which(master$f1 == hy_names & master$stage == stage),]
  p1_inds <- p1_inds$sample
  p2_inds <- p2_inds$sample
  hy_inds <- hy_inds$sample
  #
  
  # find alternate homozygous snps in p1 and p2 (100% genotyping rate)
  # find genes showing ase in p1 and p2
  
  # no ase in all parents
  
  p1Ase <- data.frame(related_accession=character(), 
                      mrnaID=character(),
                      stringsAsFactors=FALSE)
  keeps <- c("CHROM","POS")
  p1_snps <- all_snps[keeps]
  
  ind <- p1_inds[[1]]
  #gAse <- read.table(paste(geneiase_dir,ind, "_geneiase_ase.txt.static.gene.pval.tab",sep = ""), header = TRUE, stringsAsFactors = FALSE) 
  mAse <- read.table(paste(mbased_dir,ind, "_mbased_ase_unphased.txt",sep = ""), header = TRUE, stringsAsFactors = FALSE) 
  
  #gAse <- gAse[which(gAse$p.nom <= 0.05),]
  mAse <- mAse[which(mAse$pValueASE <= 0.05),]
  #gAse <- cbind(gAse, colsplit(gAse$feat, ";", c("related_accession", "gene_name")))
  mAse <- cbind(mAse, colsplit(mAse$mrnaID, ";", c("related_accession", "gene_name"))) 
  #gAse$mrnaID <- gAse$feat
  #ase <- merge(gAse,mAse, by = c("related_accession"))
  ase <- mAse
  keeps <- c("related_accession", "mrnaID")
  p1Ase <- ase[keeps]
  
  for (ind in p1_inds)
  {
    #gAse <- read.table(paste(geneiase_dir,ind, "_geneiase_ase.txt.static.gene.pval.tab",sep = ""), header = TRUE, stringsAsFactors = FALSE) 
    mAse <- read.table(paste(mbased_dir,ind, "_mbased_ase_unphased.txt",sep = ""), header = TRUE, stringsAsFactors = FALSE) 
    
    #gAse <- gAse[which(gAse$p.nom <= 0.05),]
    mAse <- mAse[which(mAse$pValueASE <= 0.05),]
    #gAse <- cbind(gAse, colsplit(gAse$feat, ";", c("related_accession", "gene_name")))
    mAse <- cbind(mAse, colsplit(mAse$mrnaID, ";", c("related_accession", "gene_name"))) 
    #ase <- merge(gAse,mAse, by = c("related_accession"))
    #gAse$mrnaID <- gAse$feat
    ase <- mAse
    keeps <- c("related_accession", "mrnaID")
    ase <- ase[keeps]
    
    p1Ase <- merge(p1Ase, ase, by = c("related_accession", "mrnaID"))
    keeps <- c("CHROM","POS", paste(ind, ".GT", sep = ""))
    snps <- all_snps[keeps]
    snps <- cbind(snps, colsplit(snps[,3], "/", c(paste(ind,"_a1", sep = ""), paste(ind,"_a2", sep = ""))))
    snps <- snps[, -3]
    p1_snps <- merge(p1_snps, snps, by = c("CHROM", "POS"))
    
  }
  
  p2Ase <- data.frame(related_accession=character(), 
                      mrnaID=character(),
                      stringsAsFactors=FALSE)
  keeps <- c("CHROM","POS")
  p2_snps <- all_snps[keeps]
  
  ind <- p2_inds[[1]]
  #gAse <- read.table(paste(geneiase_dir,ind, "_geneiase_ase.txt.static.gene.pval.tab",sep = ""), header = TRUE, stringsAsFactors = FALSE) 
  mAse <- read.table(paste(mbased_dir,ind, "_mbased_ase_unphased.txt",sep = ""), header = TRUE, stringsAsFactors = FALSE) 
  
  #gAse <- gAse[which(gAse$p.nom <= 0.05),]
  mAse <- mAse[which(mAse$pValueASE <= 0.05),]
  #gAse <- cbind(gAse, colsplit(gAse$feat, ";", c("related_accession", "gene_name")))
  mAse <- cbind(mAse, colsplit(mAse$mrnaID, ";", c("related_accession", "gene_name"))) 
  #gAse$mrnaID <- gAse$feat
  #ase <- merge(gAse,mAse, by = c("related_accession"))
  ase <- mAse
  keeps <- c("related_accession", "mrnaID")
  p2Ase <- ase[keeps]
  
  for (ind in p2_inds)
  {
    #gAse <- read.table(paste(geneiase_dir,ind, "_geneiase_ase.txt.static.gene.pval.tab",sep = ""), header = TRUE, stringsAsFactors = FALSE) 
    mAse <- read.table(paste(mbased_dir,ind, "_mbased_ase_unphased.txt",sep = ""), header = TRUE, stringsAsFactors = FALSE) 
    
    #gAse <- gAse[which(gAse$p.nom <= 0.05),]
    mAse <- mAse[which(mAse$pValueASE <= 0.05),]
    #gAse <- cbind(gAse, colsplit(gAse$feat, ";", c("related_accession", "gene_name")))
    mAse <- cbind(mAse, colsplit(mAse$mrnaID, ";", c("related_accession", "gene_name"))) 
    #ase <- merge(gAse,mAse, by = c("related_accession"))
    #gAse$mrnaID <- gAse$feat
    ase <- mAse
    keeps <- c("related_accession", "mrnaID")
    ase <- ase[keeps]
    
    p2Ase <- merge(p2Ase, ase, by = c("related_accession", "mrnaID"))
    keeps <- c("CHROM","POS", paste(ind, ".GT", sep = ""))
    snps <- all_snps[keeps]
    snps <- cbind(snps, colsplit(snps[,3], "/", c(paste(ind,"_a1", sep = ""), paste(ind,"_a2", sep = ""))))
    snps <- snps[, -3]
    p2_snps <- merge(p2_snps, snps, by = c("CHROM", "POS"))
    
  }
  
  ase_in_pure_cross_f1 <- merge(p1Ase,p2Ase, by = c("related_accession", "mrnaID"))
  ase_in_pure_cross_f1 <- ase_in_pure_cross_f1$related_accession
  length(ase_in_pure_cross_f1)
  
  # no ase in any parent
  #p1Ase <- c()
  #keeps <- c("CHROM","POS")
  #p1_snps <- all_snps[keeps]
  #for (ind in p1_inds)
  #  #ind <- "CAE1"
  #{
  #  #gAse <- read.table(paste(geneiase_dir,ind, "_geneiase_ase.txt.static.gene.pval.tab",sep = ""), header = TRUE, stringsAsFactors = FALSE) 
  #  mAse <- read.table(paste(mbased_dir,ind, "_mbased_ase_unphased.txt",sep = ""), header = TRUE, stringsAsFactors = FALSE) 
  #  
  #  #gAse <- gAse[which(gAse$p.nom <= 0.05),]
  #  mAse <- mAse[which(mAse$pValueASE <= 0.05),]
  #  #gAse <- cbind(gAse, colsplit(gAse$feat, ";", c("related_accession", "gene_name")))
  #  mAse <- cbind(mAse, colsplit(mAse$mrnaID, ";", c("related_accession", "gene_name"))) 
  #  #gAse$mrnaID <- gAse$feat
  #  #ase <- merge(gAse,mAse, by = ("related_accession"))
  #  ase <- mAse
  #  p1Ase <- c(p1Ase, ase$related_accession)
  #  
  #  keeps <- c("CHROM","POS", paste(ind, ".GT", sep = ""))
  #  snps <- all_snps[keeps]
  #  snps <- cbind(snps, colsplit(snps[,3], "/", c(paste(ind,"_a1", sep = ""), paste(ind,"_a2", sep = ""))))
  #  snps <- snps[, -3]
  #  p1_snps <- merge(p1_snps, snps, by = c("CHROM", "POS"))
  #  head(p1_snps)
  #  
  #}
  #p2Ase <- c()
  #keeps <- c("CHROM","POS")
  #p2_snps <- all_snps[keeps]
  #for (ind in p2_inds)
  #{
  #  #gAse <- read.table(paste(geneiase_dir,ind, "_geneiase_ase.txt.static.gene.pval.tab",sep = ""), header = TRUE, stringsAsFactors = FALSE) 
  #  mAse <- read.table(paste(mbased_dir,ind, "_mbased_ase_unphased.txt",sep = ""), header = TRUE, stringsAsFactors = FALSE) 
  #  
  #  #gAse <- gAse[which(gAse$p.nom <= 0.05),]
  #  mAse <- mAse[which(mAse$pValueASE <= 0.05),]
  #  #gAse <- cbind(gAse, colsplit(gAse$feat, ";", c("related_accession", "gene_name")))
  #  mAse <- cbind(mAse, colsplit(mAse$mrnaID, ";", c("related_accession", "gene_name"))) 
  #  #gAse$mrnaID <- gAse$feat
  #  #ase <- merge(gAse,mAse, by = ("related_accession"))
  #  ase <- mAse
  #  
  #  p2Ase <- c(p2Ase, ase$related_accession)
  #  
  #  keeps <- c("CHROM","POS", paste(ind, ".GT", sep = ""))
  #  snps <- all_snps[keeps]
  #  snps <- cbind(snps, colsplit(snps[,3], "/", c(paste(ind,"_a1", sep = ""), paste(ind,"_a2", sep = ""))))
  #  snps <- snps[, -3]
  #  p2_snps <- merge(p2_snps, snps, by = c("CHROM", "POS"))
  #}
  #
  #ase_in_pure_cross_f1 <- unique(c(p1Ase,p2Ase))
  #length(ase_in_pure_cross_f1)
  
  for (col_num in c(4:(length(p1_snps))))
  {
    p1_snps <- p1_snps[which(p1_snps[,3] == p1_snps[,col_num]),]
  }
  for (col_num in c(4:(length(p2_snps))))
  {
    p2_snps <- p2_snps[which(p2_snps[,3] == p2_snps[,col_num]),]
  }
  p1_snps <- p1_snps[which(p1_snps[,3] != "."),]
  p2_snps <- p2_snps[which(p2_snps[,3] != "."),]
  p_snps <- merge(p1_snps, p2_snps, by = c("CHROM", "POS"))
  p_snps <- p_snps[which(p_snps[,3] != p_snps[,(length(p1_snps)+1)]),]
  p_snps <- merge(p_snps, snp_features, by = c("CHROM", "POS"))
  keeps <- c("mrnaID")
  alterntae_homozygous_in_parents <- p_snps[keeps]
  keeps <- c("CHROM","POS","mrnaID", "geneID")
  alterntae_homozygous_in_parents_positions <- p_snps[keeps]
  alterntae_homozygous_in_parents_positions <-cbind(alterntae_homozygous_in_parents_positions, colsplit(alterntae_homozygous_in_parents_positions$geneID, ";", c("Geneid", "gene")))
  
  hyAse <- data.frame(related_accession=character(), 
                      mrnaID=character(),
                      stringsAsFactors=FALSE)
  hyAse_any <- data.frame(related_accession=character(), 
                          mrnaID=character(),
                          stringsAsFactors=FALSE)
  ind <- hy_inds[[1]]
  #gAse <- read.table(paste(geneiase_dir,ind, "_geneiase_ase.txt.static.gene.pval.tab",sep = ""), header = TRUE, stringsAsFactors = FALSE) 
  mAse <- read.table(paste(mbased_dir,ind, "_mbased_ase_unphased.txt",sep = ""), header = TRUE, stringsAsFactors = FALSE) 
  
  #gAse <- gAse[which(gAse$p.nom <= 0.05),]
  mAse <- mAse[which(mAse$pValueASE <= 0.05),]
  #gAse <- cbind(gAse, colsplit(gAse$feat, ";", c("related_accession", "gene_name")))
  mAse <- cbind(mAse, colsplit(mAse$mrnaID, ";", c("related_accession", "gene_name"))) 
  #gAse$mrnaID <- gAse$feat
  #ase <- merge(gAse,mAse, by = c("related_accession"))
  ase <- mAse
  keeps <- c("related_accession", "mrnaID")
  hyAse <- ase[keeps]
  
  for (ind in hy_inds)
  {
    #gAse <- read.table(paste(geneiase_dir,ind, "_geneiase_ase.txt.static.gene.pval.tab",sep = ""), header = TRUE, stringsAsFactors = FALSE) 
    mAse <- read.table(paste(mbased_dir,ind, "_mbased_ase_unphased.txt",sep = ""), header = TRUE, stringsAsFactors = FALSE) 
    
    #gAse <- gAse[which(gAse$p.nom <= 0.05),]
    mAse <- mAse[which(mAse$pValueASE <= 0.05),]
    #gAse <- cbind(gAse, colsplit(gAse$feat, ";", c("related_accession", "gene_name")))
    mAse <- cbind(mAse, colsplit(mAse$mrnaID, ";", c("related_accession", "gene_name"))) 
    #ase <- merge(gAse,mAse, by = c("related_accession"))
    #gAse$mrnaID <- gAse$feat
    ase <- mAse
    keeps <- c("related_accession", "mrnaID")
    ase <- ase[keeps]
    
    hyAse_any <- rbind(hyAse_any, ase)
    hyAse <- merge(hyAse, ase, by = c("related_accession", "mrnaID"))
    
  }
  
  hyAse_all_h_no_p <- data.frame(related_accession=character(), mrnaID=character(),hyAse_all_h_no_p=character(),stringsAsFactors=FALSE)
  hyAse_any_h_no_p <- data.frame(related_accession=character(), mrnaID=character(),hyAse_any_h_no_p=character(),stringsAsFactors=FALSE)
  hyAse_any_h_no_p_alt_hom_in_p <- data.frame(related_accession=character(), mrnaID=character(),hyAse_any_h_no_p_alt_hom_in_p=character(),stringsAsFactors=FALSE)
  hyAse_all_h_no_p_alt_hom_in_p <- data.frame(related_accession=character(), mrnaID=character(),hyAse_all_h_no_p_alt_hom_in_p=character(),stringsAsFactors=FALSE)
  
  hyAse_all_h_no_p <- unique(hyAse[! hyAse$related_accession %in% ase_in_pure_cross_f1,])
  hyAse_any_h_no_p <- unique(hyAse_any[! hyAse_any$related_accession %in% ase_in_pure_cross_f1,])
  hyAse_any_h_no_p_alt_hom_in_p <- merge(hyAse_any_h_no_p, alterntae_homozygous_in_parents, by = c("mrnaID"))
  hyAse_all_h_no_p_alt_hom_in_p <- merge(hyAse_all_h_no_p, alterntae_homozygous_in_parents, by = c("mrnaID"))
  
  if(length(hyAse_all_h_no_p$mrnaID) == 0)
  {
    hyAse_all_h_no_p <- data.frame(related_accession=character(), mrnaID=character(),hyAse_all_h_no_p=character(),stringsAsFactors=FALSE)
  }else
  {
    hyAse_all_h_no_p$hyAse_all_h_no_p <- "yes"
  }
  if(length(hyAse_any_h_no_p$mrnaID) == 0)
  {
    hyAse_any_h_no_p <- data.frame(related_accession=character(), mrnaID=character(),hyAse_any_h_no_p=character(),stringsAsFactors=FALSE)
  }else
  {
    hyAse_any_h_no_p$hyAse_any_h_no_p <- "yes"
  }
  if(length(hyAse_any_h_no_p_alt_hom_in_p$mrnaID) == 0)
  {
    hyAse_any_h_no_p_alt_hom_in_p <- data.frame(related_accession=character(), mrnaID=character(),hyAse_any_h_no_p_alt_hom_in_p=character(),stringsAsFactors=FALSE)
  }else
  {
    hyAse_any_h_no_p_alt_hom_in_p$hyAse_any_h_no_p_alt_hom_in_p <- "yes"
  }
  if(length(hyAse_all_h_no_p_alt_hom_in_p$mrnaID) == 0)
  {
    hyAse_all_h_no_p_alt_hom_in_p <- data.frame(related_accession=character(), mrnaID=character(),hyAse_all_h_no_p_alt_hom_in_p=character(),stringsAsFactors=FALSE)
  }else
  {
    hyAse_all_h_no_p_alt_hom_in_p$hyAse_all_h_no_p_alt_hom_in_p <- "yes"
  }
  
  hyAse1 <- merge(unique(hyAse_all_h_no_p),unique(hyAse_any_h_no_p),all = TRUE, by = c("mrnaID", "related_accession"))
  hyAse1 <- merge(hyAse1,unique(hyAse_all_h_no_p_alt_hom_in_p),all = TRUE, by = c("mrnaID", "related_accession"))
  hyAse1 <- merge(hyAse1,unique(hyAse_any_h_no_p_alt_hom_in_p),all = TRUE, by = c("mrnaID", "related_accession"))
  hyAse1$hyAse_all_h_no_p[is.na(hyAse1$hyAse_all_h_no_p)] <- "no"
  hyAse1$hyAse_any_h_no_p[is.na(hyAse1$hyAse_any_h_no_p)] <- "no"
  hyAse1$hyAse_all_h_no_p_alt_hom_in_p[is.na(hyAse1$hyAse_all_h_no_p_alt_hom_in_p)] <- "no"
  hyAse1$hyAse_any_h_no_p_alt_hom_in_p[is.na(hyAse1$hyAse_any_h_no_p_alt_hom_in_p)] <- "no"
  
  inheritance_and_ase <- merge(unique(hyAse1), all_genes, by = c("related_accession"))
  nrow(inheritance_and_ase[which (inheritance_and_ase$hyAse_any_h_no_p == "yes"),])
  nrow(inheritance_and_ase[which (inheritance_and_ase$hyAse_any_h_no_p_alt_hom_in_p == "yes"),])
  nrow(inheritance_and_ase[which (inheritance_and_ase$hyAse_all_h_no_p == "yes"),])
  nrow(inheritance_and_ase[which (inheritance_and_ase$hyAse_all_h_no_p_alt_hom_in_p == "yes"),])
  
  #trans test (lenient - analyze genes with parental counts at het sites)
  
  keeps <- c("Geneid", p1_inds)
  p1_cts <- norm_cts[keeps]
  keeps <- c("Geneid", p2_inds)
  p2_cts <- norm_cts[keeps]
  p1_cts <- data.frame(Geneid=p1_cts[,1], p1_mean_ct=rowMeans(p1_cts[,1:length(p1_inds)+1]))
  p2_cts <- data.frame(Geneid=p2_cts[,1], p2_mean_ct=rowMeans(p2_cts[,1:length(p2_inds)+1]))
  p_cts <- merge(p1_cts,p2_cts, by = c("Geneid"))
  
  head(p2_cts)
  length(unique(alterntae_homozygous_in_parents_positions$Geneid))
  
  hy_cts <- data.frame(CHROM=character(), 
                       POS=numeric(),
                       mrnaID=character(),
                       geneID=character(),
                       stringsAsFactors=FALSE)
  for (ind in hy_inds)
    #ind <- "OVE1"
  {
    hy_mom_name <- paste("hy_mom_cts",ind, sep = "_")
    hy_dad_name <- paste("hy_dad_cts",ind, sep = "_")
    hy <- read.csv(paste(maternal_counts_dir,ind, "_parental_counts_features.txt",sep = ""),na.strings=c("","NA"), header = TRUE, stringsAsFactors = FALSE, sep = "\t")
    hy$hy_mom_cts <- hy$momCount
    hy$hy_dad_cts <- hy$dadCount
    hy$CHROM <- hy$chrom
    hy$POS <- hy$position
    colnames(hy)[colnames(hy)=="hy_mom_cts"] <- hy_mom_name
    colnames(hy)[colnames(hy)=="hy_dad_cts"] <- hy_dad_name
    #keeps <- c("CHROM", "POS","mrnaID","geneID",hy_mom_name,hy_dad_name)
    ## separate hybrids for transTest (instead of using mean counts across unnormalized counts)
    hy_mapa_name <- paste("mapa_cts",ind, sep = "_")
    hy$mapa_cts <- paste(hy$momCount,hy$dadCount, sep = ";")
    colnames(hy)[colnames(hy)=="mapa_cts"] <- hy_mapa_name
    keeps <- c("CHROM", "POS","mrnaID","geneID",hy_mapa_name)
    hy <- hy[keeps]
    hy_cts <- merge(hy_cts, hy,all = TRUE, by = c("CHROM", "POS","mrnaID","geneID"))
  }
  #trans test (strict - only analyze genes with sites alternatively homozygous in parents and het in hybrids)
  #mean_cts <- merge(hy_cts, alterntae_homozygous_in_parents_positions, by = c("CHROM", "POS","mrnaID","geneID"))
  #trans test (lenient - analyze genes with any het sites in hybrids)
  #mean_cts <- hy_cts
  for(hy_df in c(5:(length(hy_inds)+4)))
  {
    fish_table <- hy_cts[,c(1:4,hy_df)]
    fish_table$Geneid <- fish_table$mrnaID
    fish_table <- merge(fish_table, p_cts, by = c("Geneid"))
    fish_table <- cbind(fish_table, colsplit(fish_table[,6], ";", c("mom_cts", "dad_cts")))
    ind <- strsplit(colnames(fish_table)[6],"_")[[1]][3]
    Geneids <- c()
    CHROMs <- c()
    POSs <- c()
    mrnaIDs <- c()
    geneIDs <- c()
    Hs <- c()
    Ps <- c()
    P_H_ratios <- c()
    fishers_exacts <- c()
    for(row_num in c(1:nrow(fish_table)))
      #row_num <- 1
    {
      Geneids <- c(Geneids,fish_table$Geneid[row_num])
      CHROMs  <- c(CHROMs ,fish_table$CHROM[row_num] )
      POSs    <- c(POSs   ,fish_table$POS[row_num]   )
      mrnaIDs <- c(mrnaIDs,fish_table$mrnaID[row_num])
      geneIDs <- c(geneIDs,fish_table$geneID[row_num])
      Hs <- c(Hs,(fish_table$mom_cts[row_num] /fish_table$dad_cts[row_num]))
      Ps <- c(Ps,(fish_table$p1_mean_ct[row_num] /fish_table$p2_mean_ct[row_num]))
      P_H_ratios <- c(P_H_ratios, ((fish_table$p1_mean_ct[row_num] /fish_table$p2_mean_ct[row_num]))/((fish_table$mom_cts[row_num] /fish_table$dad_cts[row_num])))
      test <- matrix(c(fish_table$mom_cts[row_num],fish_table$p1_mean_ct[row_num],fish_table$dad_cts[row_num],fish_table$p2_mean_ct[row_num]), ncol=2)
      if(any(is.na(test)))
      {
        fishers_exacts <- c(fishers_exacts,"NA")
      }else
      {
        fish <- fisher.test(test)
        fishers_exacts <- c(fishers_exacts,fish$p.value)
      }
    }
    if(hy_df == 5)
    {
      hy_h_name <- paste("Hr",ind, sep = "_")
      hy_p_name <- paste("Pr",ind, sep = "_")
      hy_ph_name <- paste("P_H_ratio",ind, sep = "_")
      hy_fp_name <- paste("fishers_exact_p",ind, sep = "_")
      fisher1 <- data.frame (Geneid        = Geneids,
                             CHROM         = CHROMs,
                             POS           = POSs ,
                             mrnaID        = mrnaIDs ,
                             geneID        = geneIDs,
                             Hr             = Hs ,
                             Pr             = Ps ,
                             P_H_ratio     = P_H_ratios ,
                             fishers_exact_p = fishers_exacts,
                             stringsAsFactors=FALSE)
      colnames(fisher1)[colnames(fisher1)=="Hr"] <- hy_h_name
      colnames(fisher1)[colnames(fisher1)=="Pr"] <- hy_p_name
      colnames(fisher1)[colnames(fisher1)=="P_H_ratio"] <- hy_ph_name
      colnames(fisher1)[colnames(fisher1)=="fishers_exact_p"] <- hy_fp_name
    }else
    {
      hy_h_name <- paste("Hr",ind, sep = "_")
      hy_p_name <- paste("Pr",ind, sep = "_")
      hy_ph_name <- paste("P_H_ratio",ind, sep = "_")
      hy_fp_name <- paste("fishers_exact_p",ind, sep = "_")
      fisher2 <- data.frame (Geneid        = Geneids,
                             CHROM         = CHROMs,
                             POS           = POSs ,
                             mrnaID        = mrnaIDs ,
                             geneID        = geneIDs,
                             Hr             = Hs ,
                             Pr             = Ps ,
                             P_H_ratio     = P_H_ratios ,
                             fishers_exact_p = fishers_exacts,
                             stringsAsFactors=FALSE)
      colnames(fisher2)[colnames(fisher2)=="Hr"] <- hy_h_name
      colnames(fisher2)[colnames(fisher2)=="Pr"] <- hy_p_name
      colnames(fisher2)[colnames(fisher2)=="P_H_ratio"] <- hy_ph_name
      colnames(fisher2)[colnames(fisher2)=="fishers_exact_p"] <- hy_fp_name
      fisher1 <- merge(fisher1, fisher2, by = c("CHROM", "POS","mrnaID","geneID","Geneid"))
    }
  }
  head(fisher1)
  fisher1$tag = paste(fisher1$CHROM, fisher1$mrnaID, sep = ";")
  make.true.NA <- function(x) if(is.character(x)||is.factor(x)){
    is.na(x) <- x=="NA"; x} else {
      x}
  fisher1[] <- lapply(fisher1, make.true.NA)
  tags <- c()
  H_sign <- c()
  P_sign <- c()
  sig_fisher <- c()
  #vector.is.empty <- function(x) return(length(x) ==0 )
  for(tag in unique(fisher1$tag))
    #tag <- "NW_015150454.1;XM_015398283.1;LOC107099927"
  {
    #H_sign <- c()
    chr_table <- fisher1[which(fisher1$tag == tag),]
    H_t <- chr_table[ , grepl( "Hr_" , names( chr_table ) ) ]
    H_t <- H_t[ colSums(is.na(H_t)) == 0]
    P_t <- chr_table[ , grepl( "Pr_" , names( chr_table ) ) ]
    P_t <- P_t[ colSums(is.na(P_t)) == 0]
    fishp <- chr_table[ , grepl( "fishers_exact_p_" , names( chr_table ) ) ]
    fishp <- fishp[ colSums(is.na(fishp)) == 0]
    if(nrow(H_t) >0 & nrow(P_t) >0 & nrow(fishp) >0 &length(colnames(H_t)) >0 & length(colnames(P_t))  >0 & length(colnames(fishp))  >0)
    {
      tags <- c(tags,tag)
      if (all(H_t>1))
      {
        H_sign <- c(H_sign, "+")
      }
      if (all(H_t<1))
      {
        H_sign <- c(H_sign, "-")
      }
      if (any(H_t>1) & any(H_t<1) | any(H_t==1))
      {
        H_sign <- c(H_sign, "ambiguous")
      }
      if (all(P_t>1))
      {
        P_sign <- c(P_sign, "+")
      }
      if (all(P_t<1))
      {
        P_sign <- c(P_sign, "-")
      }
      if (any(P_t>1) & any(P_t<1) | any(P_t==1))
      {
        P_sign <- c(P_sign, "ambiguous")
      }
      if (any(fishp>0.05))
      {
        sig_fisher <- c(sig_fisher, "no")
      }
      if (all(fishp<=0.05))
      {
        sig_fisher <- c(sig_fisher, "yes")
      }
    }
    #if(length(H_sign) != length(tag))
    #{
    #  break
    #}
  }
  trans_test <- data.frame (tag= tags,
                            H_sign = H_sign,
                            P_sign =P_sign,
                            sig_fisher=sig_fisher,
                            stringsAsFactors=FALSE)
  #mean_cts$hy_mom_mean <- rowMeans(subset(mean_cts, select = c(grep("mom", names(mean_cts), value = TRUE))), na.rm =TRUE)
  #mean_cts$hy_dad_mean <- rowMeans(subset(mean_cts, select = c(grep("dad", names(mean_cts), value = TRUE))), na.rm =TRUE)
  #mean_cts$Geneid <- mean_cts$mrnaID
  #keeps <- c("CHROM", "POS","mrnaID","geneID","Geneid","hy_mom_mean","hy_dad_mean")
  #mean_cts <- mean_cts[keeps]
  #mean_cts <- merge(mean_cts, p_cts, by = c("Geneid"))
  #
  #Geneids <- c()
  #CHROMs <- c()
  #POSs <- c()
  #mrnaIDs <- c()
  #geneIDs <- c()
  #Hs <- c()
  #Ps <- c()
  #P_H_ratios <- c()
  #fishers_exacts <- c()
  #
  #for(row_num in c(1:nrow(mean_cts)))
  #  #row_num <- 1
  #{
  #  Geneids <- c(Geneids,mean_cts$Geneid[row_num])
  #  CHROMs  <- c(CHROMs ,mean_cts$CHROM[row_num] )
  #  POSs    <- c(POSs   ,mean_cts$POS[row_num]   )
  #  mrnaIDs <- c(mrnaIDs,mean_cts$mrnaID[row_num])
  #  geneIDs <- c(geneIDs,mean_cts$geneID[row_num])
  #  Hs <- c(Hs,(mean_cts$hy_mom_mean[row_num] /mean_cts$hy_dad_mean[row_num]))
  #  Ps <- c(Ps,(mean_cts$p1_mean_ct[row_num] /mean_cts$p2_mean_ct[row_num]))
  #  P_H_ratios <- c(P_H_ratios, ((mean_cts$p1_mean_ct[row_num] /mean_cts$p2_mean_ct[row_num]))/((mean_cts$hy_mom_mean[row_num] /mean_cts$hy_dad_mean[row_num])))
  #
  #  test <- matrix(c(mean_cts$hy_mom_mean[row_num],mean_cts$p1_mean_ct[row_num],mean_cts$hy_dad_mean[row_num],mean_cts$p2_mean_ct[row_num]), ncol=2)
  #  if(any(is.na(test)))
  #  {
  #    fishers_exacts <- c(fishers_exacts,"NA")
  #  }else
  #  {
  #    fish <- fisher.test(test)
  #    fishers_exacts <- c(fishers_exacts,fish$p.value)
  #  }
  #
  #}
  #
  #
  #mean_cts <- data.frame (Geneid        = Geneids,
  #                        CHROM         = CHROMs,
  #                        POS           = POSs ,
  #                        mrnaID        = mrnaIDs ,
  #                        geneID        = geneIDs,
  #                        H             = Hs ,
  #                        P             = Ps ,
  #                        P_H_ratio     = P_H_ratios ,
  #                        fishers_exact_p = fishers_exacts,
  #                        stringsAsFactors=FALSE)
  #
  #mean_cts$tag = paste(mean_cts$CHROM, mean_cts$mrnaID, sep = ";")
  #mean_cts <- na.omit(mean_cts)
  #tags <- c()
  #H_sign <- c()
  #P_sign <- c()
  #sig_fisher <- c()
  #
  ##vector.is.empty <- function(x) return(length(x) ==0 )
  #
  #for(tag in unique(mean_cts$tag))
  #  #tag <- "NW_015152288.1;XM_015404710.1;ptrh2"
  #{
  #  #H_sign <- c()
  #  chr_table <- mean_cts[which(mean_cts$tag == tag),]
  #  tags <- c(tags,tag)
  #
  #  if (all(chr_table$H>1))
  #  {
  #    H_sign <- c(H_sign, "+")
  #  }
  #  if (all(chr_table$H<1))
  #  {
  #    H_sign <- c(H_sign, "-")
  #  }
  #  if (any(chr_table$H>1) & any(chr_table$H<1) | any(chr_table$H==1))
  #  {
  #    H_sign <- c(H_sign, "ambiguous")
  #  }
  #  if (all(chr_table$P>1))
  #  {
  #    P_sign <- c(P_sign, "+")
  #  }
  #  if (all(chr_table$P<1))
  #  {
  #    P_sign <- c(P_sign, "-")
  #  }
  #  if (any(chr_table$P>1) & any(chr_table$P<1) | any(chr_table$P==1))
  #  {
  #    P_sign <- c(P_sign, "ambiguous")
  #  }
  #  if (any(chr_table$fishers_exact_p>0.05))
  #  {
  #    sig_fisher <- c(sig_fisher, "no")
  #  }
  #  if (all(chr_table$fishers_exact_p<=0.05))
  #  {
  #    sig_fisher <- c(sig_fisher, "yes")
  #  }
  #
  #
  #  #if (any(is.na(chr_table$P)))
  #  #{
  #  #  P_sign <- c(P_sign, "NA")
  #  #}
  #  #if (any(is.na(chr_table$H)))
  #  #{
  #  #  H_sign <- c(H_sign, "NA")
  #  #}
  #  #
  #  #if (any(is.na(chr_table$fishers_exact_p)))
  #  #{
  #  #  sig_fisher <- c(sig_fisher, "NA")
  #  #}
  #
  #}
  #
  #
  #trans_test <- data.frame (tag= tags,
  #                          H_sign = H_sign,
  #                          P_sign =P_sign,
  #                          sig_fisher=sig_fisher,
  #                          stringsAsFactors=FALSE)
  #
  trans_test <-cbind(trans_test, colsplit(trans_test$tag, ";", c("CHROM", "related_accession","gene")))
  head(trans_test)
  head(inheritance_and_ase)
  # number of alternate homozygous sites between p1 and p2 (if using strict trans test)
  nrow(trans_test)
  # number of genes showing ase in hybrids not parents
  nrow(inheritance_and_ase[which (inheritance_and_ase$hyAse_any_h_no_p == "yes"),])
  nrow(inheritance_and_ase[which (inheritance_and_ase$hyAse_any_h_no_p_alt_hom_in_p == "yes"),])
  inheritance_and_ase <- merge(trans_test, inheritance_and_ase, all = TRUE, by = c("related_accession"))
  nrow(inheritance_and_ase)
  
  
  
  #### define cis/trans etc
  
  # strict (heterozygous sites in hybrids that are alternatively homozygous sites in parents)
  
  cis <- inheritance_and_ase[which(inheritance_and_ase$padj_p1p2 <= 0.05 & 
                                     inheritance_and_ase$hyAse_all_h_no_p_alt_hom_in_p == "yes" &
                                     inheritance_and_ase$sig_fisher == "no"),]
  
  trans <- inheritance_and_ase[which(inheritance_and_ase$padj_p1p2 <= 0.05 & 
                                       inheritance_and_ase$hyAse_all_h_no_p_alt_hom_in_p != "yes" &
                                       inheritance_and_ase$sig_fisher == "yes"),]
  
  cis_plus_trans <- inheritance_and_ase[which(inheritance_and_ase$padj_p1p2 <= 0.05 & 
                                                inheritance_and_ase$hyAse_all_h_no_p_alt_hom_in_p == "yes" &
                                                inheritance_and_ase$sig_fisher == "yes" &
                                                inheritance_and_ase$H_sign != "ambiguous" &
                                                inheritance_and_ase$H_sign == inheritance_and_ase$P_sign),]
  
  cis_by_trans <- inheritance_and_ase[which(inheritance_and_ase$padj_p1p2 <= 0.05 & 
                                              inheritance_and_ase$hyAse_all_h_no_p_alt_hom_in_p == "yes" &
                                              inheritance_and_ase$sig_fisher == "yes" &
                                              inheritance_and_ase$H_sign != "ambiguous" &
                                              inheritance_and_ase$H_sign != inheritance_and_ase$P_sign),]
  
  comp <-  inheritance_and_ase[which(inheritance_and_ase$padj_p1p2 > 0.05 & 
                                       inheritance_and_ase$hyAse_all_h_no_p_alt_hom_in_p == "yes" &
                                       inheritance_and_ase$sig_fisher == "yes"),]
  
  mis <- inheritance_and_ase[which(inheritance_and_ase$padj_ph <= 0.05),]
  
  mis_and_ase <- inheritance_and_ase[which(inheritance_and_ase$padj_ph <= 0.05 & 
                                             inheritance_and_ase$hyAse_all_h_no_p_alt_hom_in_p == "yes"),]
  
  
  parent_pops1_name_s <- c(parent_pops1_name_s,parent_pops1_name)
  parent_pops2_name_s <- c(parent_pops2_name_s,parent_pops2_name)
  hy_name_s           <- c(hy_name_s,hy_name)
  stage_s             <- c(stage_s,stage)
  
  st_informative_genes_s <- c(st_informative_genes_s,nrow(unique(inheritance_and_ase[which(inheritance_and_ase$hyAse_all_h_no_p_alt_hom_in_p == "yes" |inheritance_and_ase$hyAse_all_h_no_p_alt_hom_in_p == "no"),])))
  st_cis_s               <- c(st_cis_s, length(unique(cis$related_accession)))
  st_trans_s             <- c(st_trans_s, length(unique(trans$related_accession)))
  st_cis_plus_trans_s    <- c(st_cis_plus_trans_s, length(unique(cis_plus_trans$related_accession)))
  st_cis_by_trans_s      <- c(st_cis_by_trans_s, length(unique(cis_by_trans$related_accession)))
  st_comp_s              <- c(st_comp_s, length(unique(comp$related_accession)))
  st_mis_s               <- c(st_mis_s,length(unique(mis$related_accession)))
  st_mis_ase_s           <- c(st_mis_ase_s,length(unique(mis_and_ase$related_accession)))                 
  st_mis_ase_comp_s      <- c(st_mis_ase_comp_s ,length(unique(intersect(mis_and_ase$related_accession,comp$related_accession))))     
  st_mis_ase_cis_s       <- c(st_mis_ase_cis_s ,length(unique(intersect(mis_and_ase$related_accession,cis$related_accession))))      
  
  st_cis_genes               <- c(st_cis_genes, paste(unique(cis$related_accession), sep = ";",collapse=";"))
  st_trans_genes             <- c(st_trans_genes, paste(unique(trans$related_accession), sep = ";",collapse=";"))
  st_cis_plus_trans_genes    <- c(st_cis_plus_trans_genes, paste(unique(cis_plus_trans$related_accession), sep = ";",collapse=";"))
  st_cis_by_trans_genes      <- c(st_cis_by_trans_genes, paste(unique(cis_by_trans$related_accession), sep = ";",collapse=";"))
  st_comp_genes              <- c(st_comp_genes, paste(unique(comp$related_accession), sep = ";",collapse=";"))
  st_mis_ase_comp_genes      <- c(st_mis_ase_comp_genes ,paste(unique(intersect(mis_and_ase$related_accession,comp$related_accession)), sep = ";",collapse=";"))     
  st_mis_ase_cis_genes       <- c(st_mis_ase_cis_genes ,paste(unique(intersect(mis_and_ase$related_accession,cis$related_accession)), sep = ";",collapse=";")) 
  st_mis_trans_genes<- c(st_mis_trans_genes,paste(unique(intersect(mis$related_accession,trans$related_accession)), sep = ";",collapse=";"))
  
  
  # lenient (all heterozygous sites in hybrids)
  
  cis <- inheritance_and_ase[which(inheritance_and_ase$padj_p1p2 <= 0.05 & 
                                     inheritance_and_ase$hyAse_all_h_no_p == "yes" &
                                     inheritance_and_ase$sig_fisher == "no"),]
  
  trans <- inheritance_and_ase[which(inheritance_and_ase$padj_p1p2 <= 0.05 & 
                                       inheritance_and_ase$hyAse_all_h_no_p != "yes" &
                                       inheritance_and_ase$sig_fisher == "yes"),]
  
  cis_plus_trans <- inheritance_and_ase[which(inheritance_and_ase$padj_p1p2 <= 0.05 & 
                                                inheritance_and_ase$hyAse_all_h_no_p == "yes" &
                                                inheritance_and_ase$sig_fisher == "yes" &
                                                inheritance_and_ase$H_sign != "ambiguous" &
                                                inheritance_and_ase$H_sign == inheritance_and_ase$P_sign),]
  
  cis_by_trans <- inheritance_and_ase[which(inheritance_and_ase$padj_p1p2 <= 0.05 & 
                                              inheritance_and_ase$hyAse_all_h_no_p == "yes" &
                                              inheritance_and_ase$sig_fisher == "yes" &
                                              inheritance_and_ase$H_sign != "ambiguous" &
                                              inheritance_and_ase$H_sign != inheritance_and_ase$P_sign),]
  
  comp <-  inheritance_and_ase[which(inheritance_and_ase$padj_p1p2 > 0.05 & 
                                       inheritance_and_ase$hyAse_all_h_no_p == "yes" &
                                       inheritance_and_ase$sig_fisher == "yes"),]
  
  mis <- inheritance_and_ase[which(inheritance_and_ase$padj_ph <= 0.05),]
  
  mis_and_ase <- inheritance_and_ase[which(inheritance_and_ase$padj_ph <= 0.05 & 
                                             inheritance_and_ase$hyAse_all_h_no_p == "yes"),]
  
  
  informative_genes_s <- c(informative_genes_s,nrow(unique(inheritance_and_ase[which(inheritance_and_ase$hyAse_all_h_no_p == "yes" |inheritance_and_ase$hyAse_all_h_no_p == "no"),])))
  cis_s               <- c(cis_s, length(unique(cis$related_accession)))
  trans_s             <- c(trans_s, length(unique(trans$related_accession)))
  cis_plus_trans_s    <- c(cis_plus_trans_s, length(unique(cis_plus_trans$related_accession)))
  cis_by_trans_s      <- c(cis_by_trans_s, length(unique(cis_by_trans$related_accession)))
  comp_s              <- c(comp_s, length(unique(comp$related_accession)))
  mis_s               <- c(mis_s,length(unique(mis$related_accession)))
  mis_ase_s           <- c(mis_ase_s,length(unique(mis_and_ase$related_accession))) ##########                  
  mis_ase_comp_s      <- c(mis_ase_comp_s ,length(unique(intersect(mis_and_ase$related_accession,comp$related_accession))))     
  mis_ase_cis_s       <- c(mis_ase_cis_s ,length(unique(intersect(mis_and_ase$related_accession,cis$related_accession))))      
  
  cis_genes               <- c(cis_genes, paste(unique(cis$related_accession), sep = ";",collapse=";"))
  trans_genes             <- c(trans_genes, paste(unique(trans$related_accession), sep = ";",collapse=";"))
  cis_plus_trans_genes    <- c(cis_plus_trans_genes, paste(unique(cis_plus_trans$related_accession), sep = ";",collapse=";"))
  cis_by_trans_genes      <- c(cis_by_trans_genes, paste(unique(cis_by_trans$related_accession), sep = ";",collapse=";"))
  comp_genes              <- c(comp_genes, paste(unique(comp$related_accession), sep = ";",collapse=";"))
  mis_ase_comp_genes      <- c(mis_ase_comp_genes ,paste(unique(intersect(mis_and_ase$related_accession,comp$related_accession)), sep = ";",collapse=";"))     
  mis_ase_cis_genes       <- c(mis_ase_cis_genes ,paste(unique(intersect(mis_and_ase$related_accession,cis$related_accession)), sep = ";",collapse=";")) 
  mis_trans_genes     <- c(mis_trans_genes,paste(unique(intersect(mis$related_accession,trans$related_accession)), sep = ";",collapse=";"))
  
  
  
  # super-lenient (any heterozygous site, no trans test)
  
  cis <- inheritance_and_ase[which(inheritance_and_ase$padj_p1p2 <= 0.05 & 
                                     inheritance_and_ase$hyAse_all_h_no_p == "yes"),]
  
  trans <- inheritance_and_ase[which(inheritance_and_ase$padj_p1p2 <= 0.05 & 
                                       inheritance_and_ase$hyAse_all_h_no_p != "yes"),]
  
  comp <-  inheritance_and_ase[which(inheritance_and_ase$padj_p1p2 > 0.05 & 
                                       inheritance_and_ase$hyAse_all_h_no_p == "yes"),]
  
  mis <- inheritance_and_ase[which(inheritance_and_ase$padj_ph <= 0.05),]
  
  mis_and_ase <- inheritance_and_ase[which(inheritance_and_ase$padj_ph <= 0.05 & 
                                             inheritance_and_ase$hyAse_all_h_no_p == "yes"),]
  
  sl_informative_genes_s <- c(sl_informative_genes_s,nrow(unique(inheritance_and_ase[which(inheritance_and_ase$hyAse_all_h_no_p == "yes" |inheritance_and_ase$hyAse_all_h_no_p == "no"),])))
  sl_cis_s               <- c(sl_cis_s, length(unique(cis$related_accession)))
  sl_trans_s             <- c(sl_trans_s, length(unique(trans$related_accession)))
  sl_comp_s              <- c(sl_comp_s, length(unique(comp$related_accession)))
  sl_mis_s               <- c(sl_mis_s,length(unique(mis$related_accession)))
  sl_mis_ase_s           <- c(sl_mis_ase_s,length(unique(mis_and_ase$related_accession)))                   
  sl_mis_ase_comp_s      <- c(sl_mis_ase_comp_s ,length(unique(intersect(mis_and_ase$related_accession,comp$related_accession))))     
  sl_mis_ase_cis_s       <- c(sl_mis_ase_cis_s ,length(unique(intersect(mis_and_ase$related_accession,cis$related_accession))))      
  
  sl_cis_genes               <- c(sl_cis_genes, paste(unique(cis$related_accession), sep = ";",collapse=";"))
  sl_trans_genes             <- c(sl_trans_genes, paste(unique(trans$related_accession), sep = ";",collapse=";"))
  sl_comp_genes              <- c(sl_comp_genes, paste(unique(comp$related_accession), sep = ";",collapse=";"))
  sl_mis_ase_comp_genes      <- c(sl_mis_ase_comp_genes ,paste(unique(intersect(mis_and_ase$related_accession,comp$related_accession)), sep = ";",collapse=";"))     
  sl_mis_ase_cis_genes       <- c(sl_mis_ase_cis_genes ,paste(unique(intersect(mis_and_ase$related_accession,cis$related_accession)), sep = ";",collapse=";")) 
  sl_mis_trans_genes<- c(sl_mis_trans_genes,paste(unique(intersect(mis$related_accession,trans$related_accession)), sep = ";",collapse=";"))
  
  
  print(parent_pops1_name)
  print(parent_pops2_name)
  print(hy_name)
  print(stage)
  
} 

strict_table <- data.frame(parent_pops1_name_s = parent_pops1_name_s ,
                           parent_pops2_name_s = parent_pops2_name_s ,
                           hy_name_s           = hy_name_s           ,
                           stage_s             = stage_s             ,
                           informative_genes  = st_informative_genes_s    ,
                           cis                = st_cis_s                  ,
                           trans              = st_trans_s                ,
                           cis_plus_trans     = st_cis_plus_trans_s       ,
                           cis_by_trans       = st_cis_by_trans_s         ,
                           comp               = st_comp_s                 ,
                           mis                = st_mis_s                  ,
                           mis_ase            = st_mis_ase_s              ,
                           mis_ase_comp       = st_mis_ase_comp_s         ,
                           mis_ase_cis        = st_mis_ase_cis_s          ,
                           cis_genes            = st_cis_genes            ,
                           trans_genes          = st_trans_genes          ,
                           cis_plus_trans_genes = st_cis_plus_trans_genes ,
                           cis_by_trans_genes   = st_cis_by_trans_genes   ,
                           comp_genes           = st_comp_genes           ,
                           mis_ase_comp_genes   = st_mis_ase_comp_genes   ,
                           mis_ase_cis_genes    = st_mis_ase_cis_genes    ,
                           mis_trans_genes    = st_mis_trans_genes,
                           stringsAsFactors = FALSE)


lenient_table <- data.frame(parent_pops1_name_s = parent_pops1_name_s ,
                            parent_pops2_name_s = parent_pops2_name_s ,
                            hy_name_s           = hy_name_s           ,
                            stage_s             = stage_s             ,
                            informative_genes  = informative_genes_s    ,
                            cis                = cis_s                  ,
                            trans              = trans_s                ,
                            cis_plus_trans     = cis_plus_trans_s       ,
                            cis_by_trans       = cis_by_trans_s         ,
                            comp               = comp_s                 ,
                            mis                = mis_s                  ,
                            mis_ase            = mis_ase_s              ,
                            mis_ase_comp       = mis_ase_comp_s         ,
                            mis_ase_cis        = mis_ase_cis_s          ,
                            cis_genes            = cis_genes            ,
                            trans_genes          = trans_genes          ,
                            cis_plus_trans_genes = cis_plus_trans_genes ,
                            cis_by_trans_genes   = cis_by_trans_genes   ,
                            comp_genes           = comp_genes           ,
                            mis_ase_comp_genes   = mis_ase_comp_genes   ,
                            mis_ase_cis_genes    = mis_ase_cis_genes    ,
                            mis_trans_genes    = mis_trans_genes,
                            stringsAsFactors = FALSE)

super_lenient_table <- data.frame(parent_pops1_name_s = parent_pops1_name_s ,
                                  parent_pops2_name_s = parent_pops2_name_s ,
                                  hy_name_s           = hy_name_s           ,
                                  stage_s             = stage_s             ,
                                  informative_genes  = sl_informative_genes_s    ,
                                  cis                = sl_cis_s                  ,
                                  trans              = sl_trans_s                ,
                                  comp               = sl_comp_s                 ,
                                  mis                = sl_mis_s                  ,
                                  mis_ase            = sl_mis_ase_s              ,
                                  mis_ase_comp       = sl_mis_ase_comp_s         ,
                                  mis_ase_cis        = sl_mis_ase_cis_s          ,
                                  cis_genes            = sl_cis_genes            ,
                                  trans_genes          = sl_trans_genes          ,
                                  comp_genes           = sl_comp_genes           ,
                                  mis_ase_comp_genes   = sl_mis_ase_comp_genes   ,
                                  mis_ase_cis_genes    = sl_mis_ase_cis_genes    ,
                                  mis_trans_genes    = sl_mis_trans_genes,
                                  stringsAsFactors = FALSE)

#write.table(strict_table,        "C:/Users/jmcgirr/Documents/all_2018_samples/ase_data/ase_and_mse/strict_ase_unphased_mbased_transtest_by_individual_no_ase_all_parents.txt", row.names = FALSE, quote = FALSE, sep ="\t")
#write.table(lenient_table,       "C:/Users/jmcgirr/Documents/all_2018_samples/ase_data/ase_and_mse/THE_ONE_lenient_ase_unphased_mbased_transtest_by_individual_no_ase_all_parents.txt", row.names = FALSE, quote = FALSE, sep ="\t")
#write.table(super_lenient_table, "C:/Users/jmcgirr/Documents/all_2018_samples/ase_data/ase_and_mse/super_lenient_ase_unphased_mbased_transtest_by_individual_no_ase_all_parents.txt", row.names = FALSE, quote = FALSE, sep ="\t")



  
  #con = yellow ="#f6c700"
  #cis = orange = "#E77200"
  #comp = pink = "#FF00CC"
  #trans = black = "#000000"
  #over_dom = red = #bd1e24
  #under_dom = blue = #0067a7
  #
  
  cols <- 	c("#f6c700", "#bd1e24", "#0067a7","#000000","#E77200","#FF00CC")
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
  
  
  #tiff(paste("D:/Martin Lab/rna_2018/ASE/plots/",hy_name, "_",stage,"ase.tiff", sep = ""), width = 6, height = 6, units = 'in', res = 300)
  #par(mfrow=c(1,1))
  plot(ase_type$lfc_p1p2, ase_type$lfc_ph, col = cols_ase, 
       pch = c(16,17)[as.factor(ase_type$ase)],
       cex = c(1, 1, 1, 1,1,1)[ase_type$ase_type],cex.axis=1.5, ylab = "", xlab = "", ylim = c(-2,1.6),xlim = c(-2,1.6))
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
  #dev.off()
  
  cols_b <- c("#f6c700", "#FF00CC","#0067a7", "#bd1e24","#E77200","#000000")
  labs <- c("additive",paste(parent_pops1_name,"dominant", sep = " "),paste(parent_pops1_name,"dominant", sep = " "),"overdominant","underdominant", "conserved")
  bars <- c(((nrow(ase_type[which(ase_type$ase_type == 1),]) / total_genes) *100), ((nrow(ase_type[which(ase_type$ase_type == 6),]) / total_genes) *100), ((nrow(ase_type[which(ase_type$ase_type == 3),]) / total_genes) *100),((nrow(ase_type[which(ase_type$ase_type == 2),]) / total_genes) *100),((nrow(ase_type[which(ase_type$ase_type == 5),]) / total_genes) *100),((nrow(ase_type[which(ase_type$ase_type == 4),]) / total_genes) *100))
  
  #tiff("D:/Martin Lab/RNA-seq/axm/post_reviews/ASE/plots/17dpf_bar.tiff", width = 6, height = 4, units = 'in', res = 1000)
  b1 <- barplot(bars, las =1, horiz = TRUE, font.axis = 2, cex.axis=1.5, col = cols_b, xlim = c(0,b1$))
  #barplot(bars, las =1, horiz = TRUE, font.axis = 2, cex.axis=1.5, col = cols_b)
  
  offset <- 10
  text(y=b1[1,], x=((nrow(ase_type[which(ase_type$ase_type == 1),]) / total_genes) *100)+offset,con, font = 2, cex = 1)
  text(y=b1[2,], x=((nrow(ase_type[which(ase_type$ase_type == 6),]) / total_genes) *100)+offset,comp, font = 2, cex = 1)
  text(y=b1[3,], x=((nrow(ase_type[which(ase_type$ase_type == 3),]) / total_genes) *100)+offset,mis_under, font = 2, cex = 1)
  text(y=b1[4,], x=((nrow(ase_type[which(ase_type$ase_type == 2),]) / total_genes) *100)+offset,mis_over, font = 2, cex = 1)
  text(y=b1[5,], x=((nrow(ase_type[which(ase_type$ase_type == 5),]) / total_genes) *100)+offset,trans, font = 2, cex = 1)
  text(y=b1[6,], x=((nrow(ase_type[which(ase_type$ase_type == 4),]) / total_genes) *100)+offset,cis, font = 2, cex = 1)
  
  
  
  
  
 
#####
###############################################
###############################################
############ visualize ase at cool genes ######
###############################################

library(seqinr)
library(reshape2)
st <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/ase_data/ase_and_mse/super_lenient_ase_unphased_mbased_transtest_by_individual_no_ase_all_parents.txt", stringsAsFactors = FALSE, header = TRUE, sep = "\t")
#le <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/ase_data/ase_and_mse/lenient_ase.txt", stringsAsFactors = FALSE, header = TRUE, sep = "\t")
#sl <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/ase_data/ase_and_mse/super_lenient_ase.txt", stringsAsFactors = FALSE, header = TRUE, sep = "\t")

norm_cts <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/norm_cts_all_2018_no_seq_round_control.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
inheritance_comps <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/ase_data/ase_comparisons.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
master <- read.table("D:/Martin Lab/rna_2018/all_2018_samples/table_maker_master_outlier_rm.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
norm_cts <- cbind(norm_cts, colsplit(norm_cts$Geneid, ";", c("related_accession", "gene_name")))
mrna <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/ase_data/mrna.saf", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
mrna <- cbind(mrna, colsplit(mrna$GeneID, ";", c("related_accession", "gene_name")))
maternal_counts_dir <- "C:/Users/jmcgirr/Documents/all_2018_samples/ase_data/maternal_counts/"
gff <- read.delim("D:/Cyprinodon/ref_C_variegatus-1.0_scaffolds_no_header.gff3",header = FALSE, stringsAsFactors = FALSE, sep = "\t")


cis         <- c()  
comp        <- c()  
trans       <- c()  
mis_ase_comp<- c()  
mis_ase_cis <- c()  

for (row_num in c(1:nrow(st)))
{
  cis                      <- c(cis, strsplit(st$cis_genes[row_num],";")[[1]])
  comp                    <- c(comp, strsplit(st$comp_genes[row_num],";")[[1]])
  trans                  <- c(trans, strsplit(st$trans_genes[row_num],";")[[1]])
  mis_ase_comp    <- c(mis_ase_comp, strsplit(st$mis_ase_comp_genes[row_num],";")[[1]])
  mis_ase_cis     <- c(mis_ase_cis,  strsplit(st$mis_ase_cis_genes[row_num],";")[[1]])
}
#duplicated(ciss)
all <- c(cis, comp, trans, mis_ase_comp, mis_ase_cis)
dall <- data.frame(table(all))
dall <- dall[dall$Freq > 1,]

d1 <- data.frame(table(cis))
d1 <- d1[d1$Freq > 1,]
d2 <- data.frame(table(comp))
d2 <- d2[d2$Freq > 1,]
d3 <- data.frame(table(trans))
d3 <- d3[d3$Freq > 1,]
d4 <- data.frame(table(mis_ase_comp))
d4 <- d4[d4$Freq > 1,]
d5 <- data.frame(table(mis_ase_cis))
d5 <- d5[d5$Freq > 1,]


genes <- list(dall,d1, d2,d3,d4,d5)
genes
  
  
  
  
  #### visualize cool genes
  # "#E60012" red A
  # "#00A0E9" blue P
  # "#f6c700" yellow M
  # "#000000" black O
  # "#920783" purple AxP
  # "#FF6600" orange AxM
  # "#22AC38" green MxP
  # "#800000" dark red AxO
axm_cols <- c("#E60012","#f6c700","#FF6600")
mxa_cols <- c("#f6c700","#E60012","#FF6600")
axp_cols <- c("#E60012","#00A0E9","#920783")
pxa_cols <- c("#00A0E9","#E60012","#920783")
mxp_cols <- c("#f6c700","#00A0E9","#22AC38")


cool_gene <- "XM_015371309.1"
cool_p1 <- "CRPA"
cool_p2 <- "CRPP"
cool_h <- "CAxCP"
cool_stage <- "8dpf"
cool_cols <- axp_cols
pop_gen_comp <- "caxcp"

##### make cool gene plots #####
pop_gen_dna <- read.csv(paste("D:/Martin Lab/rna_2018/fst/dna/",pop_gen_comp,"_popgen_dna_stats.csv",sep = ""), header = TRUE)
fst <- read.table(paste("C:/Users/jmcgirr/Documents/all_2018_samples/fst_dna/",pop_gen_comp,"_fst",sep = ""), header = TRUE)
pop_gen_dna[pop_gen_dna<0] <- 0
pop_gen_dna <-  na.omit(pop_gen_dna)
lake <- strsplit(pop_gen_comp,split ="")[[1]][1]

taj_m <-    read.table(paste("C:/Users/jmcgirr/Documents/all_2018_samples/fst_dna/",lake,"m_taj_d.txt.Tajima.D", sep = ""), header = TRUE, stringsAsFactors = FALSE)
taj_a <-    read.table(paste("C:/Users/jmcgirr/Documents/all_2018_samples/fst_dna/",lake,"a_taj_d.txt.Tajima.D", sep = ""), header = TRUE, stringsAsFactors = FALSE)
taj_p <-    read.table(paste("C:/Users/jmcgirr/Documents/all_2018_samples/fst_dna/",lake,"p_taj_d.txt.Tajima.D", sep = ""), header = TRUE, stringsAsFactors = FALSE)
taj_m <- na.omit(taj_m)
taj_a <- na.omit(taj_a)
taj_p <- na.omit(taj_p)

p1_inds <- master[which(master$f1 == cool_p1 & master$stage == cool_stage),]
p2_inds <- master[which(master$f1 == cool_p2 & master$stage == cool_stage),]
hy_inds <- master[which(master$f1 == cool_h & master$stage == cool_stage),]
p1_inds <- p1_inds$sample
p2_inds <- p2_inds$sample
hy_inds <- hy_inds$sample
cool_genes_out_dir <- "C:/Users/jmcgirr/Documents/all_2018_samples/ase_data/plots/cool_genes/"
norm_cts_gene <- norm_cts[which(norm_cts$related_accession == cool_gene),]
ph_cts <- norm_cts_gene[c(p1_inds,p2_inds,hy_inds)]
ph_cts  <- as.numeric(ph_cts[1,])
inds <- c(p1_inds,p2_inds,hy_inds)
type_inds <- c(rep(cool_p1,length(p1_inds)),rep(cool_p2,length(p2_inds)),rep(cool_h,length(hy_inds)))
plot_counts <- data.frame(type_ind = type_inds, sample = inds, cts = ph_cts)
plot_counts$type_ind <- as.character(plot_counts$type_ind)
plot_counts$type_ind <- factor(plot_counts$type_ind, levels=unique(plot_counts$type_ind))
plot_counts$type_ind <- factor(plot_counts$type_ind, level = c(cool_p1,cool_h,cool_p2))
plot_title <- paste(cool_gene,"\n",norm_cts_gene$gene_name[1], sep = "")
cool_gene_out <- paste(cool_genes_out_dir,norm_cts_gene$gene_name[1],"_",cool_h,"_",cool_stage, sep = "")
#tiff(paste(cool_gene_out, "_gene_counts.tiff", sep = ""), width = 5, height = 6, units = 'in', res = 300)
  plot(plot_counts$type_ind,plot_counts$cts, ylab = "normalized counts", main = plot_title, border = "white")
  p1_cts <- plot_counts[which(plot_counts$type_ind == cool_p1),]
  points(rep(1,nrow(p1_cts)), p1_cts$cts, pch = 21,bg =cool_cols[1],  col = "black", cex = 2)
  p2_cts <- plot_counts[which(plot_counts$type_ind == cool_p2),]
  points(rep(3,nrow(p2_cts)), p2_cts$cts, pch = 21,bg =cool_cols[2],  col = "black", cex = 2)
  hy_cts <- plot_counts[which(plot_counts$type_ind == cool_h),]
  points(rep(2,nrow(hy_cts)), hy_cts$cts, pch = 21,bg =cool_cols[3],  col = "black", cex = 2)
#dev.off()
gene_chrom <- mrna[which(mrna$related_accession == cool_gene),]
pop_gen_plot_dna <- pop_gen_dna[which(pop_gen_dna$scaffold == gene_chrom$Chr[1] & pop_gen_dna$start > gene_chrom$Start[1] - 100000 & pop_gen_dna$end < gene_chrom$End[1] +100000),]
smooth_fst <-   smooth.spline(pop_gen_plot_dna$mid, pop_gen_plot_dna[,9], spar = .5)
smooth_dxy <-   smooth.spline(pop_gen_plot_dna$mid, pop_gen_plot_dna[,8], spar = .5)
smooth_p1_pi <- smooth.spline(pop_gen_plot_dna$mid, pop_gen_plot_dna[,6], spar = .5)
smooth_p2_pi <- smooth.spline(pop_gen_plot_dna$mid, pop_gen_plot_dna[,7], spar = .5)
#tiff(paste(cool_gene_out, "_pop_gen_dna.tiff", sep = ""), width = 5, height = 6, units = 'in', res = 300)
  plot(pop_gen_plot_dna$mid, pop_gen_plot_dna[,9],ylim = c(0,1),col = "white",xaxt = 'n', xlab = "", ylab = "", main = paste(cool_p1, "vs", cool_p2, sep = " "))
  title(xlab = paste(round(pop_gen_plot_dna$mid[nrow(pop_gen_plot_dna)] -pop_gen_plot_dna$mid[1], digits = -3)," bp",sep = ""), line=0.3, cex.lab=1)
  lines(smooth_fst, lwd = 2, col = "black")
  lines(smooth_dxy, lwd = 2, col = "grey")
  lines(smooth_p1_pi, lwd = 2, col = cool_cols[1])
  lines(smooth_p2_pi, lwd = 2, col = cool_cols[2])
  rect(gene_chrom$Start[1], -1000, gene_chrom$End[1], 1000, border = NA, col = col2alpha("darkblue",0.2))
  legend("topleft", inset=0,
  c("Fst", "Dxy", paste("Pi ",cool_p1, sep = ""),paste("Pi ",cool_p2, sep = "")),
  cex = 0.9, col=c("black","grey",cool_cols[1],cool_cols[2]), horiz=FALSE, bty="n", lty = 1, lwd=3)
#dev.off()
#tiff(paste(cool_gene_out, "_fst.tiff", sep = ""), width = 5, height = 6, units = 'in', res = 300)
  fst_plot <- fst[which(fst$CHROM == gene_chrom$Chr[1]),]
  plot(pop_gen_plot_dna$mid, pop_gen_plot_dna[,9],ylim = c(0,1),xaxt ='n',col = "white", xlab = "", ylab = "", main = paste(cool_p1, "vs", cool_p2, sep = " "))
  title(xlab = paste(round(pop_gen_plot_dna$mid[nrow(pop_gen_plot_dna)] -pop_gen_plot_dna$mid[1], digits = -3)," bp",sep = ""), line=0.3, cex.lab=1)
  points(fst_plot$POS, fst_plot$WEIR_AND_COCKERHAM_FST)
  rect(gene_chrom$Start[1], -1000, gene_chrom$End[1], 1000, border = NA, col = col2alpha("darkblue",0.2))
  best_cis_candidate <- fst_plot[which(fst_plot$CHROM == gene_chrom$Chr[1] & fst_plot$POS > (gene_chrom$Start[1] - 10000) & fst_plot$POS < (gene_chrom$End[1] + 10000) ),]
  best_cis_candidate <- best_cis_candidate[order(best_cis_candidate$WEIR_AND_COCKERHAM_FST, decreasing = TRUE),]
#dev.off()
#range01 <- function(x){(x-min(x))/(max(x)-min(x))}
#sw_m <- sweeps_m[which(sweeps_m$V1 == gene_chrom$Chr[1]),]
#smooth_m <-   smooth.spline(sw_m$V2, range01(sw_m$V4), spar = .2)
#sw_p <- sweeps_p[which(sweeps_p$V1 == gene_chrom$Chr[1]),]
#smooth_p <-   smooth.spline(sw_p$V2, range01(sw_p$V4), spar = .2)
#sw_a <- sweeps_a[which(sweeps_a$V1 == gene_chrom$Chr[1]),]
#smooth_a <-   smooth.spline(sw_a$V2, range01(sw_a$V4), spar = .2)
##tiff(paste(cool_gene_out, "_sweed.tiff", sep = ""), width = 5, height = 6, units = 'in', res = 300)
#  plot(sw_m$V2, range01(sw_m$V4),xlab = "", main = "sweeps",ylab = "CLR",xlim = c((gene_chrom$Start[1] -100000),(gene_chrom$End[1] +100000)), xaxt = 'n',col = "white")
#  lines(smooth_m, lwd = 2, col = "#f6c700")
#  lines(smooth_p, lwd = 2, col = "#00A0E9")
#  lines(smooth_a, lwd = 2, col = "#E60012")
#  legend("topleft", inset=0,
#  c("generalist","molluscivore", "scale-eater"),
#  cex = 0.9, fill=c("#E60012","#f6c700","#00A0E9"), horiz=FALSE, bty="n")
#  title(xlab = paste(round(pop_gen_plot_dna$mid[nrow(pop_gen_plot_dna)] -pop_gen_plot_dna$mid[1], digits = -3)," bp",sep = ""), line=0.3, cex.lab=1)
#  rect(gene_chrom$Start[1], -1000, gene_chrom$End[1], 1000, border = NA, col = col2alpha("darkblue",0.2))
##dev.off()
#tiff(paste(cool_gene_out, "_tajima.tiff", sep = ""), width = 5, height = 6, units = 'in', res = 300)
tj_m <- taj_m[which(taj_m$CHROM == gene_chrom$Chr[1]),]
smooth_m <-   smooth.spline(tj_m$BIN_START, tj_m$TajimaD, spar = .1)
tj_p <- taj_p[which(taj_p$CHROM == gene_chrom$Chr[1]),]
smooth_p <-   smooth.spline(tj_p$BIN_START, tj_p$TajimaD, spar = .1)
tj_a <- taj_a[which(taj_a$CHROM == gene_chrom$Chr[1]),]
smooth_a <-   smooth.spline(tj_a$BIN_START, tj_a$TajimaD, spar = .1)
plot(tj_p$BIN_START, tj_p$TajimaD,xlab = "", main = "Tajimas D",ylab = "D",xlim = c((gene_chrom$Start[1] -100000),(gene_chrom$End[1] +100000)), xaxt = 'n',col = "white")
lines(smooth_m, lwd = 2, col = "#f6c700")
lines(smooth_p, lwd = 2, col = "#00A0E9")
lines(smooth_a, lwd = 2, col = "#E60012")
legend("topleft", inset=0,
c("generalist","molluscivore", "scale-eater"),
cex = 0.9, fill=c("#E60012","#f6c700","#00A0E9"), horiz=FALSE, bty="n")
title(xlab = paste(round(pop_gen_plot_dna$mid[nrow(pop_gen_plot_dna)] -pop_gen_plot_dna$mid[1], digits = -3)," bp",sep = ""), line=0.3, cex.lab=1)
rect(gene_chrom$Start[1], -1000, gene_chrom$End[1], 1000, border = NA, col = col2alpha("darkblue",0.2))
#dev.off()
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

#tiff(paste(cool_gene_out,"_",ind, "_allele_counts.tiff", sep = ""), width = 6, height = 4.5, units = 'in', res = 300)
plot(allele_cts$position, c(1:length(allele_cts$position)),ylim = y_ax_range,
xlim = c(min(exons$V4)-300,max(exons$V5)+300),col = "white",xlab = "" ,ylab = "parental read counts", main = plot_title, xaxt = 'n')
title(xlab = paste(round(((max(exons$V5)+300)-(min(exons$V4)+300)), digits = -3)," bp", sep = ""), line=0.3, cex.lab=1)
abline(h =mid_point,col = col2alpha("darkblue",0.2), lwd = 2)
for(i in c(1:nrow(exons)))
{
rect(exons$V4[i], mid_point - ((y_ax_range[2]-y_ax_range[1]) /40), exons$V5[i], mid_point + ((y_ax_range[2]-y_ax_range[1])/40), border = NA, col = col2alpha("darkblue",0.2))
}
points(allele_cts$position, allele_cts[,4], pch = 21,bg =cool_cols[1],  col = "black", cex = 1.2)
points(allele_cts$position, allele_cts[,5], pch = 21,bg =cool_cols[2],  col = "black", cex = 1.2)
#dev.off()
}
ind <- hy_inds[[1]]
cts <- read.csv(paste(maternal_counts_dir,ind, "_parental_counts_features.txt",sep = ""),na.strings=c("","NA"), header = TRUE, stringsAsFactors = FALSE, sep = "\t") 
cts <- cbind(cts, colsplit(cts$mrnaID, ";", c("related_accession", "gene_name")))
cts <- cts[which(cts$chrom == gene_chrom$Chr[1] & cts$related_accession == cool_gene),]
keeps <- c("chrom","mrnaID","position")
momallele_cts <- cts[keeps]
for (ind in hy_inds)
{
cts <- read.csv(paste(maternal_counts_dir,ind, "_parental_counts_features.txt",sep = ""),na.strings=c("","NA"), header = TRUE, stringsAsFactors = FALSE, sep = "\t") 
cts <- cbind(cts, colsplit(cts$mrnaID, ";", c("related_accession", "gene_name")))
cts <- cts[which(cts$chrom == gene_chrom$Chr[1] & cts$related_accession == cool_gene),]
colnames(cts)[colnames(cts)=="momCount"] <- paste(ind,"_momCount", sep = "")
keeps <- c("chrom","mrnaID","position",paste(ind,"_momCount", sep = ""))
momallele_cts <- merge(momallele_cts, cts[keeps],all = TRUE, by = c("chrom","mrnaID","position"))
}
ind <- hy_inds[[1]]
cts <- read.csv(paste(maternal_counts_dir,ind, "_parental_counts_features.txt",sep = ""),na.strings=c("","NA"), header = TRUE, stringsAsFactors = FALSE, sep = "\t") 
cts <- cbind(cts, colsplit(cts$mrnaID, ";", c("related_accession", "gene_name")))
cts <- cts[which(cts$chrom == gene_chrom$Chr[1] & cts$related_accession == cool_gene),]
keeps <- c("chrom","mrnaID","position")
dadallele_cts <- cts[keeps]
for (ind in hy_inds)
{
cts <- read.csv(paste(maternal_counts_dir,ind, "_parental_counts_features.txt",sep = ""),na.strings=c("","NA"), header = TRUE, stringsAsFactors = FALSE, sep = "\t") 
cts <- cbind(cts, colsplit(cts$mrnaID, ";", c("related_accession", "gene_name")))
cts <- cts[which(cts$chrom == gene_chrom$Chr[1] & cts$related_accession == cool_gene),]
colnames(cts)[colnames(cts)=="dadCount"] <- paste(ind,"_dadCount", sep = "")
keeps <- c("chrom","mrnaID","position",paste(ind,"_dadCount", sep = ""))
dadallele_cts <- merge(dadallele_cts, cts[keeps],all = TRUE, by = c("chrom","mrnaID","position"))
}
allele_cts <- merge(momallele_cts, dadallele_cts, by = c("chrom", "mrnaID", "position"))
hy_cts <- na.omit(as.vector(as.matrix(allele_cts[,c(4:((2*length(hy_inds))+3))])))
#exons <- gff[grep(cool_gene, gff$V9),]
#exons <- exons[which(exons$V3 == "exon"),]
#exons <- exons[c("V1","V4","V5")]
mid_point <- max(hy_cts) - ((max(hy_cts) - min(hy_cts))/2)
y_ax_range <- c(min(hy_cts) -((max(hy_cts)-mid_point)*0.2),max(hy_cts) +((max(hy_cts)-mid_point)*0.1))
#tiff(paste(cool_gene_out,"_all_hybrids_", "_allele_counts.tiff", sep = ""), width = 6, height = 4.5, units = 'in', res = 300)
plot(allele_cts$position, c(1:length(allele_cts$position)),ylim = c(min(hy_cts) -((max(hy_cts)-mid_point)*0.2),max(hy_cts) +((max(hy_cts)-mid_point)*0.1)),
xlim = c(min(exons$V4)-300,max(exons$V5)+300),col = "white",xlab = "" ,ylab = "parental read counts", main = plot_title, xaxt = 'n')
title(xlab = paste(round(((max(exons$V5)+300)-(min(exons$V4)+300)), digits = -3)," bp", sep = ""), line=0.3, cex.lab=1)
abline(h =mid_point,col = col2alpha("darkblue",0.2), lwd = 2)
for(i in c(1:nrow(exons)))
{
rect(exons$V4[i], mid_point - ((y_ax_range[2]-y_ax_range[1]) /40), exons$V5[i], mid_point + ((y_ax_range[2]-y_ax_range[1]) /40), border = NA, col = col2alpha("darkblue",0.2))
}
for (ind in hy_inds)
{
h <- paste(ind,"_momCount", sep = "")
df <- momallele_cts[c("position",h)]

for(i in c(1:nrow(df)))
{
points(df$position[i], df[,2][i], pch = 21,bg =cool_cols[1],  col = "black", cex = 1.2)
}
}
for (ind in hy_inds)
{
h <- paste(ind,"_dadCount", sep = "")
df <- dadallele_cts[c("position",h)]

for(i in c(1:nrow(df)))
{
points(df$position[i], df[,2][i], pch = 21,bg =cool_cols[2],  col = "black", cex = 1.2)
}
}
legend("topleft", inset=0,
c(paste(cool_p1," allele", sep = ""),paste(cool_p2," allele", sep = "")),
cex = 0.9, fill=c(cool_cols[1],cool_cols[2]), horiz=FALSE, bty="n")
#dev.off()
######
head(best_cis_candidate)
plot_title

#####
###############################################
###############################################
############ DE vs fst/dxy ####################
###############################################
require(ggplot2)
mis <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/ase_data/fst_dxy_de_mrna.txt ", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
#mis1 <- merge(mis, sum_stats_final, by = c("parents", "f1", "stage"))
mis8 <- mis[which(mis$stage =="8dpf"),] 
mis2 <- mis[which(mis$stage =="48hpf"),]
cols <- c("#91648F","#808080","#5AA27C","#f6c700","#E77200","#FF00CC","#bd1e24","#0067a7")

#pdf("D:/Martin Lab/rna_2018/fst/dxy_vs_prop_mse.pdf",width=7,height=5)

plt <- ggplot(mis2, aes(x=dxy_dna, y=prop_misexpressed, color = cross_type)) + geom_point(aes(shape=lake),size = 4)
plt + theme_classic() +scale_color_manual(values = cols) +labs(x = "\ngenome-wide mean dxy", y = "% genes misexpressed in F1 hybrids\n") +guides(col = guide_legend(override.aes = list(shape = 19, size = 5))) +ggtitle("48 hpf misexpressed")+theme(plot.title = element_text(hjust = 0.5))
plt <- ggplot(mis8, aes(x=dxy_dna, y=prop_misexpressed, color = cross_type)) + geom_point(aes(shape=lake),size = 4)
plt + theme_classic() +scale_color_manual(values = cols) +labs(x = "\ngenome-wide mean dxy", y = "% genes misexpressed in F1 hybrids\n") +guides(col = guide_legend(override.aes = list(shape = 19, size = 5))) +ggtitle("8dpf misexpressed")+theme(plot.title = element_text(hjust = 0.5))


#dev.off()

#####
###############################################
###############################################
############ GO analyses ######################
###############################################

# shinygo http://bioinformatics.sdstate.edu/go/
library(pals)
library(RColorBrewer)
library(plotfunctions)
library(circlize)

high_go <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/GO/high_level_GO_ai_all.csv", stringsAsFactors = FALSE, header = TRUE, sep = ",")
head(high_go$High.level.GO.category)
total_genes <- sum(high_go$N)
slices <- c(high_go[which(high_go$High.level.GO.category == "Anatomical structure development "),][1,1],
            high_go[which(high_go$High.level.GO.category == "Biosynthetic process "),][1,1],
            high_go[which(high_go$High.level.GO.category == "Regulation of metabolic process "),][1,1],
            high_go[which(high_go$High.level.GO.category == "Cellular component biogenesis "),][1,1],
            high_go[which(high_go$High.level.GO.category == "Signaling "),][1,1])
            #high_go[which(high_go$High.level.GO.category == "Developmental growth "),][1,1],
            #high_go[which(high_go$High.level.GO.category == "Pigmentation "),][1,1],
            #high_go[which(high_go$High.level.GO.category == "Behavior "),][1,1],
            #high_go[which(high_go$High.level.GO.category == "Reproduction "),][1,1]
other_terms <- total_genes - sum(slices)
slices <- c(slices,other_terms)
slices_per <- paste(round((slices/other_terms) *100,digits = 0),"%", sep = "")
slices_per[6] <- ""
pie_cols <- c(red,blu,yel,gre,pur,"gray56")
tiff("C:/Users/jmcgirr/Documents/all_2018_samples/manuscript_figs/high_level_go_ai_all.tiff", width = 4.5, height = 4.5, units = 'in', res = 1000)
pie(slices, col = pie_cols, labels = slices_per,init.angle=90, cex=0.8)
dev.off()

high_go <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/GO/high_level_GO_pmp_all.csv", stringsAsFactors = FALSE, header = TRUE, sep = ",")
total_genes <- sum(high_go$N)
slices <- c(high_go[which(high_go$High.level.GO.category == "Anatomical structure development "),][1,1],
            high_go[which(high_go$High.level.GO.category == "Biosynthetic process "),][1,1],
            high_go[which(high_go$High.level.GO.category == "Regulation of metabolic process "),][1,1],
            high_go[which(high_go$High.level.GO.category == "Cellular component biogenesis "),][1,1],
            high_go[which(high_go$High.level.GO.category == "Signaling "),][1,1])
            #high_go[which(high_go$High.level.GO.category == "Developmental growth "),][1,1],
            #high_go[which(high_go$High.level.GO.category == "Pigmentation "),][1,1])
other_terms <- total_genes - sum(slices)
slices <- c(slices,other_terms)
slices_per <- paste(round((slices/other_terms) *100,digits = 0),"%", sep = "")
slices_per[6] <- ""
pie_cols <- c(red,blu,yel,gre,pur,"gray56")
tiff("C:/Users/jmcgirr/Documents/all_2018_samples/manuscript_figs/high_level_go_pmp_all.tiff", width = 4.5, height = 4.5, units = 'in', res = 1000)
pie(rev(slices), col = rev(pie_cols), labels = rev(slices_per),init.angle=90, cex=0.8)
dev.off()


high_go <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/GO/high_level_GO_ai_fst_dxy.csv", stringsAsFactors = FALSE, header = TRUE, sep = ",")
total_genes <- sum(high_go$N)
slices <- c(high_go[which(high_go$High.level.GO.category == "Anatomical structure development "),][1,1],
            high_go[which(high_go$High.level.GO.category == "Biosynthetic process "),][1,1],
            high_go[which(high_go$High.level.GO.category == "Regulation of metabolic process "),][1,1],
            high_go[which(high_go$High.level.GO.category == "Cellular component biogenesis "),][1,1]
            )
other_terms <- total_genes - sum(slices)
slices <- c(slices,other_terms)
slices_per <- paste(round((slices/other_terms) *100,digits = 0),"%", sep = "")
slices_per[5] <- ""
pie_cols <- c(red,blu,yel,gre,"gray56")
tiff("C:/Users/jmcgirr/Documents/all_2018_samples/manuscript_figs/high_level_go_ai_fst_dxy.tiff", width = 4.5, height = 4.5, units = 'in', res = 1000)
pie(rev(slices), col = rev(pie_cols), labels = rev(slices_per),init.angle=270, cex=0.8)
dev.off()

high_go <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/GO/high_level_GO_pmp_fst_dxy.csv", stringsAsFactors = FALSE, header = TRUE, sep = ",")
total_genes <- sum(high_go$N)
slices <- c(high_go[which(high_go$High.level.GO.category == "Anatomical structure development "),][1,1],
            high_go[which(high_go$High.level.GO.category == "Biosynthetic process "),][1,1],
            high_go[which(high_go$High.level.GO.category == "Regulation of metabolic process "),][1,1],
            high_go[which(high_go$High.level.GO.category == "Cellular component biogenesis "),][1,1],
            high_go[which(high_go$High.level.GO.category == "Signaling "),][1,1])
other_terms <- total_genes - sum(slices)
slices <- c(slices,other_terms)
slices_per <- paste(round((slices/other_terms) *100,digits = 0),"%", sep = "")
slices_per[6] <- ""
pie_cols <- c(red,blu,yel,gre,pur,"gray56")
tiff("C:/Users/jmcgirr/Documents/all_2018_samples/manuscript_figs/high_level_go_pmp_fst_dxy.tiff", width = 4.5, height = 4.5, units = 'in', res = 1000)
pie(slices, col = pie_cols, labels = slices_per,init.angle=270, cex=0.8)
dev.off()









go <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/GO/GO_enriched_ai_fst_dxy.csv", stringsAsFactors = FALSE, header = TRUE, sep = ",")
head(go)
go <- go[order(go$Genes.in.list, decreasing = TRUE),]


pal = colorRampPalette(c(red, blu))
go$order = findInterval(go$Enrichment.FDR, sort(go$Enrichment.FDR))

tiff("C:/Users/jmcgirr/Documents/all_2018_samples/manuscript_figs/enriched_ai.tiff", width = 5, height = 5, units = 'in', res = 1000)
par(las=2) # make label text perpendicular to axis
par(mar=c(5,15,4,2)) # increase y-axis margin.
barplot(go$Genes.in.list, col=pal(nrow(go))[go$order],horiz=TRUE, names.arg=c(tolower(go$Functional.Category)), xaxt = 'n',cex.names=0.6)
axis(1,c(0,2,4,6,8,10,12), las=1,cex.axis=0.7)
gradientLegend(valRange=range(go$Enrichment.FDR), color = pal(nrow(go)),pos = c(9,15,10,20), n.seg=2,dec = 3, coords = TRUE)
dev.off()

go <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/GO/GO_enriched_pmp_fst_dxy.csv", stringsAsFactors = FALSE, header = TRUE, sep = ",")
head(go)
go <- go[order(go$Genes.in.list, decreasing = TRUE),]
go$order = findInterval(go$Enrichment.FDR, sort(go$Enrichment.FDR))

tiff("C:/Users/jmcgirr/Documents/all_2018_samples/manuscript_figs/enriched_pmp.tiff", width = 5, height = 5, units = 'in', res = 1000)
par(las=2) # make label text perpendicular to axis
par(mar=c(5,15,4,2)) # increase y-axis margin.
barplot(go$Genes.in.list, col=pal(nrow(go))[go$order],horiz=TRUE, names.arg=c(tolower(go$Functional.Category)), xaxt = 'n',cex.names=0.58)
axis(1,c(0,2,4,6,8,10,12,14), las=1,cex.axis=0.7)
gradientLegend(valRange=range(go$Enrichment.FDR), color = pal(nrow(go)),pos = c(9,15,10,20), n.seg=2,dec = 3, coords = TRUE)
dev.off()







labs = c('alphabet', 'alphabet2', 'glasbey', 'kelly', 'polychrome')
op = par(mar = c(0, 5, 3, 1))
pal.bands(alphabet(), alphabet2(), glasbey(), kelly(), polychrome(), 
          labels = labs, show.names = FALSE)

test <- as.matrix(read.csv("C:/Users/jmcgirr/Documents/all_2018_samples/GO/GO_enriched_adaptive_misexpression_fst_dxy_sec24d.csv" ,header = TRUE,row.names=1))
grid.col <- c(c(1:26),red,blu,yel,gre,pur,grb,lir)
names(grid.col) = c(rownames(test), colnames(test))
grid.col
chordDiagram(test,annotationTrack = c("name", "grid"), grid.col = grid.col)


test <- as.matrix(read.csv("C:/Users/jmcgirr/Documents/all_2018_samples/GO/GO_enriched_parallel_mp_fst_dxy_itga8.csv" ,header = TRUE,row.names=1))
grid.col <- c(mycols,"black","grey","black","grey","black","grey","black","grey","black","grey")
names(grid.col) = c(rownames(test), colnames(test))
grid.col
chordDiagram(test,annotationTrack = c("name", "grid"), grid.col = grid.col)
mycols <- c("#0000FF","#FF0000","#00FF00","#FFA500","#FF00B6","#005300","#FFD300","#009FFF","#9A4D42","#00FFBE","#783FC1","#1F9698","#FFACFD","#B1CC71","#F1085C","#FE8F42","#222222","#201A01","#720055","#766C95","#02AD24","#C8FF00","#886C00","#FFB79F","#858567","#A10300","#14F9FF","#00479E","#DC5E93","#93D4FF","#004CFF","#006400")
kelly()


#####
###############################################
###############################################
############ correct dxy ######################
###############################################
crosses <- c("caxcm","caxcp","oaxom", "oaxop", "cmxcp", "opxom")

i <- 1
for (i in c(1:length(crosses)))
{
Bmbodxy <- read.csv(paste("D:/Martin Lab/rna_2018/fst/dna/",crosses[i],"_popgen_dna_stats.csv",sep=""), header = TRUE,stringsAsFactors = FALSE)
head(Bmbodxy)

#correct avg dxy in a window given that it should be a weighted average since we dont factor in the non-variant sites where dxy is 0

Bmbodxy$total_length<-(Bmbodxy$end-Bmbodxy$start)
Bmbodxy$dxy_sum_variants<-(Bmbodxy[,8]*Bmbodxy$sites) #weighting the dxy given the number of sites, for non variant sites would be something like (0*#sites_skipped)
Bmbodxy$corr_dxy<-(Bmbodxy$dxy_sum_variants/Bmbodxy$total_length) #+0/Bmbodxy$total_length, but that doesnt actually matter
Bmbodxy$pi_a_sum_variants<-(Bmbodxy[,6]*Bmbodxy$sites) #weighting the dxy given the number of sites, for non variant sites would be something like (0*#sites_skipped)
Bmbodxy$corr_pi_a<-(Bmbodxy$pi_a_sum_variants/Bmbodxy$total_length) #+0/Bmbodxy$total_length, but that doesnt actually matter
Bmbodxy$pi_b_sum_variants<-(Bmbodxy[,7]*Bmbodxy$sites) #weighting the dxy given the number of sites, for non variant sites would be something like (0*#sites_skipped)
Bmbodxy$corr_pi_b<-(Bmbodxy$pi_b_sum_variants/Bmbodxy$total_length) #+0/Bmbodxy$total_length, but that doesnt actually matter

write.csv(Bmbodxy,paste("D:/Martin Lab/rna_2018/fst/dna/",crosses[i],"_popgen_dna_stats_corr_dxy.csv",sep=""), quote = FALSE, row.names = FALSE)
}


#####
###############################################
###############################################
############ top sweed hits ###################
###############################################
ca_sweed <- read.csv("C:/Users/jmcgirr/Documents/all_2018_samples/interesting_genes_tables/ca_sweed_95.csv", header = TRUE, row.names = NULL, stringsAsFactors = FALSE)
cm_sweed <- read.csv("C:/Users/jmcgirr/Documents/all_2018_samples/interesting_genes_tables/cm_sweed_95.csv", header = TRUE, row.names = NULL, stringsAsFactors = FALSE)
cp_sweed <- read.csv("C:/Users/jmcgirr/Documents/all_2018_samples/interesting_genes_tables/cp_sweed_95.csv", header = TRUE, row.names = NULL, stringsAsFactors = FALSE)
oa_sweed <- read.csv("C:/Users/jmcgirr/Documents/all_2018_samples/interesting_genes_tables/oa_sweed_95.csv", header = TRUE, row.names = NULL, stringsAsFactors = FALSE)
om_sweed <- read.csv("C:/Users/jmcgirr/Documents/all_2018_samples/interesting_genes_tables/om_sweed_95.csv", header = TRUE, row.names = NULL, stringsAsFactors = FALSE)
op_sweed <- read.csv("C:/Users/jmcgirr/Documents/all_2018_samples/interesting_genes_tables/op_sweed_95.csv", header = TRUE, row.names = NULL, stringsAsFactors = FALSE)
ca_sweed <- unique(ca_sweed$related_accession)
cm_sweed <- unique(cm_sweed$related_accession)
cp_sweed <- unique(cp_sweed$related_accession)
oa_sweed <- unique(oa_sweed$related_accession)
om_sweed <- unique(om_sweed$related_accession)
op_sweed <- unique(op_sweed$related_accession)
caxcm_sweed <- unique(c(ca_sweed,cm_sweed))
caxcp_sweed<-  unique(c(ca_sweed,cp_sweed))
cmxcp_sweed<-  unique(c(cm_sweed,cp_sweed))
oaxom_sweed<-  unique(c(oa_sweed,om_sweed))
oaxop_sweed<-  unique(c(oa_sweed,op_sweed))
omxop_sweed<-  unique(c(om_sweed,op_sweed))

head(cm_sweed)
intersect(oaxom_8_ai , om_sweed)
intersect(caxcp_8_ai , cp_sweed)
intersect(oaxop_48_ai, op_sweed)
intersect(oaxop_8_ai , op_sweed)
intersect(cmxcp_48_ai, cmxcp_sweed)
intersect(cmxcp_8_ai , cmxcp_sweed)
intersect(omxop_48_ai, omxop_sweed)
intersect(omxop_8_ai , omxop_sweed)


intersect(intersect(oaxom_8_ai ,om_sweed),intersect(oaxom_dxy,oaxom_fixed))
intersect(intersect(caxcp_8_ai ,cm_sweed),intersect(caxcp_dxy,caxcp_fixed))
intersect(intersect(oaxop_48_ai,op_sweed),intersect(oaxop_dxy,oaxop_fixed))
intersect(intersect(oaxop_8_ai ,op_sweed),intersect(oaxop_dxy,oaxop_fixed))
intersect(intersect(cmxcp_48_ai,cmxcp_sweed),intersect(cmxcp_dxy,cmxcp_fixed))
intersect(intersect(cmxcp_8_ai ,cmxcp_sweed),intersect(cmxcp_dxy,cmxcp_fixed))
intersect(intersect(omxop_48_ai,omxop_sweed),intersect(omxop_dxy,omxop_fixed))
intersect(intersect(omxop_8_ai ,omxop_sweed),intersect(omxop_dxy,omxop_fixed))



intersect(oaxom_8_ai , om_taj)
intersect(caxcp_8_ai , cp_taj)
intersect(oaxop_48_ai, op_taj)
intersect(oaxop_8_ai , om_taj)
intersect(cmxcp_48_ai, cmxcp_sweed)
intersect(cmxcp_8_ai , cmxcp_taj)
intersect(omxop_48_ai, omxop_sweed)
intersect(omxop_8_ai , omxop_sweed)






