######## DEseq Pipeline ####### 
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
#install.packages("rlang")
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



#####
##### make a comparison table for a single comparison ####
##########################################################

stage            <- "48hpf"

pops1_1 <- c("CRPA")
pops1_2 <- c("OSPA")
pops1_3 <- c()
pops1_4 <- c()
pops1_5 <- c()
pops1_6 <- c()

pops2_1 <- c("CRPP")
pops2_2 <- c("OSPP")
pops2_3 <- c()
pops2_4 <- c()
pops2_5 <- c()
pops2_6 <- c()

pops1 <- c(pops1_1,	pops1_2,	pops1_3,	pops1_4,	pops1_5,	pops1_6)
pops2 <- c(pops2_1,	pops2_2,	pops2_3,	pops2_4,	pops2_5,	pops2_6)
pops1_name <- paste(pops1,collapse="_&_")
pops2_name <- paste(pops2,collapse="_&_")
pops1_name <- paste("(", pops1_name, ")", sep = "")
pops2_name <- paste("(", pops2_name, ")", sep = "")
single_run <- paste(pops1_name, "_vs_", pops2_name, "_" ,stage, sep = "")
comp_name <- paste("condition_species_",pops1_name, "_vs_", pops2_name, "_" ,stage,".txt", sep = "")

setwd("D:/Martin Lab/rna_2018/all_2018_samples/conditions/")
master <- read.table("D:/Martin Lab/rna_2018/all_2018_samples/table_maker_master.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
head(master)
pops1_table <- master[master$f1 %in% pops1, ]
pops2_table <- master[master$f1 %in% pops2, ]
pops1_table$species <- "a" 
pops2_table$species <- "b"
comp_table <- rbind(pops1_table,pops2_table)
comp_table <- comp_table[which(comp_table$stage == stage),]
comp_table
sample_size_a <- nrow(comp_table[which(comp_table$species == "a"),])
sample_size_b <- nrow(comp_table[which(comp_table$species == "b"),])
sample_size_a
sample_size_b

write.table(comp_table, comp_name, row.names = FALSE, quote= FALSE,sep="\t")
#comp_table <- comp_table[which(comp_table$sample != "OAE5"),]
#write.table(comp_table, "condition_species_ospA_vs_ospP_8dpf_no_OAE5.txt", row.names = FALSE, quote= FALSE,sep="\t")
#general split
#comp_table <- master[which(master$sequencing_round == 1 & master$stage == "8dpf" & master$cross_type == "p"),]
#comp_name <- "round_1_pure_cross_8dpf.txt"


#####
###############################################################################################
###############################################################################################
##### perform DE analyses for comparisons. Output MA plots, gene tables, and summary stats ####
###############################################################################################
###############################################################################################

#### Use comparisons.txt to loop all comparisons ###
#### To make single comparison, uncomment next line and skip to single_comp_only below ###

# compare hybrids
all_comps <- read.table("D:/Martin Lab/rna_2018/fst/hybrid_mse_comparisons.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
# compare all
#all_comps <- read.table("D:/Martin Lab/rna_2018/all_2018_samples/comparisons.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
#all_comps <- read.table("D:/Martin Lab/rna_2018/all_2018_samples/bootstrapping_and_downsampling/boot_and_downs_comparisons.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
head(all_comps)
final_features <- read.table("D:/Martin Lab/rna_2018/all_2018_samples/features_gff.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t", quote = "")
all_cts <- read.table("D:/Martin Lab/rna_2018/all_2018_samples/conditions/all_48hpf_and_8dpf_counts_my_gtf_geneid_final", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
blast_key <- read.table("D:/Martin Lab/rna_2018/all_2018_samples/blast/cyprinodon_to_danio_one_way_best_hit_symbols.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
head(blast_key)
#i <- 84
#subset_comps <- c(91:98)

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

for (i in (1:nrow(all_comps)))
#for (i in subset_comps)
{
  pops1_name <- all_comps$pops1[i]
  pops2_name <- all_comps$pops2[i]
  pops1_t <- pops1_name
  pops2_t <- pops2_name
  #pops1 <- strsplit(pops1_name,"_&_")[[1]]
  #pops2 <- strsplit(pops2_name,"_&_")[[1]]
  pops1 <- strsplit(pops1_name,"_and_")[[1]]
  pops2 <- strsplit(pops2_name,"_and_")[[1]]
  stage <- all_comps$stage[i]
  #pops1_name <- paste("(", pops1_name, ")", sep = "")
  #pops2_name <- paste("(", pops2_name, ")", sep = "")
  comp_name <- paste("condition_species_",pops1_name, "_vs_", pops2_name, "_" ,stage,".txt", sep = "")
  single_run <- paste(pops1_name, "_vs_", pops2_name, "_" ,stage, sep = "")
  
  setwd("D:/Martin Lab/rna_2018/all_2018_samples/conditions/")
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
  
  setwd("D:/Martin Lab/rna_2018/all_2018_samples/conditions/")
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
  de_plot <- paste("D:/Martin Lab/rna_2018/all_2018_samples/de_plots/boots_and_dns/",comp_file, "_de_plot.tiff", sep = "")
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
write.table(sum_stats, "D:/Martin Lab/rna_2018/plug_n_chug//summary_stats_subsetgroups.txt", row.names = FALSE, quote= FALSE,sep="\t")


#####
########################################
########################################
##### inheritance patterns #############
########################################

#### make plot like McManus 2010 4B

inheritance_comps <- read.table("D:/Martin Lab/rna_2018/all_2018_samples/inheritance_patterns/inheritance_comparisons.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
DE_genes_dir <- "D:/Martin Lab/rna_2018/all_2018_samples/conditions/"
inheritance_plots_dir <- "D:/Martin Lab/rna_2018/all_2018_samples/inheritance_patterns/plots/"
#i <- 5
subset_comps <- c(15:18)

for (i in subset_comps)
#for (i in (1:nrow(inheritance_comps)))
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
  
  cols <- 	c("#f6c700", "#E77200", "#FF00CC","#000000", "#bd1e24", "#0067a7")
  
  #con = yellow ="#f6c700"
  #add = orange = "#E77200"
  #p1_dom = pink = "#FF00CC"
  #p2_dom = black = "#000000"
  #over_dom = red = #bd1e24
  #under_dom = blue = #0067a7
  
  cols_in <- cols[inheritance$inheritance]
  x_name <- paste(parent_pops1_name, hy_name, sep = " vs. ")
  y_name <- paste(parent_pops2_name, hy_name, sep = " vs. ")
  
  tiff(outfile_plot, width = 7, height = 7, units = 'in', res = 1000)
  plot(inheritance$lfc_p1h, inheritance$lfc_p2h, pch =16, col = cols_in, 
       ylab = paste("log2 fold change", y_name),
       xlab = paste("log2 fold change", x_name),cex.axis=1.5,
       main = plot_title )
  abline(v=0, col = "black", lty = 3, lwd = 1.8)
  abline(h=0, col = "black", lty = 3, lwd = 1.8)
  legend("topleft", inset=0,
         c(paste("conserved", prop_con),
           paste("additive", prop_add),
           paste(parent_pops1_name, "dominant", prop_p1_dom),
           paste(parent_pops2_name, "dominant", prop_p2_dom), 
           paste("overdominant", prop_overdom), 
           paste("underdominant", prop_underdom)),
         cex = 0.9, fill=cols, horiz=FALSE, bty="n")
  dev.off()
  
  
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


setwd("D:/Martin Lab/rna_2018/all_2018_samples/")
#comp <- "all_samples_8dpf"
comp <- "all_samples_48hpf"
comp <- "table_maker_master"
comp_file <- paste(comp , ".txt", sep = "")
#comp_file <- "table_maker_master_outlier_rm.txt"
all_cts <- read.table("D:/Martin Lab/rna_2018/all_2018_samples/conditions/all_48hpf_and_8dpf_counts_my_gtf_geneid_final", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
sample_list <- read.table(comp_file, header = TRUE, stringsAsFactors = FALSE)
keeps <- c("Geneid", sample_list$sample)
keeper <- sample_list$sample
cts <- all_cts[keeps]
cts_out <- paste(comp, "_counts.txt", sep = "")
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
                                design= ~f1+sequencing_round)
#}else 
#{
#dds <- DESeqDataSetFromMatrix(countData = cts,
#                              colData = colData,
#                              design= ~f1)
#}

dds <- estimateSizeFactors(dds)
idx <- rowSums(counts(dds, normalized=TRUE) >= 1 ) >= (nrow(colData)/2)
length(rowSums(counts(dds, normalized=TRUE) >= 1 ) >= (nrow(colData)/2))
length(idx)
dds <- dds[idx,]
dds <- DESeq(dds)

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
plotPCA(rld.sub, intgroup=c("f1"))+geom_point(aes(shape='parents'))+scale_shape_manual(values=c(25))

pcaData <- plotPCA(rld.sub, intgroup=c("f1",'sequencing_round'), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=f1, shape=sequencing_round)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
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


##### WGCNA 

#source("http://bioconductor.org/biocLite.R")
#biocLite("WGCNA")
#install.packages("flashClust")
library(WGCNA)
library(flashClust)

setwd("D:/Martin Lab/rna_2018/conditions/")
cts <- read.table("counts_round_1_rna_2018_final", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
sample_list <- read.table("round_1_pure_cross_8dpf.txt", header = TRUE, stringsAsFactors = FALSE)
keeps <- c("Geneid", sample_list$sample)
keeper <- sample_list$sample
cts <- cts[keeps]
#write.table(cts, "round_1_8dpf_counts.txt", row.names = FALSE, quote= FALSE,sep="\t")
cts <- as.matrix(read.table("round_1_8dpf_counts.txt" ,sep = "\t",header = TRUE,row.names=1))

dim(cts)
datExpr = as.data.frame(t(cts))
head(datExpr)
dim(datExpr)

# Run this to check if there are gene outliers
gsg = goodSamplesGenes(datExpr, verbose = 3)
gsg$allOK


#If the last statement returns TRUE, all genes have passed the cuts. If not, we remove the offending genes and samples from the data with the following:
if (!gsg$allOK)
{if (sum(!gsg$goodGenes)>0)
  printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse= ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse=", ")))
  datExpr= datExpr[gsg$goodSamples, gsg$goodGenes]
}



#build a adjacency "correlation" matrix
enableWGCNAThreads()
softPower = 18
adjacency = adjacency(datExpr, power = softPower, type = "signed") #specify network type
head(adjacency)
TOM = TOMsimilarity(adjacency, TOMType="signed") # specify network type
dissTOM = 1-TOM




#####
#####################################################
#####################################################
#### filter genes DE between outgroups and lakes ####
#####################################################
#####################################################

## hybrid filter compraisons
DE_genes_dir <- "D:/Martin Lab/rna_2018/all_2018_samples/conditions/"
stage <- "8dpf"


  ssi_a_nca_vs_hybrids <- read.csv(paste(DE_genes_dir,"DE_(CRPA)_vs_(NAxCA)_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  ssi_a_cun_vs_hybrids <- read.csv(paste(DE_genes_dir,"DE_(CRPA)_vs_(UPxCA)_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  pure_a_vs_lake_cross <- read.csv(paste(DE_genes_dir,"DE_(CRPA_&_OSPA)_vs_(CAxOA_&_OAxCA)_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  ssi_a_nca_vs_hybrids <- ssi_a_nca_vs_hybrids[which(ssi_a_nca_vs_hybrids$padj < 0.05),]
  ssi_a_cun_vs_hybrids <- ssi_a_cun_vs_hybrids[which(ssi_a_cun_vs_hybrids$padj < 0.05),]
  pure_a_vs_lake_cross <- pure_a_vs_lake_cross[which(pure_a_vs_lake_cross$padj < 0.05),]
  ssi_a_nca_vs_hybrids$zeb_tag <- paste(ssi_a_nca_vs_hybrids$tag, ssi_a_nca_vs_hybrids$zeb_gene_symbol, sep = ";")
  ssi_a_cun_vs_hybrids$zeb_tag <- paste(ssi_a_cun_vs_hybrids$tag, ssi_a_cun_vs_hybrids$zeb_gene_symbol, sep = ";")
  pure_a_vs_lake_cross$zeb_tag <- paste(pure_a_vs_lake_cross$tag, pure_a_vs_lake_cross$zeb_gene_symbol, sep = ";")
  ssi_a_nca_vs_hybrids <- ssi_a_nca_vs_hybrids$zeb_tag
  ssi_a_cun_vs_hybrids <- ssi_a_cun_vs_hybrids$zeb_tag
  pure_a_vs_lake_cross <- pure_a_vs_lake_cross$zeb_tag
  
  # snail
  am_vs_hybrids   <- read.csv(paste(DE_genes_dir,"DE_(CRPA_&_OSPA_&_CRPM_&_OSPM)_vs_(CAxCM_&_CMxCA_&_OAxOM_&_OMxOA)_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  m_vs_lake_cross <- read.csv(paste(DE_genes_dir,"DE_(CRPM_&_OSPM)_vs_(CMxOM)_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  am_vs_hybrids   <- am_vs_hybrids[which(am_vs_hybrids$padj < 0.05),]
  m_vs_lake_cross <- m_vs_lake_cross[which(m_vs_lake_cross$padj < 0.05),]
  am_vs_hybrids$zeb_tag   <- paste(am_vs_hybrids$tag, am_vs_hybrids$zeb_gene_symbol, sep = ";")
  m_vs_lake_cross$zeb_tag <- paste(m_vs_lake_cross$tag, m_vs_lake_cross$zeb_gene_symbol, sep = ";")
  
  am_vs_hybrids   <- am_vs_hybrids$zeb_tag
  m_vs_lake_cross <- m_vs_lake_cross$zeb_tag
  
  misexpressed_axm_hybrids_only <- setdiff(am_vs_hybrids,(unique(c(m_vs_lake_cross,ssi_a_nca_vs_hybrids,ssi_a_cun_vs_hybrids,pure_a_vs_lake_cross))))
  length(misexpressed_axm_hybrids_only)
  #misexpressed_axm_hybrids_only <- setdiff(ssi_a_nca_vs_hybrids,ssi_a_cun_vs_hybrids)
  #length(misexpressed_axm_hybrids_only)
  #length(ssi_a_nca_vs_hybrids)
  #length(unique(ssi_a_cun_vs_hybrids))
  #length(ssi_a_cun_vs_hybrids)
  
  #install.packages('gplots')
  #require(gplots)
  #venn(list(GrpA=ssi_a_nca_vs_hybrids,GrpB=ssi_a_cun_vs_hybrids,GrpC=pure_a_vs_lake_cross,GrpD=m_vs_lake_cross, GrpE=am_vs_hybrids))
  #install.packages('venn')
  library(venn)
  
  # scale 
  ap_vs_hybrids   <- read.csv(paste(DE_genes_dir,"DE_(CRPA_&_OSPA_&_CRPP_&_OSPP)_vs_(CAxCP_&_OAxOP_&_OPxOA)_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  p_vs_lake_cross <- read.csv(paste(DE_genes_dir,"DE_(CRPP_&_OSPP)_vs_(OPxCP)_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  ap_vs_hybrids   <- ap_vs_hybrids[which(ap_vs_hybrids$padj < 0.05),]
  p_vs_lake_cross <- p_vs_lake_cross[which(p_vs_lake_cross$padj < 0.05),]
  
  ap_vs_hybrids$zeb_tag   <- paste(ap_vs_hybrids$tag, ap_vs_hybrids$zeb_gene_symbol, sep = ";")
  p_vs_lake_cross$zeb_tag <- paste(p_vs_lake_cross$tag, p_vs_lake_cross$zeb_gene_symbol, sep = ";")
  
  ap_vs_hybrids   <- ap_vs_hybrids$zeb_tag
  p_vs_lake_cross <- p_vs_lake_cross$zeb_tag
  
  misexpressed_axp_hybrids_only <- setdiff(ap_vs_hybrids,(unique(c(p_vs_lake_cross,ssi_a_nca_vs_hybrids,ssi_a_cun_vs_hybrids,pure_a_vs_lake_cross))))
  length(misexpressed_axp_hybrids_only)
  length(intersect(misexpressed_axp_hybrids_only,ssi_a_vs_p ))
  length(intersect(misexpressed_axm_hybrids_only,ssi_a_vs_m ))
  
  #tiff(paste("D:/Martin Lab/rna_2018/all_2018_samples/filter_comparisons/hybrid_cross_filter/axm_hybrid_", stage, "_venn.tiff", sep = ""), width = 7, height = 7, units = 'in', res = 1000)
  venn(list(NCA_x_CRPA=ssi_a_nca_vs_hybrids,CUN_x_CRPA=ssi_a_cun_vs_hybrids,A_lake_cross=pure_a_vs_lake_cross,M_lake_cross=m_vs_lake_cross, A_x_M_hybrids=am_vs_hybrids), 
       ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
  #dev.off()
  
  #tiff(paste("D:/Martin Lab/rna_2018/all_2018_samples/filter_comparisons/hybrid_cross_filter/axp_hybrid_", stage, "_venn.tiff", sep = ""), width = 7, height = 7, units = 'in', res = 1000)
  venn(list(NCA_x_CRPA=ssi_a_nca_vs_hybrids,CUN_x_CRPA=ssi_a_cun_vs_hybrids,A_lake_cross=pure_a_vs_lake_cross,P_lake_cross=p_vs_lake_cross, A_x_P_hybrids=ap_vs_hybrids), 
       ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
  #dev.off()
  
  #write.csv(misexpressed_axm_hybrids_only, paste("D:/Martin Lab/rna_2018/all_2018_samples/filter_comparisons/hybrid_cross_filter/misexpressed_axm_hybrids_only_", stage,".csv", sep = ""), row.names = FALSE, quote = FALSE)
  #write.csv(misexpressed_axp_hybrids_only, paste("D:/Martin Lab/rna_2018/all_2018_samples/filter_comparisons/hybrid_cross_filter/misexpressed_axp_hybrids_only_", stage,".csv", sep = ""), row.names = FALSE, quote = FALSE)
  
  #biocLite("eulerr")
  #library(eulerr)
  #install_github("jolars/eulerr")
  #cols <- 	c("#f6c700", "#E77200", "#FF00CC","#000000", "#bd1e24", "#0067a7")
  #VennDiag <- euler(list(NCA_x_CRPA=ssi_a_nca_vs_hybrids,CUN_x_CRPA=ssi_a_cun_vs_hybrids,A_lake_cross=pure_a_vs_lake_cross,P_lake_cross=p_vs_lake_cross, A_x_P_hybrids=ap_vs_hybrids)) 
  #plot(VennDiag, counts = TRUE, font=1, cex=0.5, alpha=0.6,
  #     fill=cols,quantities = list(fontsize = 8),
  #     labels=c("","","","",""))
  #VennDiag <- euler(list(NCA_x_CRPA=unique(ssi_a_nca_vs_hybrids),CUN_x_CRPA=unique(ssi_a_cun_vs_hybrids),A_lake_cross=unique(pure_a_vs_lake_cross),M_lake_cross=unique(m_vs_lake_cross), A_x_M_hybrids=unique(am_vs_hybrids))) 
  #plot(VennDiag, counts = FALSE, font=1, cex=0.1, alpha=0.6,
  #     fill=cols,quantities = list(fontsize = 8),
  #     labels=c("","","","",""))
  #library(venneuler)
  #vd <- venneuler(NCA_x_CRPA=ssi_a_nca_vs_hybrids,CUN_x_CRPA=ssi_a_cun_vs_hybrids,A_lake_cross=pure_a_vs_lake_cross,P_lake_cross=p_vs_lake_cross, A_x_P_hybrids=ap_vs_hybrids)
  #plot(vd)
  #
  #
  #biocLite(c("RBGL","graph"))
  #library(vennerable)
  #install.packages("devtools")
  #library(devtools)
  #install_github("js229/Vennerable")
  #library(Vennerable)
  #vignette("Venn")
  #
  #Vcombo <- Venn(SetNames = c("1","2","3","4","5"),Weight = list(NCA_x_CRPA=ssi_a_nca_vs_hybrids,CUN_x_CRPA=ssi_a_cun_vs_hybrids,A_lake_cross=pure_a_vs_lake_cross,P_lake_cross=p_vs_lake_cross, A_x_P_hybrids=ap_vs_hybrids))
  #plot(Vcombo)
  #
  #install.packages("rJava")
  #library(rJava)
  #Sys.setenv(JAVA_HOME="C:/Program Files/Java/jdk-10.0.1/")



## pure population filter comparisons
DE_genes_dir <- "D:/Martin Lab/rna_2018/all_2018_samples/conditions/"
stage <- "8dpf"

  ssi_a_vs_nca <- read.csv(paste(DE_genes_dir,"DE_(CRPA_&_OSPA)_vs_(NCA)_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  ssi_a_vs_cun <- read.csv(paste(DE_genes_dir,"DE_(CRPA_&_OSPA)_vs_(UPxUA)_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  a_lake_cross <- read.csv(paste(DE_genes_dir,"DE_(CRPA)_vs_(OSPA)_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  ssi_a_vs_m   <- read.csv(paste(DE_genes_dir,"DE_(CRPA_&_OSPA)_vs_(CRPM_&_OSPM)_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  ssi_a_vs_p   <- read.csv(paste(DE_genes_dir,"DE_(CRPA_&_OSPA)_vs_(CRPP_&_OSPP)_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  ssi_a_vs_nca  <- ssi_a_vs_nca [which(ssi_a_vs_nca$padj < 0.05),]
  ssi_a_vs_cun  <- ssi_a_vs_cun [which(ssi_a_vs_cun$padj < 0.05),]
  a_lake_cross  <- a_lake_cross [which(a_lake_cross$padj < 0.05),]
  ssi_a_vs_m    <- ssi_a_vs_m   [which(ssi_a_vs_m$padj < 0.05),]
  ssi_a_vs_p    <- ssi_a_vs_p   [which(ssi_a_vs_p$padj < 0.05),]
  ssi_a_vs_nca$zeb_tag  <- paste(ssi_a_vs_nca$tag, ssi_a_vs_nca$zeb_gene_symbol, sep = ";")  
  ssi_a_vs_cun$zeb_tag  <- paste(ssi_a_vs_cun$tag, ssi_a_vs_cun$zeb_gene_symbol, sep = ";")
  a_lake_cross$zeb_tag  <- paste(a_lake_cross$tag, a_lake_cross$zeb_gene_symbol, sep = ";")
  ssi_a_vs_m$zeb_tag  <-   paste(ssi_a_vs_m$tag, ssi_a_vs_m$zeb_gene_symbol, sep = ";")
  ssi_a_vs_p$zeb_tag  <-   paste(ssi_a_vs_p$tag, ssi_a_vs_p$zeb_gene_symbol, sep = ";")
  
  ssi_a_vs_nca  <- ssi_a_vs_nca$zeb_tag
  ssi_a_vs_cun  <- ssi_a_vs_cun$zeb_tag
  a_lake_cross  <- a_lake_cross$zeb_tag
  ssi_a_vs_m    <- ssi_a_vs_m$zeb_tag
  ssi_a_vs_p    <- ssi_a_vs_p$zeb_tag
  length(ssi_a_vs_p)
  
  de_am_only <- setdiff(ssi_a_vs_m,(unique(c(a_lake_cross,ssi_a_vs_cun,ssi_a_vs_nca))))
  de_ap_only <- setdiff(ssi_a_vs_p,(unique(c(a_lake_cross,ssi_a_vs_cun,ssi_a_vs_nca))))
  length(de_am_only)
  length(de_ap_only)
  
  #write.table(de_am_only, paste("D:/Martin Lab/rna_2018/all_2018_samples/filter_comparisons/pure_cross_filter/de_am_only_", stage,".csv", sep = ""), row.names = FALSE, quote = FALSE, col.names = FALSE)
  #write.table(de_ap_only, paste("D:/Martin Lab/rna_2018/all_2018_samples/filter_comparisons/pure_cross_filter/de_ap_only_", stage,".csv", sep = ""), row.names = FALSE, quote = FALSE, col.names = FALSE)
  
  #library(venn)
  #tiff(paste("D:/Martin Lab/rna_2018/all_2018_samples/filter_comparisons/pure_cross_filter/a_vs_m_", stage, "_venn.tiff", sep = ""), width = 7, height = 7, units = 'in', res = 1000)
  venn(list(A_vs_NCA=ssi_a_vs_nca,A_vs_CUN=ssi_a_vs_cun,A_lake_cross=a_lake_cross,A_vs_M=ssi_a_vs_m), 
       ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
  #dev.off()
  #tiff(paste("D:/Martin Lab/rna_2018/all_2018_samples/filter_comparisons/pure_cross_filter/a_vs_p_", stage, "_venn.tiff", sep = ""), width = 7, height = 7, units = 'in', res = 1000)
  venn(list(A_vs_NCA=ssi_a_vs_nca,A_vs_CUN=ssi_a_vs_cun,A_lake_cross=a_lake_cross,A_vs_P=ssi_a_vs_p), 
       ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
  #dev.off()
  
  
## Evolution letters pattern? Parallel evolution of expression?  
  library(eulerr)
  am   <- read.csv(paste(DE_genes_dir,"DE_(CRPA_&_OSPA)_vs_(CRPM_&_OSPM)_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  ap   <- read.csv(paste(DE_genes_dir,"DE_(CRPA_&_OSPA)_vs_(CRPP_&_OSPP)_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  am    <- am   [which(am$padj < 0.05),]
  ap    <- ap   [which(ap$padj < 0.05),]
  both <- merge(ap,am, by = c("tag"))
  head(both)
  nrow(both)
  nrow(both[which(both$log2FoldChange.x > 0 & both$log2FoldChange.y > 0),]) + nrow(both[which(both$log2FoldChange.x < 0 & both$log2FoldChange.y < 0),])
  nrow(both[which(both$log2FoldChange.x > 0 & both$log2FoldChange.y < 0),])
  nrow(both[which(both$log2FoldChange.x < 0 & both$log2FoldChange.y > 0),])
  
  
  #cols <- 	c("#f6c700", "#E77200", "#FF00CC","#000000", "#bd1e24", "#0067a7")
  cols <- c("green4", "blue2", "darkblue")
  #tiff(paste("D:/Martin Lab/rna_2018/all_2018_samples/filter_comparisons/pure_cross_filter/parallel_expression_venn_",stage,".tiff", sep =""), width = 7, height = 7, units = 'in', res = 1000)
  VennDiag <- euler(list(A_vs_M=ssi_a_vs_m,A_vs_P=ssi_a_vs_p)) 
  plot(VennDiag, counts = TRUE, font=1, cex=0.5, alpha=0.5,
       fill=cols,quantities = list(fontsize = 8))
  #dev.off()
       

#####
######################################################
######################################################
##### proportion DE vs FST (PLUG_N_CHUG) #############
######################################################

# error bars bar plot
{
de <- read.table("D:/Martin Lab/rna_2018/all_2018_samples/summary_stats_DE_boots_downs.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
head(de)
de$range_bs <- de$bs_1000_90q - de$bs_1000_10q
de$range_ds <- de$ds_3v3_90q - de$ds_3v3_10q
de$prop_de <- (de$DE / de$informative_transcripts)*100

max(de$range_bs)
max(de$range_ds)
de[70,]
de$names <- paste(de$pop1, de$pop2, sep = " VS ")

de1 <- de[which(de$names == "CRPA_&_OSPA VS CRPM_&_OSPM"),]
de2 <- de[which(de$names == "CRPA_&_OSPA VS CRPP_&_OSPP"),]
de3 <- de[which(de$names == "CRPA_&_OSPA VS UPxUA"),]
de4 <- de[which(de$names == "CRPA_&_OSPA VS NCA"),]
de <- rbind(de1, de2, de3, de4)
de <- de[which(de$stage == "8dpf"),]
de1 <- de[which(de$names == "CRPA_&_OSPA_&_OYSA VS CRPM_&_OSPM"),]
de2 <- de[which(de$names == "CRPA_&_OSPA_&_OYSA VS CRPP_&_OSPP_&_LILP"),]
de3 <- de[which(de$names == "CRPA_&_OSPA VS UPxUA"),]
de4 <- de[which(de$names == "CRPA_&_OSPA VS NCA"),]
de <- rbind(de1, de2, de3, de4)
de <- de[which(de$stage == "8dpf"),]


barCenters <- barplot(de$prop_de,names.arg = de$names,
                      beside = true, las = 2,
                      cex.names = 0.75, xaxt = "n",
                      ylim = c(0,(max(de$prop_de+3))),
                      border = "black", axes = TRUE)

segments(barCenters, de$ds_3v3_10q, barCenters, de$ds_3v3_90q)
}

library(DESeq2)
library(pheatmap)
library(vegan)
library(rgl)
library(ape)
library(dplyr)

sum_stats_final <- data.frame(pops1=character(), 
                              pops2=character(),
                              stage=character(),
                              n_pops1=numeric(),
                              n_pops2=numeric(),
                              informative_transcripts=numeric(),
                              de_up=numeric(),      
                              de_dn=numeric(),      
                              prop_de=numeric(),
                              de_total=numeric(),
                              median_SE=numeric(),
                              overdom_outside_parental_1_way_95=numeric(),
                              underdom_outside_parental_1_way_95=numeric(),
                              total_misexpressed_1_way_95=numeric(),
                              prop_misexpressed_1_way_95=numeric(),
                              overdom_outside_parental_2_way_95=numeric(),
                              underdom_outside_parental_2_way_95=numeric(),
                              total_misexpressed_2_way_95=numeric(),
                              prop_misexpressed_2_way_95=numeric(),
                              stringsAsFactors=FALSE)
stages <- c("48hpf", "8dpf")
for (stage in stages)
{
setwd("D:/Martin Lab/rna_2018/all_2018_samples/")
comp <- paste("all_samples_", stage, sep = "")
comp_file <- paste("D:/Martin Lab/rna_2018/all_2018_samples/", comp , ".txt", sep = "")
all_cts <- read.table("D:/Martin Lab/rna_2018/all_2018_samples/conditions/all_48hpf_and_8dpf_counts_my_gtf_geneid_final", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
sample_list <- read.table(comp_file, header = TRUE, stringsAsFactors = FALSE)
keeps <- c("Geneid", sample_list$sample)
keeper <- sample_list$sample
cts <- all_cts[keeps]
cts_out <- paste(comp, "_counts.txt", sep = "")
setwd("D:/Martin Lab/rna_2018/plug_n_chug/conditions/")
write.table(cts, cts_out, row.names = FALSE, quote= FALSE,sep="\t")
cts <- as.matrix(read.table(cts_out ,sep = "\t",header = TRUE,row.names=1)) 
colData <- as.matrix(read.table(comp_file ,header = TRUE,row.names=1))

ncol(cts)
head(cts)
nrow(colData)
# only include sequencing round in model for 48hpf
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = colData,
                              design= ~f1)
dds <- estimateSizeFactors(dds)
idx <- rowSums(counts(dds, normalized=TRUE) >= 10 ) >= (nrow(colData))
length(rowSums(counts(dds, normalized=TRUE) >= 10 )) >= (nrow(colData))
length(idx)
dds <- dds[idx,]
dds <- DESeq(dds)
norm_cts <- data.frame(counts(dds, normalized=TRUE))
head(norm_cts)
norm_cts$genes <- rownames(norm_cts)

# compare hybrids
all_comps <- read.table("D:/Martin Lab/rna_2018/fst/hybrid_mse_comparisons.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
#compare all (will take a while)
#all_comps <- read.table("D:/Martin Lab/rna_2018/all_2018_samples/comparisons.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
final_features <- read.table("D:/Martin Lab/rna_2018/all_2018_samples/features_gff.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t", quote = "")
blast_key <- read.table("D:/Martin Lab/rna_2018/all_2018_samples/blast/cyprinodon_to_danio_one_way_best_hit_symbols.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")

all_comps <- all_comps[which(all_comps$stage == stage),]
underdom <- c()
overdom <- c()
sum_stats <- data.frame(pops1=character(), 
                        pops2=character(),
                        stage=character(),
                        n_pops1=numeric(),
                        n_pops2=numeric(),
                        informative_transcripts=numeric(),
                        de_up=numeric(),      
                        de_dn=numeric(),      
                        prop_de=numeric(),
                        de_total=numeric(),
                        median_SE=numeric(),
                        overdom_outside_parental_1_way_95=numeric(),
                        underdom_outside_parental_1_way_95=numeric(),
                        total_misexpressed_1_way_95=numeric(),
                        prop_misexpressed_1_way_95=numeric(),
                        overdom_outside_parental_2_way_95=numeric(),
                        underdom_outside_parental_2_way_95=numeric(),
                        total_misexpressed_2_way_95=numeric(),
                        prop_misexpressed_2_way_95=numeric(),
                        stringsAsFactors=FALSE)
for (i in (1:nrow(all_comps)))
#subset_comps <- c(1:2)
#for (i in subset_comps)
{
  # set sample names and make count table
  pops1_name <- all_comps$pops1[i]
  pops2_name <- all_comps$pops2[i]
  pops1 <- strsplit(pops1_name,"_and_")[[1]]
  pops2 <- strsplit(pops2_name,"_and_")[[1]]
  comp_name <- paste(pops1_name, "_vs_", pops2_name, sep = "")
  master <- read.table("D:/Martin Lab/rna_2018/all_2018_samples/table_maker_master_outlier_rm.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
  head(master)
  pops1_table <- master[master$f1 %in% pops1, ]
  pops2_table <- master[master$f1 %in% pops2, ]
  pops1_table <- pops1_table[which(pops1_table$stage == stage),]
  pops2_table <- pops2_table[which(pops2_table$stage == stage),]
  keeps <- c(pops1_table$sample)
  pops1_cts <- norm_cts[keeps]
  keeps <- c(pops2_table$sample)
  pops2_cts <- norm_cts[keeps]
  #pops1_table$f1 <- pops1_name
  #pops2_table$f1 <- pops2_name
  all_contrast_pops <- rbind(pops1_table,pops2_table)
  colDataframe <- as.data.frame(colData)
  colData_with_contrast_grp_1 <- colDataframe[colDataframe$f1 %in% pops1_table$f1, ] 
  colData_with_contrast_grp_2 <- colDataframe[colDataframe$f1 %in% pops2_table$f1, ] 
  colData_without_contrast_grps <- colDataframe[! colDataframe$f1 %in% all_contrast_pops$f1, ] 
  #colData_without_contrast_grps$contrast_groups <- "not_contrasted"
  colData_with_contrast_grp_1$f1 <- pops1_name 
  colData_with_contrast_grp_2$f1 <- pops2_name 
  colDataframe <- rbind(colData_without_contrast_grps,colData_with_contrast_grp_1,colData_with_contrast_grp_2)
  colDataframe <- colDataframe[ order(row.names(colDataframe)), ]
  write.table(colDataframe, "D:/Martin Lab/rna_2018/fst/temp_table.txt")
  colData_contrasts <- as.matrix(read.table("D:/Martin Lab/rna_2018/fst/temp_table.txt",header = TRUE,row.names=1))
  cts_contrasts <- cts[ , order(colnames(cts))]
  
  # DEseq method with contrast
  #shithead
  #if (stage == "48hpf")
  #{
  #dds <- DESeqDataSetFromMatrix(countData = cts_contrasts,
  #                              colData = colData_contrasts,
  #                              design= ~f1+sequencing_round)
  #}else
  #{
  dds <- DESeqDataSetFromMatrix(countData = cts_contrasts,
                                colData = colData_contrasts,
                                design= ~f1) 
  #}
  dds <- estimateSizeFactors(dds)
  idx <- rowSums(counts(dds, normalized=TRUE) >= 10 ) >= (nrow(colData))
  #length(rowSums(counts(dds, normalized=TRUE) >= 10 )) >= (nrow(colData))
  #length(idx)
  dds <- dds[idx,]
  dds$f1 <- relevel(dds$f1, ref = pops1_name)
  dds <- DESeq(dds)
  res <- results(dds, contrast=c("f1",pops2_name,pops1_name))
  res$pvalue[res$baseMean < 10] <- NA
  res$padj <- p.adjust(res$pvalue, method="BH")
  resOrdered <- res[order(res$padj),]
  #summary(res)
  resLFC <- lfcShrink(dds, coef=2, res=res)
  res_ordered <- as.data.frame(resOrdered)
  
  total_genes <- nrow(res_ordered)
  de_total    <- nrow(res_ordered[which(res_ordered$padj < 0.05),])
  de_up       <- (nrow(res_ordered[which(res_ordered$log2FoldChange > 0 & res_ordered$padj < 0.05),]))
  de_dn       <- (nrow(res_ordered[which(res_ordered$log2FoldChange < 0 & res_ordered$padj < 0.05),]))
  prop_de     <- de_total/total_genes
  
  # calculate 95% CI for parental pops
  means <- rowMeans(pops1_cts)
  std_devs <- rowSds(as.matrix(pops1_cts))
  pops1_cts$mean <- means
  pops1_cts$std_dev <- std_devs
  j1 <- 1.96 * (pops1_cts$std_dev / sqrt(length(pops1_table$sample)))
  pops1_cts$low_95 <- pops1_cts$mean - j1
  pops1_cts$high_95 <- pops1_cts$mean + j1
  
  # semi-conservative one way - keep genes that do not overlap parental 95% CI in any single hybrid
  underdom <- c()
  overdom <- c()
  for (k in c(1:nrow(pops2_cts)))
  {
    hyb_cts <- as.numeric(pops2_cts[k,])
    gene_name <- row.names(pops2_cts[k,])
    low <- pops1_cts$low_95[k]
    hi <-  pops1_cts$high_95[k]
    if (length(hyb_cts[hyb_cts < low]) == length(hyb_cts))
    {
      underdom <- c(underdom, gene_name)
    }
    if (length(hyb_cts[hyb_cts > hi]) == length(hyb_cts))
    {
      overdom <- c(overdom, gene_name)
    }
  }
  
  # most conservative two way - keep genes that do not overlap 95% CIs between parents and hybrids
  means <- rowMeans(pops2_cts)
  std_devs <- rowSds(as.matrix(pops2_cts))
  pops2_cts$mean <- means
  pops2_cts$std_dev <- std_devs
  j2 <- 1.96 * (pops2_cts$std_dev / sqrt(length(pops2_table$sample)))
  pops2_cts$low_95 <- pops2_cts$mean - j2
  pops2_cts$high_95 <- pops2_cts$mean + j2
  underdom2 <- c()
  overdom2 <- c()
  for (k in c(1:nrow(pops2_cts)))
  {
    hyb_low <- pops2_cts$low_95[k]
    gene_name <- row.names(pops2_cts[k,])
    hyb_hi <- pops2_cts$high_95[k]
    p_low <- pops1_cts$low_95[k]
    p_hi <-  pops1_cts$high_95[k]
    if (hyb_low > p_hi)
    {
      overdom2 <- c(overdom2, gene_name)
    }
    if (hyb_hi < p_low)
    {
      underdom2 <- c(underdom2, gene_name)
    }
  }
  
  
  # MA plot
  
  sample_sizes_plot <- paste(nrow(colDataframe[which(colDataframe$f1 == pops1_name),]), " vs ", nrow(colDataframe[which(colDataframe$f1 == pops2_name),]), sep = "")
  total_genes_plot <- paste(total_genes, "transcripts", sep = " ")
  de_total_plot <- paste(de_total, "DE", sep = " ")
  prop_de_plt <- paste((100*(round(prop_de, digits = 3))), "% DE", sep = "")
  plotMA(resLFC, ylim=c(-6,5), main = comp_name)
  legend("topright", legend=c(sample_sizes_plot,total_genes_plot, de_total_plot, prop_de_plt),cex=1.0, bty = 'n')
  legend("bottomleft", legend=c(paste((100*(round(de_dn, digits = 3))), "% DE down", sep = "")),cex=0.8, bty = 'n')
  legend("topleft", legend=c(paste((100*(round(de_up, digits = 3))), "% DE up", sep = "")),cex=0.8, bty = 'n')
  #legend('bottomright', legend=c(boot_ci, dns_ci),cex=0.8, bty = 'n')
  de_plot <- paste("D:/Martin Lab/rna_2018/plug_n_chug/hybrid_mse/de_plots/",comp_name, "_",stage,"_de_plot.tiff", sep = "")
  tiff(de_plot, width = 7, height = 6, units = 'in', res = 600)
  plotMA(resLFC, ylim=c(-6,5), main = comp_name)
  legend("topright", legend=c(sample_sizes_plot,total_genes_plot, de_total_plot, prop_de_plt),cex=1.0, bty = 'n')
  legend("bottomleft", legend=c(paste((100*(round(de_dn, digits = 3))), "% DE down", sep = "")),cex=0.8, bty = 'n')
  legend("topleft", legend=c(paste((100*(round(de_up, digits = 3))), "% DE up", sep = "")),cex=0.8, bty = 'n')
  #legend('bottomright', legend=c(boot_ci, dns_ci),cex=0.8, bty = 'n')
  dev.off()
  SE_plot <- paste("D:/Martin Lab/rna_2018/plug_n_chug/hybrid_mse/de_plots/",comp_name,"_",stage, "_SE_plot.tiff", sep = "")
  SEs <- mcols(dds,use.names=TRUE)
  #mcols(mcols(dds))
  #head(SEs)
  SE_name <- paste("SE_f1_", pops2_name, "_vs_", pops1_name, sep = "")
  SE_name <- SEs[SE_name]
  tiff(SE_plot, width = 7, height = 6, units = 'in', res = 400)
  hist(SE_name[,1],main = comp_name)#, ylim = c(0,20000))
  dev.off()
  
  #### overlap with genes ###
  
  setwd("D:/Martin Lab/rna_2018/plug_n_chug/hybrid_mse/genes")
  res_ordered$Geneid <- rownames(res_ordered)
  final <- merge(res_ordered, final_features, by = c("Geneid"))
  #length(unique(final$Geneid))
  #length(unique(final$Geneid)) / length(unique(res_ordered$Geneid))
  #head(final)
  #nrow(final)
  #head(blast_key)
  final <- merge(final, blast_key,all.x = TRUE, by = c("product_accession"))
  final <- final[order(final$padj, decreasing = FALSE),]
  write.csv(final, paste(comp_name,"_",stage, "_genes.csv", sep = ""), row.names = FALSE)
  
  i_stats <- data.frame(parents=pops1_name, 
                        f1=pops2_name,
                        stage=stage,
                        n_parents=nrow(pops1_table),
                        n_f1=nrow(pops2_table),
                        informative_transcripts=nrow(pops2_cts),
                        de_up=de_up,      
                        de_dn=de_dn ,     
                        prop_de=prop_de,
                        de_total=de_total,
                        median_SE=median(SE_name[,1]),
                        overdom_outside_parental_1_way_95=length(overdom),
                        underdom_outside_parental_1_way_95=length(underdom),
                        total_misexpressed_1_way_95=(length(overdom)+length(underdom)),
                        prop_misexpressed_1_way_95=((length(overdom)+length(underdom))/nrow(pops2_cts)),
                        overdom_outside_parental_2_way_95=length(overdom2),
                        underdom_outside_parental_2_way_95=length(underdom2),
                        total_misexpressed_2_way_95=(length(overdom2)+length(underdom2)),
                        prop_misexpressed_2_way_95=((length(overdom2)+length(underdom2))/nrow(pops2_cts)),
                        stringsAsFactors=FALSE)
  sum_stats <- rbind(sum_stats, i_stats)
  
}
sum_stats_final <- rbind(sum_stats_final,sum_stats)
write.table(sum_stats_final, "D:/Martin Lab/rna_2018/plug_n_chug/summary_stats.txt", row.names = FALSE, quote= FALSE,sep="\t")
}











# incorporate downsamples
all_comps <- read.table("D:/Martin Lab/rna_2018/fst/comparisons.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
stage <- "48hpf"
all_comps <- all_comps[which(all_comps$stage == stage),]
for(i in c(1:nrow(all_comps)))
{
  # set sample names 
  pops1_name <- all_comps$pops1[i]
  pops2_name <- all_comps$pops2[i]
  pops1 <- paste("(", pops1_name, ")", sep = "")
  pops2 <- paste("(", pops2_name, ")", sep = "")
  comp <- paste("D:/Martin Lab/rna_2018/all_2018_samples/bootstrapping_and_downsampling/",pops1, "_vs_", pops2,"_",stage,"_","downsamples.txt", sep = "")
  
  read.table(comp, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
  
  
  # calculate 95% CI for %DE in hybrids vs parents
  means <- rowMeans(pops1_cts)
  std_devs <- rowSds(as.matrix(pops1_cts))
  pops1_cts$mean <- means
  pops1_cts$std_dev <- std_devs
  j1 <- 1.96 * (pops1_cts$std_dev / sqrt(length(pops1_table$sample)))
  pops1_cts$low_95 <- pops1_cts$mean - j1
  pops1_cts$high_95 <- pops1_cts$mean + j1
  
}

#mis <- read.table("D:/Martin Lab/rna_2018/fst/prop_DE_misexpression_fst_dxy.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
mis <- read.table("D:/Martin Lab/rna_2018/plug_n_chug/hybrid_mse/summary_stats_contrasts_and_subsets_popgen.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
#mis1 <- merge(mis, sum_stats_final, by = c("parents", "f1", "stage"))
mis8 <- mis[which(mis$stage =="8dpf"),] 
mis2 <- mis[which(mis$stage =="48hpf"),]
cols <- c("#91648F","#808080","#5AA27C","#f6c700","#E77200","#FF00CC","#bd1e24","#0067a7")

pdf("D:/Martin Lab/rna_2018/fst/dxy_vs_prop_mse.pdf",width=7,height=5)
plt <- ggplot(mis2, aes(x=dxy_dna, y=prop_DE_contrasts, color = cross_type)) + geom_point(aes(shape=lake),size = 4)
plt + theme_classic() +scale_color_manual(values = cols) +labs(x = "\ngenome-wide mean dxy", y = "% genes misexpressed in F1 hybrids\n") +guides(col = guide_legend(override.aes = list(shape = 19, size = 5))) +ggtitle("48 hpf misexpressed contrasts")+theme(plot.title = element_text(hjust = 0.5))
plt <- ggplot(mis8, aes(x=dxy_dna, y=prop_DE_contrasts, color = cross_type)) + geom_point(aes(shape=lake),size = 4)
plt + theme_classic() +scale_color_manual(values = cols) +labs(x = "\ngenome-wide mean dxy", y = "% genes misexpressed in F1 hybrids\n") +guides(col = guide_legend(override.aes = list(shape = 19, size = 5))) +ggtitle("8dpf misexpressed contrasts")+theme(plot.title = element_text(hjust = 0.5))

plt <- ggplot(mis2, aes(x=dxy_dna, y=median_SE, color = cross_type)) + geom_point(aes(shape=lake),size = 4)
plt + theme_classic() +scale_color_manual(values = cols) +labs(x = "\ngenome-wide mean dxy", y = "median standard error of gene expression\n") +guides(col = guide_legend(override.aes = list(shape = 19, size = 5))) +ggtitle("48 hpf SE")+theme(plot.title = element_text(hjust = 0.5))
plt <- ggplot(mis8, aes(x=dxy_dna, y=median_SE, color = cross_type)) + geom_point(aes(shape=lake),size = 4)
plt + theme_classic() +scale_color_manual(values = cols) +labs(x = "\ngenome-wide mean dxy", y = "median standard error of gene expression\n") +guides(col = guide_legend(override.aes = list(shape = 19, size = 5))) +ggtitle("8 dpf SE")+theme(plot.title = element_text(hjust = 0.5))

plt <- ggplot(mis2, aes(x=dxy_dna, y=prop_DE_subset, color = cross_type)) + geom_point(aes(shape=lake),size = 4)
plt + theme_classic() +scale_color_manual(values = cols) +labs(x = "\ngenome-wide mean dxy", y = "% genes misexpressed in F1 hybrids\n") +guides(col = guide_legend(override.aes = list(shape = 19, size = 5))) +ggtitle("48 hpf misexpressed subset")+theme(plot.title = element_text(hjust = 0.5))
plt <- ggplot(mis8, aes(x=dxy_dna, y=prop_DE_subset, color = cross_type)) + geom_point(aes(shape=lake),size = 4)
plt + theme_classic() +scale_color_manual(values = cols) +labs(x = "\ngenome-wide mean dxy", y = "% genes misexpressed in F1 hybrids\n") +guides(col = guide_legend(override.aes = list(shape = 19, size = 5))) +ggtitle("8dpf misexpressed subset")+theme(plot.title = element_text(hjust = 0.5))

plt <- ggplot(mis2, aes(x=dxy_dna, y=prop_misexpressed_2_way_95, color = cross_type)) + geom_point(aes(shape=lake),size = 4)
plt + theme_classic() +scale_color_manual(values = cols) +labs(x = "\ngenome-wide mean dxy", y = "% genes misexpressed in F1 hybrids\n") +guides(col = guide_legend(override.aes = list(shape = 19, size = 5))) +ggtitle("48 hpf misexpressed no 95 CI overlap")+theme(plot.title = element_text(hjust = 0.5))
plt <- ggplot(mis8, aes(x=dxy_dna, y=prop_misexpressed_2_way_95, color = cross_type)) + geom_point(aes(shape=lake),size = 4)
plt + theme_classic() +scale_color_manual(values = cols) +labs(x = "\ngenome-wide mean dxy", y = "% genes misexpressed in F1 hybrids\n") +guides(col = guide_legend(override.aes = list(shape = 19, size = 5))) +ggtitle("8dpf misexpressed no 95 CI overlap")+theme(plot.title = element_text(hjust = 0.5))


dev.off()



#all_comps <- read.table("D:/Martin Lab/rna_2018/all_2018_samples/bootstrapping_and_downsampling/boot_and_downs_comparisons.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
head(all_comps)
sum_stats <- data.frame(pop1=character(), 
                        pop2=character(),
                        stage=character(),
                        n_pop1=numeric(),
                        n_pop2=numeric(),
                        informative_transcripts=numeric(),
                        DE=numeric(),
                        up=numeric(),
                        down=numeric(),
                        stringsAsFactors=FALSE)
for (i in (1:nrow(all_comps)))
#for (i in subset_comps)
{
  pops1_name <- all_comps$pops1[i]
  pops2_name <- all_comps$pops2[i]
  pops1_t <- pops1_name
  pops2_t <- pops2_name
  pops1 <- strsplit(pops1_name,"_&_")[[1]]
  pops2 <- strsplit(pops2_name,"_&_")[[1]]
  stage <- all_comps$stage[i]
  pops1_name <- paste("(", pops1_name, ")", sep = "")
  pops2_name <- paste("(", pops2_name, ")", sep = "")
  comp_name <- paste("condition_species_",pops1_name, "_vs_", pops2_name, "_" ,stage,".txt", sep = "")
  single_run <- paste(pops1_name, "_vs_", pops2_name, "_" ,stage, sep = "")
  
  setwd("D:/Martin Lab/rna_2018/fst/conditions/")
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
  
  setwd("D:/Martin Lab/rna_2018/fst/conditions/")
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
  # most stringent -- filter out genes if read count is zero for more than half of individuals.
  
  dds <- estimateSizeFactors(dds)
  idx <- rowSums(counts(dds, normalized=TRUE) >= 10 ) >= (nrow(colData))
  #length(rowSums(counts(dds, normalized=TRUE) >= 1 ) >= (nrow(colData)/2))
  length(idx)
  dds <- dds[idx,]
  dds <- DESeq(dds)
  

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
  
  #test differences in variance
  vars <- mcols(dds,use.names=TRUE)
  mcols(mcols(dds))
  head(vars)
  hist(vars$baseVar, breaks = 100000, xlim = c(0, 200000), ylim = c(0,20000))

  
  total_genes <- nrow(res_ordered)
  de_total <- nrow(res_ordered[which(res_ordered$padj < 0.05),])
  de_up <- (nrow(res_ordered[which(res_ordered$log2FoldChange > 0 & res_ordered$padj < 0.05),]))/total_genes
  de_dn <- (nrow(res_ordered[which(res_ordered$log2FoldChange < 0 & res_ordered$padj < 0.05),]))/total_genes
  prop_de <- de_total/total_genes
  sample_sizes_plot <- paste(sample_size_a, " vs ", sample_size_b, sep = "")
  

  
  i_stats <- data.frame(pop1=pops1_t, 
                        pop2=pops2_t,
                        stage=stage,
                        n_pop1=sample_size_a,
                        n_pop2=sample_size_b,
                        informative_transcripts=total_genes,
                        DE=de_total,
                        up=nrow(res_ordered[which(res_ordered$log2FoldChange < 0 & res_ordered$padj < 0.05),]),
                        down=nrow(res_ordered[which(res_ordered$log2FoldChange > 0 & res_ordered$padj < 0.05),]),
                        stringsAsFactors=FALSE)
  sum_stats <- rbind(sum_stats, i_stats)

  
}
sum_stats$prop_DE <- sum_stats$DE / sum_stats$informative_transcripts
#write.table(sum_stats, "D:/Martin Lab/rna_2018/fst/summary_stats.txt", row.names = FALSE, quote= FALSE,sep="\t")
sum_stats

library(ggplot2)
mis <- read.table("D:/Martin Lab/rna_2018/fst/prop_DE_misexpression_fst_dxy.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
mis8 <- mis[which(mis$stage =="8dpf"),] 
mis2 <- mis[which(mis$stage =="48hpf"),]

cols <- c("#91648F","#808080","#5AA27C","#f6c700","#E77200","#FF00CC","#bd1e24","#0067a7")

# purple =  "#91648F"
# grey =    "#808080"
# green =   "#5AA27C"
# yellow =  "#f6c700"
# orange =  "#E77200"
# pink =    "#FF00CC"
# black =   "#000000"
# red =     "#bd1e24"
# blue =    "#0067a7"

#pdf("D:/Martin Lab/rna_2018/fst/fst_vs_prop_mse.pdf",width=7,height=5)
{
plt <- ggplot(mis8, aes(x=fst, y=prop_DE, color = cross_type)) + geom_point(aes(shape=lake),size = 4)
plt +   theme_classic() +scale_color_manual(values = cols) +labs(x = "\ntranscriptome-wide mean Fst", y = "% genes misexpressed in F1 hybrids\n") +
  guides(col = guide_legend(override.aes = list(shape = 19, size = 5)))+ggtitle("8 dpf") +
  theme(plot.title = element_text(hjust = 0.5))#+
  #geom_segment(aes(x = mis8$fst, y = mis8$prop_DE, xend = mis2$fst, yend = mis2$prop_DE),show.legend = FALSE)
plt <- ggplot(mis2, aes(x=fst, y=prop_DE, color = cross_type)) + geom_point(aes(shape=lake),size = 4)
plt + theme_bw() +scale_color_manual(values = cols) +labs(x = "\ntranscriptome-wide mean Fst", y = "% genes misexpressed in F1 hybrids\n") +
  guides(col = guide_legend(override.aes = list(shape = 19, size = 5))) +ggtitle("48 hpf")+
  theme(plot.title = element_text(hjust = 0.5))
#plt <- ggplot(mis8, aes(x=dxy, y=prop_DE, color = cross_type)) + geom_point(aes(shape=lake),size = 4)
#plt + theme_bw() +scale_color_manual(values = cols) +labs(x = "\ntranscriptome-wide mean dxy", y = "% genes misexpressed in F1 hybrids\n") +
#  guides(col = guide_legend(override.aes = list(shape = 19, size = 5)))+ggtitle("8 dpf") +
#  theme(plot.title = element_text(hjust = 0.5))
#plt <- ggplot(mis2, aes(x=dxy, y=prop_DE, color = cross_type)) + geom_point(aes(shape=lake),size = 4)
#plt + theme_bw() +scale_color_manual(values = cols) +labs(x = "\ntranscriptome-wide mean dxy", y = "% genes misexpressed in F1 hybrids\n") +
#  guides(col = guide_legend(override.aes = list(shape = 19, size = 5))) +ggtitle("48 hpf")+
#  theme(plot.title = element_text(hjust = 0.5))

plt <- ggplot(mis8, aes(x=fst, y=prop_misexpressed_outside_95, color = cross_type)) + geom_point(aes(shape=lake),size = 4)
plt + theme_classic() +scale_color_manual(values = cols) +labs(x = "\ntranscriptome-wide mean Fst", y = "% genes misexpressed in F1 hybrids\n") +
  guides(col = guide_legend(override.aes = list(shape = 19, size = 5)))+ggtitle("8 dpf misexpressed outside parental 95 CI") +
  theme(plot.title = element_text(hjust = 0.5))#+
  #geom_segment(aes(x = mis8$fst, y = mis8$prop_misexpressed_outside_95, xend = mis2$fst, yend = mis2$prop_misexpressed_outside_95),show.legend = FALSE)

plt <- ggplot(mis2, aes(x=fst, y=prop_misexpressed_outside_95, color = cross_type)) + geom_point(aes(shape=lake),size = 4)
plt + theme_classic() +scale_color_manual(values = cols) +labs(x = "\ntranscriptome-wide mean Fst", y = "% genes misexpressed in F1 hybrids\n") +
  guides(col = guide_legend(override.aes = list(shape = 19, size = 5))) +ggtitle("48 hpf misexpressed outside parental 95 CI")+
  theme(plot.title = element_text(hjust = 0.5))
#plt <- ggplot(mis8, aes(x=dxy, y=prop_misexpressed_outside_95, color = cross_type)) + geom_point(aes(shape=lake),size = 4)
#plt + theme_bw() +scale_color_manual(values = cols) +labs(x = "\ntranscriptome-wide mean dxy", y = "% genes misexpressed in F1 hybrids\n") +
#  guides(col = guide_legend(override.aes = list(shape = 19, size = 5)))+ggtitle("8 dpf mis outside parental 95") +
#  theme(plot.title = element_text(hjust = 0.5))
#plt <- ggplot(mis2, aes(x=dxy, y=prop_misexpressed_outside_95, color = cross_type)) + geom_point(aes(shape=lake),size = 4)
#plt + theme_bw() +scale_color_manual(values = cols) +labs(x = "\ntranscriptome-wide mean dxy", y = "% genes misexpressed in F1 hybrids\n") +
#  guides(col = guide_legend(override.aes = list(shape = 19, size = 5))) +ggtitle("48 hpf mis outside parental 95")+
#  theme(plot.title = element_text(hjust = 0.5))
}
#dev.off()


#####
######################################################
######################################################
##### ASE ############################################
######################################################

library(reshape2)
inheritance_comps <- read.table("D:/Martin Lab/rna_2018/all_2018_samples/inheritance_patterns/inheritance_comparisons.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
DE_genes_dir <- "D:/Martin Lab/rna_2018/all_2018_samples/conditions/"
ase_plots_dir <- "D:/Martin Lab/rna_2018/ASE/plots/"
geneiase_dir <- "C:/Users/jmcgirr/Documents/all_2018_samples/ase_data/geneiase/"
mbased_dir <- "C:/Users/jmcgirr/Documents/all_2018_samples/ase_data/mbased/"
maternal_counts_dir <- "C:/Users/jmcgirr/Documents/all_2018_samples/ase_data/maternal_counts/"
master <- read.table("D:/Martin Lab/rna_2018/all_2018_samples/table_maker_master_outlier_rm.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
norm_cts <- read.table("D:/Martin Lab/rna_2018/all_2018_samples/WGCNA/all_samples_outlier_rm_size_factor_normalized_counts.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
all_snps <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/ase_data/CAE1_snp_table.txt", header = TRUE, stringsAsFactors = FALSE)
all_snps <- all_snps[, -grep("HP", colnames(all_snps))]
head(all_snps)
#i <- 1
subset_comps <- c(15:18)

for (i in subset_comps)
  #for (i in (1:nrow(inheritance_comps)))
{
  parent_pops1_name <- inheritance_comps$parent_pops1[i]
  parent_pops2_name <- inheritance_comps$parent_pops2[i]
  hy_name <- inheritance_comps$hybrids[i]
  stage <- inheritance_comps$stage[i]
  parent_pops_v_hy  <- paste(DE_genes_dir, "DE_(", parent_pops1_name,"_&_",parent_pops2_name, ")_vs_(",hy_name, ")_",stage,"_genes",".csv",sep = "")
  parent_pops1_v_hy <- paste(DE_genes_dir, "DE_(", parent_pops1_name, ")_vs_(",hy_name, ")_",stage,"_genes",".csv",sep = "")
  parent_pops2_v_hy <- paste(DE_genes_dir, "DE_(", parent_pops2_name, ")_vs_(",hy_name, ")_",stage,"_genes",".csv",sep = "")
  parent_pops1_v_parent_pops2 <- paste(DE_genes_dir, "DE_(", parent_pops1_name, ")_vs_(",parent_pops2_name, ")_",stage,"_genes",".csv",sep = "")
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
  #shitter
  
  p1Ase <- c()
  for (ind in p1_inds)
    #ind <- "CAE1"
  {
    gAse <- read.table(paste(geneiase_dir,ind, "_geneiase_ase.txt.static.gene.pval.tab",sep = ""), header = TRUE, stringsAsFactors = FALSE) 
    mAse <- read.table(paste(mbased_dir,ind, "_mbased_ase.txt",sep = ""), header = TRUE, stringsAsFactors = FALSE) 
    
    gAse <- gAse[which(gAse$p.nom <= 0.05),]
    mAse <- mAse[which(mAse$pValueASE <= 0.05),]
    gAse <- cbind(gAse, colsplit(gAse$feat, ";", c("related_accession", "gene_name")))
    mAse <- cbind(mAse, colsplit(mAse$mrnaID, ";", c("related_accession", "gene_name"))) 
    ase <- merge(gAse,mAse, by = ("related_accession"))
    nrow(ase)
    nrow(mAse)
    nrow(gAse)
    p1Ase <- c(p1Ase, ase$related_accession)
  }
  p2Ase <- c()
  for (ind in p2_inds)
  {
    gAse <- read.table(paste(geneiase_dir,ind, "_geneiase_ase.txt.static.gene.pval.tab",sep = ""), header = TRUE, stringsAsFactors = FALSE) 
    mAse <- read.table(paste(mbased_dir,ind, "_mbased_ase.txt",sep = ""), header = TRUE, stringsAsFactors = FALSE) 
    
    gAse <- gAse[which(gAse$p.nom <= 0.05),]
    mAse <- mAse[which(mAse$pValueASE <= 0.05),]
    gAse <- cbind(gAse, colsplit(gAse$feat, ";", c("related_accession", "gene_name")))
    mAse <- cbind(mAse, colsplit(mAse$mrnaID, ";", c("related_accession", "gene_name"))) 
    ase <- merge(gAse,mAse, by = ("related_accession"))

    p2Ase <- c(p2Ase, ase$related_accession)
  }
  hyAse <- data.frame(related_accession=character(), 
                          mrnaID=character(),
                          stringsAsFactors=FALSE)
  ind <- hy_inds[[1]]
  gAse <- read.table(paste(geneiase_dir,ind, "_geneiase_ase.txt.static.gene.pval.tab",sep = ""), header = TRUE, stringsAsFactors = FALSE) 
  mAse <- read.table(paste(mbased_dir,ind, "_mbased_ase.txt",sep = ""), header = TRUE, stringsAsFactors = FALSE) 
  
  gAse <- gAse[which(gAse$p.nom <= 0.05),]
  mAse <- mAse[which(mAse$pValueASE <= 0.05),]
  gAse <- cbind(gAse, colsplit(gAse$feat, ";", c("related_accession", "gene_name")))
  mAse <- cbind(mAse, colsplit(mAse$mrnaID, ";", c("related_accession", "gene_name"))) 
  ase <- merge(gAse,mAse, by = c("related_accession"))
  keeps <- c("related_accession", "mrnaID")
  hyAse <- ase[keeps]
  
  for (ind in hy_inds)
  {
    gAse <- read.table(paste(geneiase_dir,ind, "_geneiase_ase.txt.static.gene.pval.tab",sep = ""), header = TRUE, stringsAsFactors = FALSE) 
    mAse <- read.table(paste(mbased_dir,ind, "_mbased_ase.txt",sep = ""), header = TRUE, stringsAsFactors = FALSE) 
    
    gAse <- gAse[which(gAse$p.nom <= 0.05),]
    mAse <- mAse[which(mAse$pValueASE <= 0.05),]
    gAse <- cbind(gAse, colsplit(gAse$feat, ";", c("related_accession", "gene_name")))
    mAse <- cbind(mAse, colsplit(mAse$mrnaID, ";", c("related_accession", "gene_name"))) 
    ase <- merge(gAse,mAse, by = ("related_accession"))
    keeps <- c("related_accession", "mrnaID")
    ase <- ase[keeps]
    
    hyAse <- merge(hyAse, ase, by = c("related_accession", "mrnaID"))
    #hyAse <- c(hyAse, ase$related_accession)
  }
  
  hyAse$ase <- "yes"
  inheritance_and_ase <- merge(hyAse, all_genes,all = TRUE, by = c("related_accession"))
  inheritance_and_ase$ase[is.na(inheritance_and_ase$ase)] <- "no"
  #nrow(inheritance_and_ase[which (inheritance_and_ase$ase == "yes"),])
  
  trans <- inheritance_and_ase[which(inheritance_and_ase$padj_p1p2 <= 0.05 & inheritance_and_ase$ase == "no"),]
  nrow(trans)
  length(unique(trans$tag))
  trans$ase_type <- 5
  mis_over <- inheritance_and_ase[which(inheritance_and_ase$padj_ph <= 0.05 & inheritance_and_ase$lfc_ph > 0 & inheritance_and_ase$ase == "no"),]
  mis_over$ase_type <- 2
  nrow(mis_over)
  mis_under <- inheritance_and_ase[which(inheritance_and_ase$padj_ph <= 0.05 & inheritance_and_ase$lfc_ph < 0 & inheritance_and_ase$ase == "no"),]
  mis_under$ase_type <- 3
  nrow(mis_under)
  #mis_ase_over <- inheritance_and_ase[which(inheritance_and_ase$padj_p1p2 > 0.05 & inheritance_and_ase$padj_ph <= 0.05 & inheritance_and_ase$lfc_ph > 0 & inheritance_and_ase$ase_in_hybrids_not_parents == "yes"),]
  #nrow(mis_ase_over)
  #mis_ase_over$ase_type <- 6
  #mis_ase_under <- inheritance_and_ase[which(inheritance_and_ase$padj_p1p2 > 0.05 & inheritance_and_ase$padj_ph <= 0.05 & inheritance_and_ase$lfc_ph < 0 & inheritance_and_ase$ase_in_hybrids_not_parents == "yes"),]
  #nrow(mis_ase_under)
  #mis_ase_under$ase_type <- 6
  mis_ase_over <- inheritance_and_ase[which(inheritance_and_ase$padj_p1p2 > 0.05  & inheritance_and_ase$lfc_ph > 0 & inheritance_and_ase$ase == "yes"),]
  nrow(mis_ase_over)
  mis_ase_over$ase_type <- 4
  mis_ase_under <- inheritance_and_ase[which(inheritance_and_ase$padj_p1p2 > 0.05 & inheritance_and_ase$lfc_ph < 0 & inheritance_and_ase$ase == "yes"),]
  nrow(mis_ase_under)
  mis_ase_under$ase_type <- 4
  
  comp <- inheritance_and_ase[which(inheritance_and_ase$padj_p1p2 > 0.05 & inheritance_and_ase$ase == "yes"),]
  nrow(comp)
  
  con <- inheritance_and_ase[which(inheritance_and_ase$padj_ph > 0.05 & inheritance_and_ase$padj_p1p2 > 0.05 & inheritance_and_ase$ase == "no"),]
  nrow(con)
  con$ase_type <- 1
  cis <- inheritance_and_ase[which(inheritance_and_ase$padj_p1p2 <= 0.05 & inheritance_and_ase$ase == "yes"),]
  nrow(cis)
  cis$ase_type <- 7
  
  ase_type <- rbind(con, mis_over, mis_under,cis, trans, mis_ase_over, mis_ase_under)
  head(ase_type)
  nrow(ase_type)
  unique(ase_type$ase_type)
  length(unique(ase_type$tag))
  
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
  
  
  
  
  
}  

#####
########################################
########################################
##### blast geneID conversion ##########
########################################

#source("http://bioconductor.org/biocLite.R")
#biocLite("biomaRt")
##install.packages("biomaRt")
#library(biomaRt)
##listMarts()
##attributes = listAttributes(ensembl)
#
#ensembl = useMart("ensembl", dataset = "drerio_gene_ensembl")
#dat = getBM(attributes = c("zfin_id_symbol", "refseq_peptide", "refseq_peptide_predicted"), values = "*", mart = ensembl)
#head(dat,50)
#nrow(dat)
#
#blast_hits <- read.table("D:/Martin Lab/rna_2018/all_2018_samples/blast/cyprinodon_to_danio_one_way_best_hit.txt", header = FALSE, stringsAsFactors = FALSE, sep = "\t")
#head(blast_hits)
#blast_hits$cyp_peptide <- blast_hits$V1
#blast_hits$zeb_peptide <- blast_hits$V2
#blast_hits$evalue <- blast_hits$V11
#keeps <- c("cyp_peptide", "zeb_peptide", "evalue")
#blast_hits <- blast_hits[keeps]
#
#for (i in nrow(blast_hits))
#{
#  
#  cyp <- blast_hits$cyp_peptide[i]
#  zeb <- blast_hits$zeb_peptide[i]
#  e <- blast_hits$evalue[i]
#  
#}

blast_hits <- read.table("D:/Martin Lab/rna_2018/all_2018_samples/blast/cyprinodon_to_danio_one_way_best_hit.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
nrow(blast_hits)
head(blast_hits)
#blast_hits <- blast_hits[!duplicated(blast_hits$V1),]
zeb_feats <- read.table("D:/Martin Lab/rna_2018/all_2018_samples/blast/ref_GRCz11_top_level.gff3", header = FALSE, stringsAsFactors = FALSE, sep = "\t")
blast_hits$product_accession <- blast_hits$qseqid
blast_hits$zeb_peptide <- blast_hits$sseqid
keeps <- c("product_accession", "zeb_peptide", "evalue")
blast_hits <- blast_hits[keeps]
blast_hits$zeb_peptide <- sub("\\..*","",blast_hits$zeb_peptide)
zeb_feats1 <- zeb_feats[which(zeb_feats$V3 == "CDS"),] 
zeb_feats1 <- cbind(zeb_feats1, colsplit(zeb_feats1$V9, "protein_id=", c("junk", "zeb_peptide")))
zeb_feats1 <- cbind(zeb_feats1, colsplit(zeb_feats1$V9, "gene=", c("junk", "gene_junk")))
zeb_feats1 <- cbind(zeb_feats1, colsplit(zeb_feats1$gene_junk, ";", c("zeb_gene_symbol", "shit")))
zeb_feats1$zeb_peptide <- sub("\\..*","",zeb_feats1$zeb_peptide)
keeps <- c("zeb_peptide", "zeb_gene_symbol")
zeb_feats1 <- zeb_feats1[keeps]
zeb_feats <- unique(zeb_feats1)
head(zeb_feats)
head(blast_hits)
blast_key <- merge(blast_hits, zeb_feats, by = c("zeb_peptide"))
blast_key <- unique(blast_key)
nrow(blast_key)
head(blast_key)
nrow(unique(blast_key))
blast_key$zeb_peptide_one_way <- blast_key$zeb_peptide
blast_key$evalue_one_way <- blast_key$evalue
blast_key$zeb_gene_symbol_one_way <- blast_key$zeb_gene_symbol
keeps <- c("zeb_peptide_one_way", "product_accession","evalue_one_way", "zeb_gene_symbol_one_way")
blast_key <- blast_key[keeps]

#write.table(blast_key, "D:/Martin Lab/rna_2018/all_2018_samples/blast/cyprinodon_to_danio_one_way_best_hit_symbols.txt", row.names = FALSE, quote = FALSE, sep = "\t")

blast_key$match <- paste(blast_key$product_accession, blast_key$zeb_peptide_one_way, sep = ":")
zeb_to_cyp <- read.table("D:/Martin Lab/rna_2018/all_2018_samples/blast/danio_to_cyprinodon_one_way_best_hit.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
head(zeb_to_cyp)
zeb_to_cyp$qseqid <- sub("\\..*","",zeb_to_cyp$qseqid)
zeb_to_cyp$match <- paste(zeb_to_cyp$sseqid, zeb_to_cyp$qseqid, sep = ":")
zeb_to_cyp$evalue_zeb_to_cyp <- zeb_to_cyp$evalue
keeps <- c("match","evalue_zeb_to_cyp")
zeb_to_cyp <- zeb_to_cyp[keeps]
nrow(blast_key)
nrow(zeb_to_cyp)
head(zeb_to_cyp)
head(blast_key)
two_way <- merge(blast_key, zeb_to_cyp, by = c("match"))
two_way$zeb_gene_symbol_two_way <- two_way$zeb_gene_symbol_one_way
keeps <- c("match","evalue_zeb_to_cyp", "zeb_gene_symbol_two_way")
two_way <- two_way[keeps]
head(two_way)
two_way <- merge(blast_key, two_way,all.x = TRUE, by = c("match"))

nrow(two_way)
head(two_way)
two_way$match <- NULL
#write.table(two_way, "D:/Martin Lab/rna_2018/all_2018_samples/blast/cyprinodon_to_danio_two_way_best_hit_symbols.txt", row.names = FALSE, quote = FALSE, sep = "\t")


#####
########################################
########################################
##### transcript integrity numbers #####
########################################

setwd("D:/Martin Lab/RNA-seq/axm/post_reviews/tins/")
tin_files <- c()
for (tin_file in tin_files)
{
  
  
  
}


#####
########################################
########################################
##### seqmonk DNA contamination qc #####
########################################
########################################

cts <- read.table("D:/Martin Lab/rna_2018/seqmonk/counts_round_1_confirm.txt", header = FALSE, stringsAsFactors = FALSE, sep = "\t")
head(cts)

cts[] <- lapply(cts, gsub, pattern = "Mean ", replacement = "", fixed = TRUE)
cts[] <- lapply(cts, gsub, pattern = ".rna.sort.dedup.bam", replacement = "", fixed = TRUE)
colnames(cts) <-  cts[1, ] # the first row will be the header
cts <-  cts[-1, ]          # removing the first row.


#round 1
#cts <- cts[, -grep("StDev", colnames(cts))]
#cts <- cts[, -grep("Average.No.value", colnames(cts))]
#cts <- cts[, -grep("Features", colnames(cts))]
#cts <- cts[, -grep("Descriptions", colnames(cts))]
#cts <- cts[, -grep("Probes", colnames(cts))]
#cts$Length <- as.numeric(cts$End) - as.numeric(cts$Start)
#cts$genes <- "gene"
#cts$number <- c(1:(nrow(cts)))
#cts$Geneid <- paste(cts$genes, cts$number, sep = "")
#
#cts1 <- data.frame(cts$Start, cts$End, cts$Length)
#colnames(cts1) <-  c("Start", "End", "Length")
#head(cts1)
#head(cts2)
#nrow(cts2)
#drops <- c("genes","number", "Geneid", "Chr", "Start", "End", "Length")
#cts2 <- cts[ , !(names(cts) %in% drops)]
#cts2 <- data.frame(sapply(cts2, function(x) as.integer(as.character(x))))
#head(cts2)
#cts3 <- data.frame(cts$Geneid, cts$Chr)
#colnames(cts3) <-  c("Geneid","Chr")
#cts1 <- data.frame(sapply(cts1, function(x) as.integer(as.character(x))))
#cts_final <- cbind(cts3, cts1, cts2)
#head(cts_final)

#round 1 confirm
head(cts)
cts <- cts[, -grep("Feature", colnames(cts))]
cts <- cts[, -grep("Description", colnames(cts))]
cts <- cts[, -grep("ID", colnames(cts))]
cts <- cts[, -grep("No value", colnames(cts))]
cts <- cts[, -grep("Probe", colnames(cts))]
cts <- cts[, -grep("Type", colnames(cts))]
cts <- cts[, -grep("Orientation", colnames(cts))]
cts <- cts[, -grep("Distance", colnames(cts))]
cts$Length <- as.numeric(cts$End) - as.numeric(cts$Start)
cts$genes <- "gene"
cts$number <- c(1:(nrow(cts)))
cts$Geneid <- paste(cts$genes, cts$number, sep = "")

cts1 <- data.frame(cts$Start, cts$End, cts$Length)
colnames(cts1) <-  c("Start", "End", "Length")
head(cts1)
head(cts2)
nrow(cts2)
drops <- c("genes","number", "Geneid", "Chr", "Start", "End", "Length")
cts2 <- cts[ , !(names(cts) %in% drops)]
cts2 <- data.frame(sapply(cts2, function(x) as.integer(as.character(x))))
head(cts2)
cts3 <- data.frame(cts$Geneid, cts$Chr)
colnames(cts3) <-  c("Geneid","Chr")
cts1 <- data.frame(sapply(cts1, function(x) as.integer(as.character(x))))
cts_final <- cbind(cts3, cts1, cts2)
cts_final <- cts_final[, -grep("Chromosome", colnames(cts_final))]
head(cts_final)

total_genes_lists <- c()
de_total_lists <- c()
prop_de_lists <- c()
sample_order_lists <- c()

comp_files <- c("crpA_&_ncA_vs_ncAxcrpA_48hpf",
                "crpA_vs_cunAxcrpA_48hpf",
                "ospA_&_crpA_&_oysA_&_ncA_vs_ospM_8dpf",
                "ospA_&_crpA_&_oysA_&_ncA_vs_ospM_48hpf",
                "ospA_&_crpA_&_oysA_&_ncA_vs_ospP_&_crpP_&_lilP_8dpf",
                "ospA_&_crpA_&_oysA_vs_ncA_8dpf",
                "ospA_&_crpA_&_oysA_vs_ncA_48hpf",
                "ospA_&_crpA_&_oysA_vs_ospM_48hpf",
                "ospA_&_crpA_&_oysA_vs_ospP_&_crpP_&_lilP_8dpf",
                "ospA_&_crpA_&_oysA_vs_ospP_&_crpP_48hpf",
                "ospA_&_crpA_&_vs_ospP_&_crpP_48hpf",
                "ospA_&_crpA_oysA_&_ncA_vs_ospP_&_crpP_48hpf",
                "ospA_&_crpA_vs_crpAxospA_&_ospAxcrpA_48hpf",
                "ospA_&_crpA_vs_ospAxospP_&_ospPxospA_crpAxcrpP_48hpf",
                "ospA_&_crpA_vs_ospM_8dpf",
                "ospA_&_crpA_vs_ospP_&_crpP_8dpf",
                "ospA_&_ospP_vs_ospPxospA_8dpf",
                "ospA_vs_crpA_8dpf",
                "ospA_vs_crpA_48hpf")

#comp_file <- "ospA_vs_crpA_48hpf"


for(comp_file in comp_files)
{  
  
  #comp_file <- "crpA_&_ncA_vs_ncAxcrpA_48hpf"
  setwd("D:/Martin Lab/rna_2018/seqmonk/conditions/")
  
  
  col_data <-             paste("condition_species_", comp_file, ".txt", sep = "")
  cts_data <-             paste("DESeq_counts_", comp_file, ".txt", sep = "")
  cts_data_genes <-       paste("DE_", comp_file, ".csv", sep = "")
  genes_out <-            paste("DE_", comp_file, "_genes.csv", sep = "")
  #ham_mbe_cans_out <-                       "_mbe_cans.csv"
  #ham_indles_out <-                         "_indels.csv"
  #ham_snps_out <-                           "_snps.csv"
  
  #
  sample_list <- read.table(col_data, header = TRUE, stringsAsFactors = FALSE)
  keeps <- c("Geneid", sample_list$sample)
  keeper <- sample_list$sample
  cts <- cts_final
  
  cts <- cts[keeps]
  head(cts)
  #write.table(cts, cts_data, row.names = FALSE, quote= FALSE,sep="\t") 
  
  
  
  #### DIFFERENTIAL EXPRESSION ANALYSIS ####
  #vignette("DESeq2")
  
  cts <- as.matrix(read.table(cts_data ,sep = "\t",header = TRUE,row.names=1))
  head(cts)
  nrow(cts)
  
  colData <- as.matrix(read.table(col_data ,header = TRUE,row.names=1))
  head(colData)
  
  
  ncol(cts)
  nrow(colData)
  #comp <- "condition"
  dds <- DESeqDataSetFromMatrix(countData = cts,
                                colData = colData,
                                design= ~species)
  
  ### count number of unused genes
  #unused <- dds[ rowSums(counts(dds)) < 1, ]
  #unused
  # most stringent -- filter out genes if read count is zero for more than half of individuals.
  
  dds <- estimateSizeFactors(dds)
  idx <- rowSums(counts(dds, normalized=TRUE) >= 1 ) >= (nrow(colData)/2)
  length(rowSums(counts(dds, normalized=TRUE) >= 1 ) >= (nrow(colData)/2))
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
  
  total_genes <- nrow(res_ordered)
  de_total <- nrow(res_ordered[which(res_ordered$padj < 0.05),])
  prop_de <- de_total/total_genes
  total_genes <- paste("of", total_genes, "informative genes", sep = " ")
  de_total <- paste(de_total, "DE", sep = " ")
  prop_de <- paste((100*(round(prop_de, digits = 3))), "% DE", sep = "")
  plotMA(resLFC, ylim=c(-5,5), main = comp_file)
  legend(10, -4, legend=c(total_genes, de_total, prop_de),cex=1.0)
  #head(resLFC)
  de_plot <- paste("D:/Martin Lab/rna_2018/seqmonk/plots/",comp_file, "_de_plot.tiff", sep = "")
  #tiff(de_plot, width = 7, height = 6, units = 'in', res = 1000)
  plotMA(resLFC, ylim=c(-7,5), main = comp_file)
  legend(1, -4.5, legend=c(total_genes, de_total, prop_de),cex=1.0)
  #dev.off()
  
  total_genes_lists <- c(total_genes, total_genes_lists)
  de_total_lists <- c(de_total, de_total_lists)
  prop_de_lists <- c(prop_de, prop_de_lists)
  sample_order_lists <- c(comp_file, sample_order_lists)
  
  #}
  prop_de_table_seq_monk_confirm <- data.frame(comparison=sample_order_lists,
                                               genes_analyzed=total_genes_lists,
                                               genes_diff_expressed=de_total_lists,
                                               proportion_diff_expressed=prop_de_lists,
                                               stringsAsFactors=FALSE) 
  
  
  #### DE genes ####
  
  #write.csv(res_ordered, file= cts_data_genes)
  
  
  genes <- read.csv(cts_data_genes, header = TRUE, stringsAsFactors = FALSE)
  head(genes)
  
  features <-read.table("C:/Users/Joseph McGirr Lab/Desktop/Cyprinodon/my.gtf.geneid", na.strings=c("", "NA"), sep="\t", header = FALSE, stringsAsFactors = FALSE)
  features["Chr"] <- features$V1
  features["Start"] <- features$V4
  features["End"] <- features$V5
  features["feature_gtf"] <- features$V3
  features <- features[which(features$feature_gtf == "CDS"),]
  features <- cbind(features, colsplit(features$V9, "gene_name ", c("junk", "symboljunk")))
  features <- cbind(features, colsplit(features$symbol, ";", c("symbol", "morejunk")))
  features <- cbind(features, colsplit(features$junk, "gene_id ", c("juuunk", "geneidjunk")))
  features <- cbind(features, colsplit(features$geneidjunk, ";", c("Geneid", "jnk")))
  keeps <- c("Chr","symbol", "Geneid")
  features <- features[keeps]
  features["tag"] <- paste(features$symbol, features$Geneid, sep = ":")
  features <- subset(features, !duplicated(features[,4])) 
  features_gff <- read.csv("D:/Martin Lab/RNA-seq/feature_table.csv", header = TRUE, stringsAsFactors = FALSE)
  features_gff["Chr"] <- features_gff$genomic_accession
  keeps <- c("Chr", "symbol", "related_accession", "product_accession", "name")
  features_gff <- features_gff[which(features_gff$feature == "CDS"),]
  features_gff["tag_f"] <- paste(features_gff$symbol, features_gff$Chr)
  features_gff <- subset(features_gff, !duplicated(features_gff[22]))
  features_gff <- features_gff[keeps]
  final_features <- merge(features, features_gff, by = c("Chr", "symbol"))
  nrow(final_features) 
  head(final_features)
  final <- merge(genes, final_features, by = c("Geneid"))
  final <- final[order(final$padj, decreasing = FALSE),]
  head(final)
  nrow(final)
  length(unique(final$Geneid))
  length(unique(final$Geneid)) / length(unique(res_ordered$Geneid))
  #write.csv(final, genes_out, row.names = FALSE)
  
}




#####
