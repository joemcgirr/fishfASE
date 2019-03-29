
library(DESeq2)
library(reshape2)
library(seqinr)
library(plyr)
library(MASS)
library(AER)
library(rlang)
library(venn)




# hybrid expression significantly different from all caribbean samples

final_features <- read.table("D:/Martin Lab/rna_2018/all_2018_samples/features_gff.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t", quote = "")
# mrna counts
all_cts <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/ase_data/all_samples_2018_counts_mrna_saf.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")

blast_key <- read.table("D:/Martin Lab/rna_2018/all_2018_samples/blast/cyprinodon_to_danio_one_way_best_hit_symbols.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
head(blast_key)
#i <- 84
#subset_comps <- c(31:34)
hybs <- c("CAxCM",
          "CMxCA",
          "OAxOM",
          "OMxOA",
          "CAxCP",
          "OAxOP",
          "OPxOA",
          "CMxCP",
          "OPxOM")
stage <- "48hpf"
pures <- c("NCA",
             "OSPA",
             "OSPM",
             "OSPP",
             "CRPA",
             "CRPM",
             "CRPP",
             "UPxUA")

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
for (i in c(1:length(hybs)))
{
  pops1_name <- hybs[i]
  pops2_name <- "all_pure"

  comp_name <- paste("condition_species_",pops1_name, "_vs_", pops2_name, "_" ,stage,".txt", sep = "")
  single_run <- paste(pops1_name, "_vs_", pops2_name, "_" ,stage, sep = "")
  
  setwd("C:/Users/jmcgirr/Documents/all_2018_samples/super_misexpression/conditions/")
  master <- read.table("D:/Martin Lab/rna_2018/all_2018_samples/table_maker_master_outlier_rm.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
  head(master)
  pops1_table <- master[master$f1 %in% pops1_name, ]
  pops2_table <- master[master$f1 %in% pures, ]
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
  cts_data <-             "DESeq_counts_delete.txt"
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
  

  i_stats <- data.frame(pop1=pops1_name, 
                        pop2=pops2_name,
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
  de_plot <- paste("C:/Users/jmcgirr/Documents/all_2018_samples/super_misexpression/de_plots/",comp_file, "_de_plot.tiff", sep = "")
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
#write.table(sum_stats, "C:/Users/jmcgirr/Documents/all_2018_samples/super_misexpression/summary_stats_48hpf.txt", row.names = FALSE, quote= FALSE,sep="\t")


#load super misexpressed
{
  DE_genes_dir <- "C:/Users/jmcgirr/Documents/all_2018_samples/super_misexpression/conditions/"
  stage <- "8dpf"
  # load misexpressed
  {
    # am crosses
    caxcm <- read.csv(paste(DE_genes_dir,"DE_CAxCM_vs_all_pure_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    cmxca <- read.csv(paste(DE_genes_dir,"DE_CMxCA_vs_all_pure_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    oaxom <- read.csv(paste(DE_genes_dir,"DE_OAxOM_vs_all_pure_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    omxoa <- read.csv(paste(DE_genes_dir,"DE_OMxOA_vs_all_pure_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    # ap crosses
    caxcp <- read.csv(paste(DE_genes_dir,"DE_CAxCP_vs_all_pure_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    oaxop <- read.csv(paste(DE_genes_dir,"DE_OAxOP_vs_all_pure_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    opxoa <- read.csv(paste(DE_genes_dir,"DE_OPxOA_vs_all_pure_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    # specialist crosses
    cmxcp <- read.csv(paste(DE_genes_dir,"DE_CMxCP_vs_all_pure_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    opxom <- read.csv(paste(DE_genes_dir,"DE_OPxOM_vs_all_pure_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    
    # sig super misexpressed
    caxcm <- caxcm[which(caxcm$padj <= 0.05),] 
    cmxca <- cmxca[which(cmxca$padj <= 0.05),] 
    oaxom <- oaxom[which(oaxom$padj <= 0.05),] 
    omxoa <- omxoa[which(omxoa$padj <= 0.05),] 
    caxcp <- caxcp[which(caxcp$padj <= 0.05),] 
    oaxop <- oaxop[which(oaxop$padj <= 0.05),] 
    opxoa <- opxoa[which(opxoa$padj <= 0.05),] 
    cmxcp <- cmxcp[which(cmxcp$padj <= 0.05),] 
    opxom <- opxom[which(opxom$padj <= 0.05),] 

    
    caxcm_8_SME <- caxcm$related_accession
    cmxca_8_SME <- cmxca$related_accession
    oaxom_8_SME <- oaxom$related_accession
    omxoa_8_SME <- omxoa$related_accession
    caxcp_8_SME <- caxcp$related_accession
    oaxop_8_SME <- oaxop$related_accession
    opxoa_8_SME <- opxoa$related_accession
    cmxcp_8_SME <- cmxcp$related_accession
    opxom_8_SME <- opxom$related_accession

  }
  stage <- "48hpf"
  # load misexpressed
  {
    # am crosses
    caxcm <- read.csv(paste(DE_genes_dir,"DE_CAxCM_vs_all_pure_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    cmxca <- read.csv(paste(DE_genes_dir,"DE_CMxCA_vs_all_pure_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    oaxom <- read.csv(paste(DE_genes_dir,"DE_OAxOM_vs_all_pure_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    omxoa <- read.csv(paste(DE_genes_dir,"DE_OMxOA_vs_all_pure_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    # ap crosses
    caxcp <- read.csv(paste(DE_genes_dir,"DE_CAxCP_vs_all_pure_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    oaxop <- read.csv(paste(DE_genes_dir,"DE_OAxOP_vs_all_pure_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    opxoa <- read.csv(paste(DE_genes_dir,"DE_OPxOA_vs_all_pure_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    # specialist crosses
    cmxcp <- read.csv(paste(DE_genes_dir,"DE_CMxCP_vs_all_pure_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    opxom <- read.csv(paste(DE_genes_dir,"DE_OPxOM_vs_all_pure_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    
    # sig super misexpressed
    caxcm <- caxcm[which(caxcm$padj <= 0.05),] 
    cmxca <- cmxca[which(cmxca$padj <= 0.05),] 
    oaxom <- oaxom[which(oaxom$padj <= 0.05),] 
    omxoa <- omxoa[which(omxoa$padj <= 0.05),] 
    caxcp <- caxcp[which(caxcp$padj <= 0.05),] 
    oaxop <- oaxop[which(oaxop$padj <= 0.05),] 
    opxoa <- opxoa[which(opxoa$padj <= 0.05),] 
    cmxcp <- cmxcp[which(cmxcp$padj <= 0.05),] 
    opxom <- opxom[which(opxom$padj <= 0.05),] 
    
    
    caxcm_48_SME <- caxcm$related_accession
    cmxca_48_SME <- cmxca$related_accession
    oaxom_48_SME <- oaxom$related_accession
    omxoa_48_SME <- omxoa$related_accession
    caxcp_48_SME <- caxcp$related_accession
    oaxop_48_SME <- oaxop$related_accession
    opxoa_48_SME <- opxoa$related_accession
    cmxcp_48_SME <- cmxcp$related_accession
    opxom_48_SME <- opxom$related_accession
    
  }
}


venn(list(caxcm_48_SME=caxcm_48_SME,cmxca_48_SME=cmxca_48_SME), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
venn(list(caxcm_8_SME=caxcm_8_SME,cmxca_8_SME=cmxca_8_SME), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
venn(list(oaxom_48_SME=oaxom_48_SME,omxoa_48_SME=omxoa_48_SME), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
venn(list(oaxom_8_SME=oaxom_8_SME,omxoa_8_SME=omxoa_8_SME), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)

venn(list(caxcp_48_SME=caxcp_48_SME), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
venn(list(caxcp_8_SME=caxcp_8_SME), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
venn(list(oaxop_48_SME=oaxop_48_SME,opxoa_48_SME=opxoa_48_SME), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
venn(list(oaxop_8_SME=oaxop_8_SME,opxoa_8_SME=opxoa_8_SME), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)

venn(list(cmxcp_48_SME=cmxcp_48_SME), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
venn(list(cmxcp_8_SME=cmxcp_8_SME), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
venn(list(opxom_48_SME=opxom_48_SME), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
venn(list(omxom_8_SME=opxom_8_SME), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)

venn(list(caxcm_48_SME=caxcm_48_SME,cmxca_48_SME=cmxca_48_SME,oaxom_48_SME=oaxom_48_SME,omxoa_48_SME=omxoa_48_SME), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
venn(list(caxcm_8_SME=caxcm_8_SME,cmxca_8_SME=cmxca_8_SME,oaxom_8_SME=oaxom_8_SME,omxoa_8_SME=omxoa_8_SME), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)

venn(list(caxcp_48_SME=caxcp_48_SME,oaxop_48_SME=oaxop_48_SME,opxoa_48_SME=opxoa_48_SME), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
venn(list(caxcp_8_SME=caxcp_8_SME,oaxop_8_SME=oaxop_8_SME,opxoa_8_SME=opxoa_8_SME), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)

venn(list(cmxcp_48_SME=cmxcp_48_SME,opxom_48_SME=opxom_48_SME), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
venn(list(cmxcp_8_SME=cmxcp_8_SME,opxom_8_SME=opxom_8_SME), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)

