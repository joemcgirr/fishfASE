master <- read.table("D:/Martin Lab/rna_2018/table_maker_master.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
head(master)
all <- read.table("C:/Users/jmcgirr/Desktop/all_samples.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
head(all)
final <- merge(master, all, by = c("sample"))
nrow(final)
head(final)
setwd("C:/Users/jmcgirr/Desktop/")
write.table(final, "table_maker_master.txt", row.names = FALSE, quote= FALSE,sep="\t")





cts <- counts(dds, normalized=TRUE)
head(cts)
counts(dds)
head(cts)
test <- as.data.frame(cts)
head(test)
test$Geneid <- rownames(test)
test2 <- test[,c(which(colnames(test)=="Geneid"),which(colnames(test)!="Geneid"))]
head(test2)
write.table(ase_8, "C:/Users/jmcgirr/Desktop/test_table.txt", row.names = FALSE, quote= FALSE,sep="\t")
cts <- read.table("C:/Users/jmcgirr/Desktop/test_table.txt" ,sep = "\t",header = TRUE)
head(cts)



keeps <- c("Geneid", datTraits$sample)
cts <- cts[keeps]
write.table(cts, "C:/Users/jmcgirr/Desktop/test_table_wgcna.txt", row.names = FALSE, quote= FALSE,sep="\t")
cts <- as.matrix(read.table("C:/Users/jmcgirr/Desktop/test_table_wgcna.txt" ,sep = "\t",header = TRUE,row.names=1))
head(cts)

hits <- read.table("C:/Users/jmcgirr/Desktop/cyprinodon_to_danio_one_way_best_hit.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
f <- read.table("C:/Users/jmcgirr/Desktop/features_gff.txt", header = FALSE, stringsAsFactors = FALSE, sep = "\t")
head(f)
head(hits)

gff <- read.table("D:/Cyprinodon/gff3_to.bed", header = FALSE, stringsAsFactors = FALSE, sep = "\t")
head(gff)
genes <- gff[which(gff$V8 == 'CDS'),]
nrow(genes)


old <- read.csv("E:/embryo_DE_reciprocal_hits_delete.csv")
head(old)
nrow(old)
head(two_way)
two_way$protein_accession <- two_way$product_accession
new <- merge(old, two_way,all.x = TRUE, by = c("protein_accession"))
nrow(new)
tail(new)
write.csv(new, "E:/embryo_DE_reciprocal_hits_blast_fixed.csv", row.names = FALSE, quote= FALSE)
nrow(new)
head(new)
#prop genes with one way hit
nrow(new[which(new$zeb_gene_symbol_one_way != "NA"),]) / nrow(new)
#prop genes with two way hit
nrow(new[which(new$zeb_gene_symbol_two_way != "NA"),]) / nrow(new)


cols_plot_counts <- c("blue4","blue4","blue4","blue4", "red","red","red","green4","green4","green4", "red","red","red","green4","green4","green4")
#8
cols_plot_counts <- c("red","red","red","green4","green4","green4", "red","red","red","green4","green4","green4","blue4","blue4","blue4","red","red","red","green4","green4","green4")

plotCounts(dds.sub, gene="gene15697", intgroup= "species", pch = 16, main = "",col = cols_plot_counts)

de <- read.csv("D:/Martin Lab/RNA-seq/axm/post_reviews/GO_enrichment/DE_(CRPA1_&_LILA_&_CRPM1_&_LILM)_vs_(LAxLM)_17dpf_genes_sig.csv", header = TRUE, stringsAsFactors = FALSE)
de <- read.csv("D:/Martin Lab/RNA-seq/axm/post_reviews/conditions/DE_(CRPA2_&_OSPA)_vs_(CRPM2_&_OSPM)_8dpf_genes.csv", header = TRUE, stringsAsFactors = FALSE)

go <- read.table("D:/Martin Lab/RNA-seq/axm/post_reviews/GO_enrichment/a_vs_m_8dpf_skeletal_morphogenesis.txt", header = TRUE, stringsAsFactors = FALSE)
go <- read.table("D:/Martin Lab/RNA-seq/axm/post_reviews/GO_enrichment/axm_mis_17_embryonic_skeletal_morphogenesis.txt", header = TRUE, stringsAsFactors = FALSE)

head(de)
head(go)
final <- merge(go, de, by = c("symbol"))
final
write.table(final,"D:/Martin Lab/RNA-seq/axm/post_reviews/GO_enrichment/axm_mis_17_embryonic_skeletal_morphogenesis_table.txt",quote = FALSE, row.names = FALSE, sep = "\t")
write.table(final,"D:/Martin Lab/RNA-seq/axm/post_reviews/GO_enrichment/a_vs_m_8dpf_skeletal_morphogenesis_table.txt",quote = FALSE, row.names = FALSE, sep = "\t")


i <- 3
subset_comps <- c(1,2,3)
subset_comps <- c(1,2,4,15)

for (i in subset_comps)
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
  sample_sizes_plot <- paste(sample_size_a, " parents vs. ", sample_size_b," hybrids", sep = "")
  
  total_genes_plot <- paste(total_genes, "transcripts", sep = " ")
  de_total_plot <- paste(de_total, "DE", sep = " ")
  prop_de <- paste((100*(round(prop_de, digits = 3))), "% DE", sep = "")
  plotMA(resLFC, ylim=c(-6,5), main = comp_file)
  #legend("bottomright", legend=c(sample_sizes_plot,total_genes_plot, de_total_plot, prop_de),cex=1.0, bty = 'n')
  #legend("bottomleft", legend=c(paste((100*(round(de_dn, digits = 3))), "% DE down", sep = "")),cex=0.8, bty = 'n')
  #legend("topleft", legend=c(paste((100*(round(de_up, digits = 3))), "% DE up", sep = "")),cex=0.8, bty = 'n')
  
  #head(resLFC)
  de_plot <- paste("D:/Martin Lab/RNA-seq/axm/post_reviews/figures/",comp_file, "_de_plot.tiff", sep = "")
  tiff(de_plot, width = 5, height = 5, units = 'in', res = 1000)
  plotMA(resLFC, ylim=c(-5,4), main = comp_file, xaxt= "n", cex.axis =1.5)
  axis(1, at=c(10,100,1000,100000), cex.axis = 1.5)
  #legend("bottomright", legend=c(total_genes_plot, de_total_plot, prop_de),cex=0.8, bty = 'n')
  #legend("bottomleft", legend=c(paste((100*(round(de_dn, digits = 3))), "% DE down", sep = "")),cex=0.8, bty = 'n')
  #legend("topleft", legend=c(paste((100*(round(de_up, digits = 3))), "% DE up", sep = "")),cex=0.8, bty = 'n')
  dev.off()
  
  #### overlap with genes ###
  
  #write.csv(res_ordered, file= cts_data_genes)
  genes <- read.csv(cts_data_genes, header = TRUE, stringsAsFactors = FALSE)
  final <- merge(genes, final_features, by = c("Geneid"))
  length(unique(final$Geneid))
  length(unique(final$Geneid)) / length(unique(res_ordered$Geneid))
  head(final)
  nrow(final)
  head(blast_key)
  final <- merge(final, blast_key,all.x = TRUE, by = c("product_accession"))
  final <- final[order(final$padj, decreasing = FALSE),]
  #write.csv(final, genes_out, row.names = FALSE)
  
}



library(ggplot2)

plt <- ggplot(mis8, aes(x=dxy, y=prop_DE, color = cross_type)) + geom_point(aes(shape=lake),size = 4)
plt +   theme_classic() +scale_color_manual(values = cols) +labs(x = "\ntranscriptome-wide mean Fst", y = "% genes misexpressed in F1 hybrids\n") +
  guides(col = guide_legend(override.aes = list(shape = 19, size = 5)))+ggtitle("8 dpf") +
  theme(plot.title = element_text(hjust = 0.5))#+
  geom_segment(aes(x = mis8$fst, y = mis8$prop_DE, xend = mis2$fst, yend = mis2$prop_DE),show.legend = FALSE)

plt <- ggplot(mis2, aes(x=dxy_dna, y=prop_DE_subset, color = cross_type))
plt + theme_classic() +scale_color_manual(values = cols) +labs(x = "\ngenome-wide mean dxy", y = "% genes misexpressed in F1 hybrids\n") +guides(col = guide_legend(override.aes = list(shape = 19, size = 5))) +ggtitle("48 hpf misexpressed contrasts")+theme(plot.title = element_text(hjust = 0.5))+
geom_segment(aes(x = mis2$dxy_dna, y = prop_DE_subset, xend = mis2$dxy_dna, yend = mis2$prop_misexpressed_2_way_95),show.legend = FALSE)

plt <- ggplot(mis8, aes(x=dxy_dna, y=prop_DE_subset, color = cross_type)) + geom_point(aes(shape=lake),size = 4)
plt + theme_classic() +scale_color_manual(values = cols) +labs(x = "\ngenome-wide mean dxy", y = "% genes misexpressed in F1 hybrids\n") +guides(col = guide_legend(override.aes = list(shape = 19, size = 5))) +ggtitle("8 dpf misexpressed contrasts vs 2 way 95 CI overlap")+theme(plot.title = element_text(hjust = 0.5))+
geom_segment(aes(x = mis8$dxy_dna, y = mis8$prop_DE_subset, xend = mis8$dxy_dna, yend = mis8$prop_misexpressed_2_way_95),show.legend = FALSE)


plt <- ggplot(mis2, aes(x=dxy_dna, y=prop_DE_subset, color = cross_type)) + geom_point(aes(shape=lake),size = 4)
plt + theme_classic() +scale_color_manual(values = cols) +labs(x = "\ngenome-wide mean dxy", y = "% genes misexpressed in F1 hybrids\n") +guides(col = guide_legend(override.aes = list(shape = 19, size = 5))) +ggtitle("48 hpf misexpressed contrasts vs 2 way 95 CI overlap")+theme(plot.title = element_text(hjust = 0.5))+
  geom_segment(aes(x = mis2$dxy_dna, y = mis2$prop_DE_subset, xend = mis2$dxy_dna, yend = mis2$prop_misexpressed_2_way_95),show.legend = FALSE)



go <- read.table("C:/Users/jmcgirr/Desktop/Book1.txt", header = TRUE, stringsAsFactors = FALSE)
do <- read.table("C:/Users/jmcgirr/Desktop/Book2.txt", header = TRUE, stringsAsFactors = FALSE)
merge(go,do,by = c("rna_sample"))










setwd("D:/Martin Lab/rna_2018/all_2018_samples/")
#comp <- "all_samples_8dpf"
comp <- "all_samples_48hpf"
#comp <- "table_maker_master"
comp_file <- paste(comp , ".txt", sep = "")
#comp_file <- "table_maker_master_outlier_rm.txt"
all_cts <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/ase_data/all_samples_2018_counts_mrna_saf.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
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


dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = colData,
                              design= ~f1)


dds <- estimateSizeFactors(dds)
idx <- rowSums(counts(dds, normalized=TRUE) >= 1 ) >= 1
length(rowSums(counts(dds, normalized=TRUE) >= 1 )) >= 1
length(idx)
dds <- dds[idx,]
dds <- DESeq(dds)
norm_cts <- data.frame(counts(dds, normalized=TRUE))
head(norm_cts)
nrow(norm_cts)
norm_cts$Geneid <- rownames(norm_cts)
norm_cts_48 <- norm_cts

#write.table(norm_cts_48, "C:/Users/jmcgirr/Documents/all_2018_samples/norm_cts_48hpf_no_seq_round.txt", row.names = FALSE, quote= FALSE,sep="\t")
#norm_cts_48 <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/norm_cts_48hpf_no_seq_round.txt", header=TRUE, stringsAsFactors = FALSE)
#norm_cts_8 <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/norm_cts_8dpf.txt", header=TRUE, stringsAsFactors = FALSE)

#head(norm_cts_48)


setwd("D:/Martin Lab/rna_2018/all_2018_samples/")
comp <- "all_samples_8dpf"
#comp <- "all_samples_48hpf"
#comp <- "table_maker_master"
comp_file <- paste(comp , ".txt", sep = "")
#comp_file <- "table_maker_master_outlier_rm.txt"
all_cts <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/ase_data/all_samples_2018_counts_mrna_saf.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
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


dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = colData,
                              design= ~f1)
#}

dds <- estimateSizeFactors(dds)
idx <- rowSums(counts(dds, normalized=TRUE) >= 1 ) >= 1
length(rowSums(counts(dds, normalized=TRUE) >= 1 )) >= 1
length(idx)
dds <- dds[idx,]
dds <- DESeq(dds)
norm_cts <- data.frame(counts(dds, normalized=TRUE))
head(norm_cts)
nrow(norm_cts)
norm_cts$Geneid <- rownames(norm_cts)
norm_cts_8 <- norm_cts

write.table(norm_cts_48, "C:/Users/jmcgirr/Documents/all_2018_samples/norm_cts_48hpf.txt", row.names = FALSE, quote= FALSE,sep="\t")
write.table(norm_cts_8, "C:/Users/jmcgirr/Documents/all_2018_samples/norm_cts_8dpf.txt", row.names = FALSE, quote= FALSE,sep="\t")

norm_cts_48t <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/norm_cts_48hpf.txt", header=TRUE, stringsAsFactors = FALSE)
norm_cts_8t <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/norm_cts_8dpf.txt", header=TRUE, stringsAsFactors = FALSE)

head(norm_cts_48t)
df <- merge(norm_cts_48t,norm_cts_8t, all = TRUE,by = c("Geneid"))
head(df)
nrow(df)
write.table(df, "C:/Users/jmcgirr/Documents/all_2018_samples/norm_cts_all_2018_no_seq_round_control.txt", row.names = FALSE, quote= FALSE,sep="\t")

a <- c("OMT1","OMT2","OMT3","OPE1","OPE2","OPE3","OPE4","OPE5","OPT1","OPT2","OPT3","ORE1","ORE2","ORE3","ORE4","ORE5","ORT1","ORT2","ORT3","OUE1","OUE3","OUE4","OUT1","OUT2","OUT3","OVE1","OVE4","OVE5","OVT1","OVT2","OVT3","OXE2","OXT1","OXT2","OXT3","OYE1","OYE2","OYE3","OYE4","OYE5","OYT1","OYT2","OYT3","OZE2","OZE4","OZE5","OZT1","OZT2","OZT3","PAE1","PAE2","PAE5","PAT1","PAT2","PAT3","CQT1","CQT2","CUT1","CUT2","CUT3","CVE1","CVE2","CVE5","CVT1","CVT2","CVT3","CWE2","CWE3","CWE4","CWT1","CWT2","CWT3","CXE2","CXE3","CXE4","CXT1","CXT2","CXT3","LFE2","LFE3","LFE4","LFT1","LFT2","LFT3","LGE3","LGE4","LGE5","LGT1","LGT2","LGT3","LIE2","LIE3","LIE5","LIT1","LIT2","LIT3","LKE1","LKE2","LKE3","LKT1","LKT2","LKT3","LPE1","LPE2","LPE3","LPE4","LPE5","NAE1","NAE2","NAE4","NAT1","NAT2","NAT3","NCE1","NCE2","NCE3","NCE4","NCE5","NCT1","NCT2","NCT3","OAE1","OAE2","OAE3","OAE4","OAE5","OAT1","OAT2","OAT3","OME1","OME2","OME3","OME4","OME5","CAE1","CAE2","CAE3","CAE4","CAE5","CAT1","CAT2","CAT3","CME1","CME2","CME5","CMT1","CMT2","CMT3","CPE1","CPE2","CPE3","CPE4","CPE5","CPT1","CPT2","CPT3","CPU1","CPU3","CPU5","CQE1","CQE2","CQE3")
b <- c("NAE4","LGT3","PAT1","PAE5","PAE2","PAE1","OZT3","OZT2","OZT1","OZE5","OZE4","OZE2","OYT3","OYT2","OYT1","OYE5","OYE4","OYE2","OXT2","OYE1","OXT3","OXT1","OVT3","OVT2","OXE2","OVE1","OVE4","OUT2","OVE5","OUT3","OUT1","OUE4","OUE3","ORT3","OUE1","ORT2","ORE2","ORE1","OPT3","OPE2","OPE1","OMT3","OPE3","OPE5","OPE4","OPT2","OPT1","ORE5","ORT1","ORE3","ORE4","OME4","NAE2","NCE2","CAE2","NAT1","NAT2","LPE4","LPE5","NAE1","OME5","LPE3","OAE5","LPE2","NCE5","NCE3","OMT1","OMT2","NCE4","NCT3","NCT1","NCT2","OAE1","OAE3","NCE1","NAT3","OAT1","OAE2","OAE4","OAT3","OAT2","OME2","OME3","OME1","CAE3","LPE1","CAT2","LKE2","LIT3","LKT1","LIT2","LKE1","LKT3","LKT2","LKE3","CPE4","CAT3","CQE3","CPE5","CPT2","CPT1","CXT3","CAT1","CAE5","CAE4","CUT2","CME5","CVT1","CVE5","LGE5","CQE2","CUT1","CQT2","CQT1","LFT1","LFT3","CPU1","CPT3","CMT3","LIT1","LIE3","LGT1","LIE2","LIE5","LGT2","LFT2","LGE4","LGE3","CXT2","LFE4","LFE2","LFE3","CVT2","CXE4","CPE2","CPE1","CMT1","CMT2","CPE3","CME1","CME2","CPU5","CPU3","CQE1","CWE2","CVT3","CVE1","CVE2","CUT3","CWE3","CWE4","CWT1","CWT2","CXT1","CXE3","CXE2","CWT3","CAE1")
setdiff(a,b)
b <- c("NAE4","OXT1","OUT1","OZT3","PAT2","PAE1","OZT2","LGT3","OUT3","OZT1","OUT2","LKE3","OXT3","OZE5","PAE2","OYT3","OVE4","OUE4","OYT2","OYT1","OAT3","OXT2","NAT2","OZE4","OAT1","OVT3","NAE1","OUE3","OYE4","CXE4","NAE2","OXE2","OZE2","OVT1","LGT1","PAT3","OVE5","OYE5","OYE1","CUT1","OYE3","OAT2","CUT3","CAE5","LGT2","OMT2","OYE2","OMT3","NAT3","OVT2","CVE1","CWT1","CXE2","NAT1","PAE5","OME2","OPT1","OVE1","LIT2","CAT3","CAT2","OAE1","OPE5","OPT3","CUT2","OME5","OPT2","CMT1","LGE3","OMT1","OME3","OAE3","CVT3","LFT1","OME4","CWE3","CAE1","LIT3","OPE4","LFE4","LFT3","CAT1","CXT3","LFE3","CXT2","CXT1","OAE2","CWT2","LIE5","LKT2","CQT2","LKT3","CAE4","LFT2","OAE4","LKT1","OME1","LIT1","CWE2","CWE4","CPT3","CME1","OPE3","LKE1","CAE3","CXE3","CWT3","CQE2","CME2","CQT1","OPE1","CMT3","CQE1","OPE2","CVE2","CVE5","LGE4","CPU1","CME5","CMT2","LFE2","LIE3","CQE3","CVT1","CAE2","CPU3","CPU5","CPT2","CPT1","CPE1","LIE2","CPE2","NCT1","NCT3","CPE5","LKE2","NCT2","NCE2","CPE4","NCE3","CVT2","CPE3","NCE1","NCE4","NCE5")
b <- c("PAT2","PAE5","PAT3","PAE1","OZT3","OZT2","OZT1","OZE5","OZE4","OZE2","OYT3","OYT1","OYE5","OYE4","OYE3","OYE2","OYE1","OXT3","OXT2","OXT1","OXE2","OVT3","OVT2","OVT1","OVE4","OVE5","OVE1","OUT3","OUT1","OUE4","OUE3","ORT3","ORT1","ORE5","ORE1","OPT3","OPT2","OPE5","OPE4","OPE3","OPE2","OPE1","OMT3","OMT2","OMT1","OME5","OME4","OME3","OME2","OME1","OAT3","OAT2","OAT1","OAE5","OAE3","OAE2","NCT3","NCT2","NCT1","NCE5","NCE3","NCE2","NAT3","NAT2","NAT1","NAE4","NAE2","LPE5","LPE4","LPE3","LPE2","LPE1","LKT3","LKT2","OAE1","NAE1","LKT1","OPT1","NCE1","LKE3","LKE2","LKE1","PAE2","LIT3","OYT2","LIT2","LIT1","LIE5","LIE3","LIE2","LGT3","LGT2","LGT1","LGE5","LGE4","LGE3","OUT2","LFT3","LFT2","LFT1","LFE4","LFE3","LFE2","CXT3","CXT2","CXT1","CXE3","CXE4","CXE2","CWT3","CWT2","CWT1","CWE4","OAE4","CWE2","CVT3","CVT2","CVT1","CVE5","CVE2","CUT1","CQT2","CQT1","CQE3","CQE2","CQE1","CPU5","CPU1","CPT3","CPT2","CPT1","CPE5","CPE4","CPE3","CPE2","CPE1","CMT3","CMT2","CMT1","CME5","CME2","CAT3","CAT2","CAT1","CAE4","CAE3","CAE5","CAE2","NCE4","CVE1","CUT3","CWE3","CUT2","CPU3","CME1","CAE1")
b <- c("OZT2","PAE1","OXT3","OUT1","OUT2","OZT3","OUT3","OXT2","OYT1","OVE4","OZT1","OYT2","OXT1","OYT3","LKE3","NAE1","OZE5","NAE2","PAE2","PAT1","CUT2","OVT3","OZE4","LGT1","LFT1","NAT2","LGT2","CAT2","OMT3","OAT3","OUE4","OAT1","CUT1","OZE2","LGE3","OYE4","OXE2","OVE5","CAT1","OYE1","CUT3","LIT2","CXE2","CAT3","CXT2","CVT3","LFT3","CMT3","NAT1","LFT2","CVE1","LIT3","CXT3","CAE5","CXE4","OUE3","OUE1","OVE1","OYE2","OAT2","CAE3","CWE3","OVT2","OYE5","PAE5","CXT1","CWT1","NAT3","LFE4","LKE1","LIE5","OMT2","LFE3","CMT1","OAE1","LKT3","LKT2","CQE1","CQT2","OME4","CVT1","CME5","LKT1","CAE4","CPT3","CVE5","CMT2","OPE5","CXE3","CWE2","CQE3","LIT1","OPT3","CPU1","OME2","CQE2","OPT1","LGE4","CQT1","OMT1","OME5","CAE2","CME2","LFE2","CWT3","CWE4","OAE3","CWT2","OME3","CME1","OPT2","CVE2","LIE3","OPE1","OAE2","OPE2","CPU3","OAE4","OPE4","OME1","CPU5","OPE3","CPE1","CPT1","LIE2","CPT2","CVT2","NCT1","NCE2","CPE2","NCT3","CPE4","NCE3","CPE5","LKE2","NCT2","NCE1","NCE4","CPE3","NCE5","CAE1")
length(b)
b <- c("PAT1","OUE1","NAE4","OXT1","OUT1","OZT3","PAT2","PAE1","OZT2","LGT3","OUT3","OZT1","OUT2","LKE3","OXT3","OZE5","PAE2","OYT3","OVE4","OUE4","OYT2","OYT1","OAT3","OXT2","NAT2","OZE4","OAT1","OVT3","NAE1","OUE3","OYE4","CXE4","NAE2","OXE2","OZE2","OVT1","LGT1","PAT3","OVE5","OYE5","OYE1","CUT1","OYE3","OAT2","CUT3","CAE5","LGT2","OMT2","OYE2","OMT3","NAT3","OVT2","CVE1","CWT1","CXE2","NAT1","PAE5","OME2","OPT1","OVE1","LIT2","CAT3","CAT2","OAE1","OPE5","OPT3","CUT2","OME5","OPT2","CMT1","LGE3","OMT1","OME3","OAE3","CVT3","LFT1","OME4","CWE3","CAE1","LIT3","OPE4","LFE4","LFT3","CAT1","CXT3","LFE3","CXT2","CXT1","OAE2","CWT2","LIE5","LKT2","CQT2","LKT3","CAE4","LFT2","OAE4","LKT1","OME1","LIT1","CWE2","CWE4","CPT3","CME1","OPE3","LKE1","CAE3","CXE3","CWT3","CQE2","CME2","CQT1","OPE1","CMT3","CQE1","OPE2","CVE2","CVE5","LGE4","CPU1","CME5","CMT2","LFE2","LIE3","CQE3","CVT1","CAE2","CPU3","CPU5","CPT2","CPT1","CPE1","LIE2","CPE2","NCT1","NCT3","CPE5","LKE2","NCT2","NCE2","CPE4","NCE3","CVT2","CPE3","NCE1","NCE4","NCE5")
a <- c("CAE1","CAE2","CAE3","CAE4","CAE5","CAT1","CAT2","CAT3","CME1","CME2","CME5","CMT1","CMT2","CMT3","CPE1","CPE2","CPE3","CPE4","CPE5","CPT1","CPT2","CPT3","CPU1","CPU3","CPU5","CQE1","CQE2","CQE3","CQT1","CQT2","CUT1","CUT2","CUT3","CVE1","CVE2","CVE5","CVT1","CVT2","CVT3","CWE2","CWE3","CWE4","CWT1","CWT2","CWT3","CXE2","CXE3","CXE4","CXT1","CXT2","CXT3","LFE2","LFE3","LFE4","LFT1","LFT2","LFT3","LGE3","LGE4","LGT1","LGT2","LGT3","LIE2","LIE3","LIE5","LIT1","LIT2","LIT3","LKE1","LKE2","LKE3","LKT1","LKT2","LKT3","NAE1","NAE2","NAE4","NAT1","NAT2","NAT3","NCE1","NCE2","NCE3","NCE4","NCE5","NCT1","NCT2","NCT3","OAE1","OAE2","OAE3","OAE4","OAT1","OAT2","OAT3","OME1","OME2","OME3","OME4","OME5","OMT1","OMT2","OMT3","OPE1","OPE2","OPE3","OPE4","OPE5","OPT1","OPT2","OPT3","OUE3","OUE4","OUT1","OUT2","OUT3","OVE1","OVE4","OVE5","OVT1","OVT2","OVT3","OXE2","OXT1","OXT2","OXT3","OYE1","OYE2","OYE3","OYE4","OYE5","OYT1","OYT2","OYT3","OZE2","OZE4","OZE5","OZT1","OZT2","OZT3","PAE1","PAE2","PAE5","PAT2","PAT3")
setdiff(b,a)


a <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/sweed/lab_desktop_results/sweed_key.txt", stringsAsFactors = FALSE)
head(a)
a$scaffold <- a$V2
b <- read.table("D:/Cyprinodon/scaffold_conversion.txt", header = TRUE, stringsAsFactors = FALSE)
head(b)
c <- merge(a,b,by = c("scaffold"))
head(c)
c <- c[order(c$V1, decreasing = FALSE),]
head(c)
c <- c[c("V1","chr")]
write.table(c,"C:/Users/jmcgirr/Documents/all_2018_samples/sweed/lab_desktop_results/sweed_key.txt", quote = FALSE, row.names = FALSE, sep = "\t")
