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

##### DE between crpa samples
# mrna counts
all_cts <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/ase_data/all_samples_2018_counts_mrna_saf.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")



#pops1 <- "CRPA"
#48 hpf
#comps <- c("CAT1","CAT2","CAT3")
#8 dpf
#comps <- c("CAE1_CAE2","CAE1_CAE3","CAE1_CAE4","CAE1_CAE5","CAE2_CAE3","CAE2_CAE4","CAE2_CAE5","CAE3_CAE4","CAE3_CAE5","CAE4_CAE5")

#pops1 <- "CRPM"
#48 hpf
#comps <- c("CMT1","CMT2","CMT3")
#8 dpf
#comps <- c("CME1","CME2","CME5")

#pops1 <- "CRPP"
#48 hpf
#comps <- c("CPT1","CPT2","CPT3")
#8 dpf
#comps <- c("CPE1_CPE2","CPE1_CPE3","CPE1_CPE4","CPE1_CPE5","CPE2_CPE3","CPE2_CPE4","CPE2_CPE5","CPE3_CPE4","CPE3_CPE5","CPE4_CPE5")

#pops1 <- "UPxUA"
#48 hpf
#comps <- c("CQT1","CQT2")
#8 dpf
#comps <- c("CQE1","CQE2","CQE3")

#pops1 <- "NCA"
#48 hpf
#comps <- c("NCT1","NCT2","NCT3")
#8 dpf
#comps <- c("NCE1_NCE2","NCE1_NCE3","NCE1_NCE4","NCE1_NCE5","NCE2_NCE3","NCE2_NCE4","NCE2_NCE5","NCE3_NCE4","NCE3_NCE5","NCE4_NCE5")


#48 hpf
prop_des <- c()
for (comp in comps)
{
  #48 hpf
  #####
  stage <- "48hpf"
  #stage <- "8dpf"
  comp_name <- paste("condition_species_CRPA_vs_CRPA", "_" ,stage,".txt", sep = "")

  setwd("C:/Users/jmcgirr/Documents/all_2018_samples/conditions/")
  master <- read.table("D:/Martin Lab/rna_2018/all_2018_samples/table_maker_master_outlier_rm.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
  head(master)
  pops_table <- master[master$f1 %in% pops1, ]
  pops_table <- pops_table[which(pops_table$stage == stage),]
  pops1_table <- pops_table[which(pops_table$sample == comp),]
  pops2_table <- pops_table[which(pops_table$sample != comp),]
  
  pops1_table$species <- "a" 
  pops2_table$species <- "b"
  comp_table <- rbind(pops1_table,pops2_table)
  #####
  #8 dpf
  #####
  #stage <- "8dpf"
  #comp_name <- paste("condition_species_CRPA_vs_CRPA", "_" ,stage,".txt", sep = "")
  #
  #setwd("C:/Users/jmcgirr/Documents/all_2018_samples/conditions/")
  #master <- read.table("D:/Martin Lab/rna_2018/all_2018_samples/table_maker_master_outlier_rm.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
  #head(master)
  #pops1 <- "CRPA"
  #pops_table <- master[master$f1 %in% pops1, ]
  #pops_table <- pops_table[which(pops_table$stage == stage),]
  #
  #p1 <- strsplit(comp,"_")[[1]][1]
  #p2 <- strsplit(comp,"_")[[1]][2]
  #pops1_table <- pops_table[which(pops_table$sample == p1 | pops_table$sample == p2 ),]
  #pops2_table <- pops_table[which(pops_table$sample != p1 & pops_table$sample != p2 ),]
  #
  #pops1_table$species <- "a" 
  #pops2_table$species <- "b"
  #comp_table <- rbind(pops1_table,pops2_table)
  #####
  
  write.table(comp_table, "C:/Users/jmcgirr/Documents/all_2018_samples/conditions/crpa_v_crpa.txt", row.names = FALSE, quote= FALSE,sep="\t")
  
  sample_list <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/conditions/crpa_v_crpa.txt", header = TRUE, stringsAsFactors = FALSE)
  samples_a <- sample_list[which(sample_list$species == "a"),]
  samples_b <- sample_list[which(sample_list$species == "b"),]
  add_sequencing_round_to_model <- length(unique(samples_a$sequencing_round)) + length(unique(samples_b$sequencing_round))
  keeps <- c("Geneid", sample_list$sample)
  keeper <- sample_list$sample
  
  cts <- all_cts[keeps]
  #head(cts)
  write.table(cts, "C:/Users/jmcgirr/Documents/all_2018_samples/conditions/crpa_v_crpa_counts.txt", row.names = FALSE, quote= FALSE,sep="\t") 
  
  ### DIFFERENTIAL EXPRESSION ANALYSIS ###
  #vignette("DESeq2")
  cts <- as.matrix(read.table("C:/Users/jmcgirr/Documents/all_2018_samples/conditions/crpa_v_crpa_counts.txt" ,sep = "\t",header = TRUE,row.names=1))
  head(cts)
  nrow(cts)
  colData <- as.matrix(read.table("C:/Users/jmcgirr/Documents/all_2018_samples/conditions/crpa_v_crpa.txt" ,header = TRUE,row.names=1))
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
  prop_des <- c(prop_des, prop_de)

}  

mean(prop_des)
total_genes
range(prop_des)

#8 dpf
prop_des <- c()
for (comp in comps)
{
  #48 hpf
  #####
  #stage <- "48hpf"
  ##stage <- "8dpf"
  #comp_name <- paste("condition_species_CRPA_vs_CRPA", "_" ,stage,".txt", sep = "")
  #
  #setwd("C:/Users/jmcgirr/Documents/all_2018_samples/conditions/")
  #master <- read.table("D:/Martin Lab/rna_2018/all_2018_samples/table_maker_master_outlier_rm.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
  #head(master)
  #pops_table <- master[master$f1 %in% pops1, ]
  #pops_table <- pops_table[which(pops_table$stage == stage),]
  #pops1_table <- pops_table[which(pops_table$sample == comp),]
  #pops2_table <- pops_table[which(pops_table$sample != comp),]
  #
  #pops1_table$species <- "a" 
  #pops2_table$species <- "b"
  #comp_table <- rbind(pops1_table,pops2_table)
  ######
  ##8 dpf
  ######
  stage <- "8dpf"
  comp_name <- paste("condition_species_CRPA_vs_CRPA", "_" ,stage,".txt", sep = "")
  
  setwd("C:/Users/jmcgirr/Documents/all_2018_samples/conditions/")
  master <- read.table("D:/Martin Lab/rna_2018/all_2018_samples/table_maker_master_outlier_rm.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
  head(master)
  pops_table <- master[master$f1 %in% pops1, ]
  pops_table <- pops_table[which(pops_table$stage == stage),]
  
  p1 <- strsplit(comp,"_")[[1]][1]
  p2 <- strsplit(comp,"_")[[1]][2]
  pops1_table <- pops_table[which(pops_table$sample == p1 | pops_table$sample == p2 ),]
  pops2_table <- pops_table[which(pops_table$sample != p1 & pops_table$sample != p2 ),]
  
  pops1_table$species <- "a" 
  pops2_table$species <- "b"
  comp_table <- rbind(pops1_table,pops2_table)
  #####
  
  write.table(comp_table, "C:/Users/jmcgirr/Documents/all_2018_samples/conditions/crpa_v_crpa.txt", row.names = FALSE, quote= FALSE,sep="\t")
  
  sample_list <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/conditions/crpa_v_crpa.txt", header = TRUE, stringsAsFactors = FALSE)
  samples_a <- sample_list[which(sample_list$species == "a"),]
  samples_b <- sample_list[which(sample_list$species == "b"),]
  add_sequencing_round_to_model <- length(unique(samples_a$sequencing_round)) + length(unique(samples_b$sequencing_round))
  keeps <- c("Geneid", sample_list$sample)
  keeper <- sample_list$sample
  
  cts <- all_cts[keeps]
  #head(cts)
  write.table(cts, "C:/Users/jmcgirr/Documents/all_2018_samples/conditions/crpa_v_crpa_counts.txt", row.names = FALSE, quote= FALSE,sep="\t") 
  
  ### DIFFERENTIAL EXPRESSION ANALYSIS ###
  #vignette("DESeq2")
  cts <- as.matrix(read.table("C:/Users/jmcgirr/Documents/all_2018_samples/conditions/crpa_v_crpa_counts.txt" ,sep = "\t",header = TRUE,row.names=1))
  head(cts)
  nrow(cts)
  colData <- as.matrix(read.table("C:/Users/jmcgirr/Documents/all_2018_samples/conditions/crpa_v_crpa.txt" ,header = TRUE,row.names=1))
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
  prop_des <- c(prop_des, prop_de)
  
}  

mean(prop_des)
total_genes
range(prop_des)

####### prop DE/ME vs #########
library(evolqg)
library(ape)
library(geiger)
library(nlme)
library(phytools)
#install.packages("phangorn")
library(phangorn)
library(adephylo)
#install.packages("adephylo")

comps <-read.table("C:/Users/jmcgirr/Documents/all_2018_samples/prop_ME_CRPA_crosses.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
head(comps)
all_48 <- comps[which(comps$stage == "48hpf"),]
all_8 <- comps[which(comps$stage == "8dpf"),]
generalists_48 <- comps[which(comps$stage == "48hpf" & comps$specialist == "n"),]
generalists_8 <- comps[which(comps$stage == "8dpf" & comps$specialist == "n"),]

plot(comps$distance_km, comps$prop_DE)
plot(comps$distance_branchlength, comps$prop_DE)

#plot(generalists_48$distance_km, generalists_48$prop_DE)
plot(all_48$distance_km, all_48$prop_DE)
#plot(generalists_8$distance_km, generalists_8$prop_DE)
plot(all_8$distance_km, all_8$prop_DE)

red <- "#E60012" 
blu <- "#00A0E9" 
yel <- "#f6c700"
bla <- "#000000" 
#pur <- "#920783"
pur <- "#9933FF"
ora <- "#FF6600" 
gre <- "#22AC38" 
dkr <- "#800000" 
whi <- "#ffffff"
lir <- "#FFCCCC"

red <- "#FF0000" 
blu <- "#0000FF" 
yel <- "#FFFF00"
bla <- "#000000" 
pur <- "#FF00FF"
grb <- "#00FFFF" 
gre <- "#00FF00" 
dkr <- "#8B4513" 
whi <- "#ffffff"
lir <- "#ffb6c1"

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

all_cols <- c(dkr,dkr,lir,lir,yel,yel,pur,pur,bla,bla,blu,blu,red,red,whi,whi,gre,gre,grb,grb)
sep_cols <- c(dkr,lir,yel,pur,bla,blu,red,whi,gre,grb)
#plot(generalists_48$distance_branchlength, generalists_48$prop_DE)

par(mfrow=c(2,2))

comps <-read.table("C:/Users/jmcgirr/Documents/all_2018_samples/prop_DE_CRPA_crosses.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
all_48 <- comps[which(comps$stage == "48hpf"),]
all_8 <- comps[which(comps$stage == "8dpf"),]
gen_48 <- all_48[which(all_48$specialist == "n"),]
gen_8 <- all_8[which(all_8$specialist == "n"),]

#tiff("C:/Users/jmcgirr/Documents/all_2018_samples/manuscript_figs/de_me_48_8.tiff", width = 5, height = 5, units = 'in', res = 1000)
sep_cols <- c(dkr,lir,yel,pur,bla,blu,red,whi,gre,grb)
par(mfrow=c(2,2))
par(mai=c(0.2,0.4,0.2,0))
lm_all <- lm(all_48$distance_branchlength~all_48$prop_DE)
lm_gen <- lm(gen_48$distance_branchlength~gen_48$prop_DE)

plot(jitter(all_48$distance_branchlength), jitter((all_48$prop_DE)*100),pch = c(21,24)[as.factor(all_48$point_shape)],
     bg =sep_cols,col = "black", cex = 2.3, ylim = c(0,31),xlim = c(-0.02,0.45),xlab = "", ylab = "",xaxt ="n",cex.axis=1.5, yaxt = "n")
axis(side=2, at=c(0,10,20,30), labels = TRUE, cex.axis= 1.5)
abline(lm_all)
abline(lm_gen)

lm_all <- lm(all_8$distance_branchlength~all_8$prop_DE)
lm_gen <- lm(gen_8$distance_branchlength~gen_8$prop_DE)
par(mai=c(0.2,0.2,0.2,0.2))
plot(jitter(all_8$distance_branchlength), jitter((all_8$prop_DE)*100),pch = c(21,24)[as.factor(all_8$point_shape)],
     bg =sep_cols,col = "black", cex = 2.3, ylim = c(0,31),xlim = c(-0.02,0.45),xlab = "", ylab = "", yaxt = "n",xaxt ="n")
#abline(lm_all)
#abline(lm_gen)


comps <-read.table("C:/Users/jmcgirr/Documents/all_2018_samples/prop_ME_CRPA_crosses.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
head(comps)
all_48 <- comps[which(comps$stage == "48hpf" & comps$comp_type == "ph"),]
all_8 <- comps[which(comps$stage == "8dpf"& comps$comp_type == "ph"),]
sep_cols <- c(dkr,lir,yel,pur,grb)

par(mai=c(0.4,0.4,0,0))
plot(jitter(all_48$distance_branchlength), jitter((all_48$prop_DE)*100),pch = c(21,24)[as.factor(all_48$point_shape)],
     bg =sep_cols,col = "black", cex = 2.3, ylim = c(0,31),xlim = c(-0.02,0.45),xlab = "", ylab = "",cex.axis=1.5, yaxt = "n")
axis(side=2, at=c(0,10,20,30), labels = TRUE, cex.axis= 1.5)


par(mai=c(0.4,0.2,0,0.2))
plot(jitter(all_8$distance_branchlength), jitter((all_8$prop_DE)*100),pch = c(21,24)[as.factor(all_8$point_shape)],
     bg =sep_cols,col = "black", cex = 2.3, ylim = c(0,31),xlim = c(-0.02,0.45),xlab = "", ylab = "", yaxt = "n",cex.axis=1.5)
#dev.off()



## PGLS phylogenetic least squares ##https://www.r-phylo.org/wiki/HowTo/PGLS
## DE
all_48 <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/mse_v_distance/de_v_crpa_48hpf.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
all_8 <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/mse_v_distance/de_v_crpa_8dpf.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
row.names(all_48) <- all_48$tree_tip
row.names(all_8) <- all_8$tree_tip
all_48_de <- all_48
all_8_de <- all_8

# ME
all_48 <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/mse_v_distance/mse_v_crpa_48hpf.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
all_8 <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/mse_v_distance/mse_v_crpa_8dpf.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
tree <- read.tree("D:/Martin Lab/rna_2018/raxml/second_try_dna/RAxML_bestTree.58_dna_no_missing")
distances <- c(cophenetic(tree)["CAF1","CPF1"],cophenetic(tree)["CAF1","CMF1"],0,
               cophenetic(tree)["CAF1","CUNP"],cophenetic(tree)["CAF1","NAF1"])
all_48$distance_branch <- distances
all_8$distance_branch <- distances
row.names(all_48) <- all_48$tree_tip
row.names(all_8) <- all_8$tree_tip
all_48_me <- all_48
all_8_me <- all_8

all_48_de$prop_DE_vs_CRPA_scaled <- all_48_de$prop_DE_vs_CRPA *100
all_8_de$prop_DE_vs_CRPA_scaled <- all_8_de$prop_DE_vs_CRPA *100
all_48_me$prop_DE_vs_CRPA_scaled <- all_48_me$prop_DE_vs_CRPA *100
all_8_me$prop_DE_vs_CRPA_scaled <- all_8_me$prop_DE_vs_CRPA *100

## plot km on x axis
{

par(mfrow=c(2,2))

comps <-read.table("C:/Users/jmcgirr/Documents/all_2018_samples/prop_DE_CRPA_crosses.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
all_48 <- comps[which(comps$stage == "48hpf"),]
all_8 <- comps[which(comps$stage == "8dpf"),]
points_pch <- c(21,21,21,21,1,1,1,1,1,21)
sep_cols <- c("black","black",yel,blu,"black","black","black","black","black", red)

#tiff("C:/Users/jmcgirr/Documents/all_2018_samples/manuscript_figs/de_me_48_8_shapes.tiff", width = 5, height = 5, units = 'in', res = 1000)
par(mfrow=c(2,2))
par(mai=c(0.2,0.4,0.2,0))

plot(jitter(all_48$distance_km,50), all_48$prop_DE*100,col = "black",pch = points_pch,bg = sep_cols,
     cex = 2.3, ylim = c(0,31),xlab = "", ylab = "",xaxt ="n",cex.axis=1.5, yaxt = "n", xlim = c(-50,1200))
axis(side=2, at=c(0,10,20,30), labels = TRUE, cex.axis= 1.3)
pglsModel <- gls(prop_DE_vs_CRPA_scaled ~ distance_km, correlation = corBrownian(phy = all_tree),
                 data = all_48_de, method = "ML")
#summary(pglsModel)
abline(a = coef(pglsModel)[1], b = coef(pglsModel)[2], lty = 3)

par(mai=c(0.2,0.2,0.2,0.2))
plot(jitter(all_8$distance_km,50), all_8$prop_DE*100,col = "black",pch = points_pch,bg = sep_cols,
     cex = 2.3, ylim = c(0,31),xlab = "", ylab = "",xaxt ="n",cex.axis=1.5, yaxt = "n", xlim = c(-50,1200))
pglsModel <- gls(prop_DE_vs_CRPA_scaled ~ distance_km, correlation = corBrownian(phy = all_tree),
                 data = all_8_de, method = "ML")
#summary(pglsModel)
#abline(a = coef(pglsModel)[1], b = coef(pglsModel)[2], lty = 3)

comps <-read.table("C:/Users/jmcgirr/Documents/all_2018_samples/prop_ME_CRPA_crosses.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
head(comps)
all_48 <- comps[which(comps$stage == "48hpf"),]
all_8 <- comps[which(comps$stage == "8dpf"),]
sep_cols <- c("grey","grey",gre,grb,"grey","grey","grey","grey","grey", lir)

par(mai=c(0.4,0.4,0,0))
plot(jitter(all_48$distance_km,50), all_48$prop_DE*100,col = "black",pch = points_pch,bg = sep_cols,
     cex = 2.3, ylim = c(0,31),xlab = "", ylab = "",cex.axis=1.5, yaxt = "n", xlim = c(-50,1200), xaxt = "n")
axis(side=2, at=c(0,10,20,30), labels = TRUE, cex.axis= 1.3)
axis(side=1, at=c(0,250,500,750,1000), labels = c(0,250,500,750,1000), cex.axis= 1)

pglsModel <- gls(prop_DE_vs_CRPA_scaled ~ distance_km, correlation = corBrownian(phy = all_tree),
                 data = all_48_me, method = "ML")
#summary(pglsModel)
#abline(a = coef(pglsModel)[1], b = coef(pglsModel)[2], lty = 3)

par(mai=c(0.4,0.2,0,0.2))
plot(jitter(all_8$distance_km,50), all_8$prop_DE*100,col = "black",pch = points_pch,bg = sep_cols,
     cex = 2.3, ylim = c(0,31),xlab = "", ylab = "",cex.axis=1.5, yaxt = "n", xlim = c(-50,1200), xaxt = "n")
axis(side=1, at=c(0,250,500,750,1000), labels = c(0,250,500,750,1000), cex.axis= 1)
pglsModel <- gls(prop_DE_vs_CRPA_scaled ~ distance_km, correlation = corBrownian(phy = all_tree),
                 data = all_8_me, method = "ML")
#summary(pglsModel)
#abline(a = coef(pglsModel)[1], b = coef(pglsModel)[2], lty = 3)
#dev.off()
}





generalists_48 <- comps[which(comps$stage == "48hpf" & comps$specialist == "n"),]
generalists_8 <- comps[which(comps$stage == "8dpf" & comps$specialist == "n"),]

plot(comps$distance_km, comps$prop_DE)
plot(comps$distance_branchlength, comps$prop_DE)

#plot(generalists_48$distance_km, generalists_48$prop_DE)
plot(all_48$distance_km, all_48$prop_DE)
#plot(generalists_8$distance_km, generalists_8$prop_DE)
plot(all_8$distance_km, all_8$prop_DE)


#plot(generalists_48$distance_branchlength, generalists_48$prop_DE)
#tiff("C:/Users/jmcgirr/Documents/all_2018_samples/manuscript_figs/de_48.tiff", width = 4.5, height = 5, units = 'in', res = 1000)
plot(jitter(all_48$distance_branchlength), jitter((all_48$prop_DE)*100),pch = c(21,24)[as.factor(all_48$point_shape)],
     bg =sep_cols,col = "black", cex = 2.3, ylim = c(0,21.5),xlim = c(-0.02,0.45),xlab = "", ylab = "",xaxt ="n",cex.axis=1.5)
#dev.off()
#plot(generalists_8$distance_branchlength, generalists_8$prop_DE)
#tiff("C:/Users/jmcgirr/Documents/all_2018_samples/manuscript_figs/de_8.tiff", width = 4.5, height = 5, units = 'in', res = 1000)
plot(jitter(all_8$distance_branchlength), jitter((all_8$prop_DE)*100),pch = c(21,24)[as.factor(all_8$point_shape)],
     bg =sep_cols,col = "black", cex = 2.3, ylim = c(0,21.5),xlim = c(-0.02,0.45),xlab = "", ylab = "", yaxt = "n",xaxt ="n")
#dev.off()



## PGLS phylogenetic least squares ##https://www.r-phylo.org/wiki/HowTo/PGLS
# ME
all_48 <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/mse_v_distance/mse_v_crpa_48hpf.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
all_8 <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/mse_v_distance/mse_v_crpa_8dpf.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
tree <- read.tree("D:/Martin Lab/rna_2018/raxml/second_try_dna/RAxML_bestTree.58_dna_no_missing")
distances <- c(cophenetic(tree)["CAF1","CPF1"],cophenetic(tree)["CAF1","CMF1"],0,
               cophenetic(tree)["CAF1","CUNP"],cophenetic(tree)["CAF1","NAF1"])
all_48$distance_branch <- distances
all_8$distance_branch <- distances
species<-c("CAF1",	"CMF1",	"CPF1",	"CUNP",	"NAF1")
all_tree<-drop.tip(tree,tree$tip.label[-match(species, tree$tip.label)])
species<-c("CAF1",	"CUNP",	"NAF1")
generalist_tree<-drop.tip(tree,tree$tip.label[-match(species, tree$tip.label)])
generalists_48 <- all_48[which(all_48$specialist == "n"),]
generalists_8 <- all_8[which(all_8$specialist == "n"),]
test<-c("CAF1",	"CUNP",	"NAF1")
row.names(all_48) <- all_48$tree_tip
row.names(all_8) <- all_8$tree_tip
row.names(generalists_48) <- test
row.names(generalists_8) <- test
all_48 <- all_48_me
all_8 <- all_8_me

pglsModel <- gls(prop_DE_vs_CRPA ~ distance_km, correlation = corBrownian(phy = generalist_tree),
                 data = generalists_48, method = "ML")
summary(pglsModel)
pglsModel <- gls(prop_DE_vs_CRPA ~ distance_km, correlation = corBrownian(phy = generalist_tree),
                 data = generalists_8, method = "ML")
summary(pglsModel)
pglsModel <- gls(prop_DE_vs_CRPA ~ distance_km, correlation = corBrownian(phy = all_tree),
                 data = all_48, method = "ML")
summary(pglsModel)
pglsModel <- gls(prop_DE_vs_CRPA ~ distance_km, correlation = corBrownian(phy = all_tree),
                 data = all_8, method = "ML")
summary(pglsModel)

## DE
all_48 <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/mse_v_distance/de_v_crpa_48hpf.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
all_8 <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/mse_v_distance/de_v_crpa_8dpf.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
generalists_48 <- all_48[which(all_48$specialist == "n"),]
generalists_8 <- all_8[which(all_8$specialist == "n"),]
row.names(all_48) <- all_48$tree_tip
row.names(all_8) <- all_8$tree_tip
test<-c("NAF1",	"CUNP",	"CAF1")
row.names(generalists_48) <- test
row.names(generalists_8) <- test
all_48 <- all_48_me
all_8 <- all_8_me

pglsModel <- gls(prop_DE_vs_CRPA ~ distance_km, correlation = corBrownian(phy = generalist_tree),
                 data = generalists_48, method = "ML")
summary(pglsModel)
pglsModel <- gls(prop_DE_vs_CRPA ~ distance_km, correlation = corBrownian(phy = generalist_tree),
                 data = generalists_8, method = "ML")
summary(pglsModel)
pglsModel <- gls(prop_DE_vs_CRPA ~ distance_km, correlation = corBrownian(phy = all_tree),
                 data = all_48, method = "ML")
summary(pglsModel)
pglsModel <- gls(prop_DE_vs_CRPA ~ distance_km, correlation = corBrownian(phy = all_tree),
                 data = all_8, method = "ML")
summary(pglsModel)

#plot tree 
species<-c("CAF1",	"CMF1",	"CPF1",	"CUNP",	"NAF1", "ETAA1","OAF1","OMF1","OPF1")
map_tree<-drop.tip(tree,tree$tip.label[-match(species, tree$tip.label)])
map_tree<-rotate(map_tree,15)
nodelabels()
#tiff("D:/Martin Lab/rna_2018/raxml/map_tree_cp_and_op.tiff", width = 4.5, height = 4.5, units = 'in', res = 1000)
#png("D:/Martin Lab/rna_2018/raxml/map_tree_cp_and_op.png", width = 5, height = 4.5, units = 'in', res = 1000,bg = "transparent")
plotTree(map_tree,offset=1,direction = "leftwards")
nodelabels(c("",100,100,100,100,100,100,100),frame="none",adj=c(1.2,0),cex = 0.4 )
add.scale.bar(0.6,7)
#dev.off()

#png("D:/Martin Lab/rna_2018/raxml/map_tree_cp_and_op.png", width = 6, height = 10, units = 'in', res = 1000,bg = "transparent")
plotTree(tree,offset=1)
#dev.off()



# phylogenetic mantel test 
tree <- read.tree("D:/Martin Lab/rna_2018/raxml/second_try_dna/RAxML_bestTree.58_dna_no_missing")
distance_km <-as.matrix(read.table("C:/Users/jmcgirr/Documents/all_2018_samples/mse_v_distance/distance_km_matrix.txt" ,header = TRUE,row.names=1))
mis_48 <-as.matrix(read.table("C:/Users/jmcgirr/Documents/all_2018_samples/mse_v_distance/misexpression_48hpf_matrix.txt" ,header = TRUE,row.names=1))
mis_8 <-as.matrix(read.table("C:/Users/jmcgirr/Documents/all_2018_samples/mse_v_distance/misexpression_8dpf_matrix.txt" ,header = TRUE,row.names=1))
#mis_48 <-as.matrix(read.table("C:/Users/jmcgirr/Documents/all_2018_samples/mse_v_distance/na_rm_test_matrix.txt" ,header = TRUE,row.names=1))

MyTree <- read.tree("D:/Martin Lab/rna_2018/raxml/second_try_dna/RAxML_bestTree.58_dna_no_missing")
head(MyTree)

tree = MyTree
species<-c("CAF1",	"CMF1",	"CPF1",	"CUNP",	"NAF1")
tree<-drop.tip(tree,tree$tip.label[-match(species, tree$tip.label)])

PhyloMantel(tree, mis_48, distance_km, k = 10000)
MantelCor(mis_48, distance_km)
## End(Not run)


#####
