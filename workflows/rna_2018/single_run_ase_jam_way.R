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

st_mis_ase_trans_genes<- c()
mis_ase_trans_genes<- c()
sl_mis_ase_trans_genes<- c()

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
st_mis_ase_trans_genes<- c(st_mis_ase_trans_genes,paste(unique(intersect(mis_and_ase$related_accession,trans$related_accession)), sep = ";",collapse=";"))


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
mis_ase_trans_genes     <- c(mis_ase_trans_genes,paste(unique(intersect(mis_and_ase$related_accession,trans$related_accession)), sep = ";",collapse=";"))



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
sl_mis_ase_trans_genes<- c(sl_mis_ase_trans_genes,paste(unique(intersect(mis_and_ase$related_accession,trans$related_accession)), sep = ";",collapse=";"))


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
                                  stringsAsFactors = FALSE)

#write.table(strict_table,        "C:/Users/jmcgirr/Documents/all_2018_samples/ase_data/ase_and_mse/strict_ase_unphased_mbased_transtest_by_individual_no_ase_all_parents.txt", row.names = FALSE, quote = FALSE, sep ="\t")
#write.table(lenient_table,       "C:/Users/jmcgirr/Documents/all_2018_samples/ase_data/ase_and_mse/THE_ONE_lenient_ase_unphased_mbased_transtest_by_individual_no_ase_all_parents.txt", row.names = FALSE, quote = FALSE, sep ="\t")
#write.table(super_lenient_table, "C:/Users/jmcgirr/Documents/all_2018_samples/ase_data/ase_and_mse/super_lenient_ase_unphased_mbased_transtest_by_individual_no_ase_all_parents.txt", row.names = FALSE, quote = FALSE, sep ="\t")
