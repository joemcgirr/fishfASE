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
gem <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/GEMMA/output_58/mean_gammas_up_jaw_summed_20kb.txt", header = FALSE, stringsAsFactors = FALSE)
head(gem)
gem_main <- cbind(gem, colsplit(gem[,3], "\\|", c("end", "pip")))

cool_genes_out_dir <- "C:/Users/jmcgirr/Documents/all_2018_samples/manuscript_figs/fig_5_can_genes/"
cool_gene_name <- "cyp8b1"


cool_gene <- "XM_015369861.1"
pop_gen_comp <- "omxop"
cool_stage <- "8dpf"
p1_col <- gre
p2_col <- blu

#load pop comp
{
  
  p1 <- strsplit(pop_gen_comp, "x")[[1]][1]
  p2 <- strsplit(pop_gen_comp, "x")[[1]][2]
  #p1_species <- strsplit(species_comp, "")[[1]][1]
  #p2_species <- strsplit(species_comp, "")[[1]][2]
  dxy_main <- read.csv(paste("D:/Martin Lab/rna_2018/fst/dna/",pop_gen_comp,"_popgen_dna_stats_corr_dxy.csv",sep = ""), header = TRUE, stringsAsFactors = FALSE)
  #dxy_species_main <- read.csv(paste("D:/Martin Lab/rna_2018/fst/dna/",species_comp,"_popgen_dna_stats_corr_dxy.csv",sep = ""), header = TRUE, stringsAsFactors = FALSE)
  
  fst_main <- read.table(paste("C:/Users/jmcgirr/Documents/all_2018_samples/fst_dna/",pop_gen_comp,"_fst",sep = ""), header = TRUE, stringsAsFactors = FALSE)
  fst_main <- fst_main[which(fst_main$WEIR_AND_COCKERHAM_FST >=0 & fst_main$WEIR_AND_COCKERHAM_FST <=1),]
  #fst_species_main <- read.table(paste("C:/Users/jmcgirr/Documents/all_2018_samples/fst_dna/",species_comp,"_fst",sep = ""), header = TRUE, stringsAsFactors = FALSE)
  #fst_species_main <- fst_species_main[which(fst_species_main$WEIR_AND_COCKERHAM_FST >=0 & fst_species_main$WEIR_AND_COCKERHAM_FST <=1),]
  
  lake <- strsplit(pop_gen_comp,split ="")[[1]][1]
  taj_m_main <-    read.table(paste("C:/Users/jmcgirr/Documents/all_2018_samples/fst_dna/",lake,"m_taj_d_20kb.txt.Tajima.D", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  taj_a_main <-    read.table(paste("C:/Users/jmcgirr/Documents/all_2018_samples/fst_dna/",lake,"a_taj_d_20kb.txt.Tajima.D", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  taj_p_main <-    read.table(paste("C:/Users/jmcgirr/Documents/all_2018_samples/fst_dna/",lake,"p_taj_d_20kb.txt.Tajima.D", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  taj_m_main <- na.omit(taj_m_main)
  taj_a_main <- na.omit(taj_a_main)
  taj_p_main <- na.omit(taj_p_main)
  #taj_m_species_main <-    read.table("C:/Users/jmcgirr/Documents/all_2018_samples/fst_dna/m_taj_d_20kb.txt.Tajima.D", header = TRUE, stringsAsFactors = FALSE)
  #taj_a_species_main <-    read.table("C:/Users/jmcgirr/Documents/all_2018_samples/fst_dna/a_taj_d_20kb.txt.Tajima.D", header = TRUE, stringsAsFactors = FALSE)
  #taj_p_species_main <-    read.table("C:/Users/jmcgirr/Documents/all_2018_samples/fst_dna/p_taj_d_20kb.txt.Tajima.D", header = TRUE, stringsAsFactors = FALSE)
  #taj_m_species_main <- na.omit(taj_m_species_main)
  #taj_a_species_main <- na.omit(taj_a_species_main)
  #taj_p_species_main <- na.omit(taj_p_species_main)
  
  
  sweed_p1_main <-    read.table(paste("C:/Users/jmcgirr/Documents/all_2018_samples/sweed/correct_seq/SweeD_Report.",p1,"_pop_bottle_58_grid_500_sweeps.txt", sep = ""), header = FALSE, stringsAsFactors = FALSE)
  sweed_p2_main <-    read.table(paste("C:/Users/jmcgirr/Documents/all_2018_samples/sweed/correct_seq/SweeD_Report.",p2,"_pop_bottle_58_grid_500_sweeps.txt", sep = ""), header = FALSE, stringsAsFactors = FALSE)
  #sweed_p1_species_main <-    read.table(paste("C:/Users/jmcgirr/Documents/all_2018_samples/sweed/correct_seq/SweeD_Report.",p1_species,"_pop_bottle_58_grid_500_sweeps.txt", sep = ""), header = FALSE, stringsAsFactors = FALSE)
  #sweed_p2_species_main <-    read.table(paste("C:/Users/jmcgirr/Documents/all_2018_samples/sweed/correct_seq/SweeD_Report.",p2_species,"_pop_bottle_58_grid_500_sweeps.txt", sep = ""), header = FALSE, stringsAsFactors = FALSE)
  
}

#tiff(paste(cool_genes_out_dir,cool_gene_name,"_counts.tiff", sep = ""), width = 6, height = 5, units = 'in', res = 1000)
#plot counts
{
  
  cross_counts <- c("CRPA",
                    "OSPA",
                    "CRPM",
                    "OSPM",
                    "CRPP",
                    "OSPP",
                    "CAxCM",
                    "CMxCA",
                    "OAxOM",
                    "OMxOA",
                    "CAxCP",
                    "OAxOP",
                    "OPxOA",
                    "CMxCP",
                    "OPxOM",
                    "UPxUA",
                    "NCA")
  cross_counts <- c("CRPA",
                    "CRPM",
                    "CRPP",
                    "CAxCM",
                    "CMxCA",
                    "CAxCP",
                    "CMxCP")
  
  cross_counts <- c("OSPA",
                    "OSPM",
                    "OSPP",
                    "OAxOM",
                    "OMxOA",
                    "OAxOP",
                    "OPxOA",
                    "OPxOM")
  
  
  #cool_cols <- c(red,red,gre,gre,blu,blu,yel,yel,yel,yel,pur,pur,pur,grb,grb,"white","black")
  cool_cols <- c(red,gre,blu,yel,yel,pur,pur,grb)
  
  type_inds <- c()
  inds <- c()
  p1_cts <- c()
  for (i in cross_counts)
  {
    
    p1_inds <- master[which(master$f1 == i & master$stage == cool_stage),]
    p1_inds <- p1_inds$sample
    norm_cts_gene <- norm_cts[which(norm_cts$related_accession == cool_gene),]
    p1_cts_inds <- norm_cts_gene[p1_inds]
    p1_cts  <- c(p1_cts,as.numeric(p1_cts_inds[1,]))
    inds <- c(inds,p1_inds)
    type_inds <- c(type_inds,rep(i,length(p1_inds)))
    
  }
  
  plot_counts <- data.frame(type_ind = type_inds, sample = inds, cts = p1_cts)
  plot_counts$type_ind <- as.character(plot_counts$type_ind)
  plot_counts$type_ind <- factor(plot_counts$type_ind, levels=unique(plot_counts$type_ind))
  #plot_counts$type_ind <- factor(plot_counts$type_ind, level = cross_counts)
  #plot_title <- paste(cool_gene,norm_cts_gene$gene_name[1],"\n", why_cool, sep = "")
  #cool_gene_out <- paste(cool_genes_out_dir,norm_cts_gene$gene_name[1],"_",cool_h,"_",cool_stage, sep = "")
  
  #tiff(paste(cool_gene_out, "_gene_counts.tiff", sep = ""), width = 6.2, height = 5, units = 'in', res = 1000)
  plot(plot_counts$type_ind,plot_counts$cts, ylab = "normalized counts", border = "white",
       cex.axis = 0.5, cex.names = 1.1,cex.lab = 1.2, yaxt = "n")
  axis(2, cex.axis = 1.2, at = NULL)
  
  for (j in c(1,2,3))
  {
    p1_cts <- plot_counts[which(plot_counts$type_ind == cross_counts[j]),]
    points(jitter(rep(j,nrow(p1_cts)),1), p1_cts$cts, pch = 21,bg =cool_cols[j],  col = "black", cex = 2)
  }
  for (j in c(4,5,6,7,8))
  {
    p1_cts <- plot_counts[which(plot_counts$type_ind == cross_counts[j]),]
    points(jitter(rep(j,nrow(p1_cts)),1), p1_cts$cts, pch = 22,bg =cool_cols[j],  col = "black", cex = 2)
  }
  
}
#dev.off()

#tiff(paste(cool_genes_out_dir,cool_gene_name,"_popgen.tiff", sep = ""), width = 3, height = 5, units = 'in', res = 1000)
xrange <- c(1670000,1970000)
par(mfrow = c(5,1),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
#plot pop gen
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
       xlab = "", ylab = "fst", yaxt = 'n',xlim =xrange, xaxt = "n")
  axis(2,at=c(0,0.5, 1), labels=c(0,0.5,1),cex.axis=1.2)
  lines(smooth_fst, lwd = 2, col = "black")
  rect(gene_chrom$Start[1], -10, gene_chrom$End[1], 1.01, border = NA, col = col2alpha(blu,0.2))
  
  #fst <- fst_species_main[which(fst_species_main$CHROM == gene_chrom$Chr[1]),]
  #fixed <- fst[which(fst$WEIR_AND_COCKERHAM_FST == 1),]
  #smooth_fst <-   smooth.spline(fst$POS, fst$WEIR_AND_COCKERHAM_FST, spar = .2)
  #plot(fst$POS, fst$WEIR_AND_COCKERHAM_FST,ylim = c(0,1),
  #     xlab = "", ylab = "fst", yaxt = 'n',xlim =xrange, xaxt = "n",col = blu)
  #axis(2,c(0,0.5,1),cex.axis=1.2)
  #lines(smooth_fst, lwd = 2, col = red)
  #rect(gene_chrom$Start[1], -1000, gene_chrom$End[1], 1000, border = NA, col = col2alpha(blu,0.2))
  
  plot(dxy$start, dxy$corr_dxy,col = "white",
       xlab = "", ylab = "dxy",xlim =xrange, xaxt = "n",cex.axis=1.2,ylim = c(0,0.02), yaxt="n")
  axis(2, at=c(0,0.01, 0.02), labels=c(0,0.01,0.02),cex.axis=1.2)
  lines(smooth_dxy, lwd = 2, col = "grey")
  #dxy <- dxy_species_main[which(dxy_species_main$scaffold == gene_chrom$Chr[1]),]
  #smooth_dxy <-   smooth.spline(dxy$start, dxy$corr_dxy, spar = .2)
  #lines(smooth_dxy, lwd = 2, col = "grey",lty = 2)
  rect(gene_chrom$Start[1], -10, gene_chrom$End[1], 10, border = NA, col = col2alpha(blu,0.2))
  
  dxy <- dxy_main[which(dxy_main$scaffold == gene_chrom$Chr[1]),]
  #plot(dxy$start, dxy[,14],col = "white",
  #     xlab = "", ylab = "pi",xlim =xrange, xaxt = "n",cex.axis=1.2)#, yaxt="n")
  #axis(2, at=c(0,0.005,0.01, 0.015), labels=c(0,"",0.01,""),cex.axis=1.2)
  lines(smooth_p1_pi, lwd = 2, col = p1_col)
  lines(smooth_p2_pi, lwd = 2, col = p2_col)
  #dxy <- dxy_species_main[which(dxy_species_main$scaffold == gene_chrom$Chr[1]),]
  #smooth_p1_pi <- smooth.spline(dxy$start, dxy[,14], spar = .2)
  #smooth_p2_pi <- smooth.spline(dxy$start, dxy[,16], spar = .2)
  #lines(smooth_p1_pi, lwd = 2, col = p1_col, lty = 2)
  #lines(smooth_p2_pi, lwd = 2, col = p2_col, lty = 2)
  
  plot(tj_p$BIN_START, tj_p$TajimaD,ylab = "D",ylim = c(-2,2.5),xlab = "",
       xlim = xrange,col = "white", yaxt = "n", xaxt = "n")
  axis(2,labels = c(-1,0,1,2),at = c(-1,0,1,2),cex.axis=1.2)
  lines(smooth_m, lwd = 2, col = gre)
  lines(smooth_p, lwd = 2, col = blu)
  #lines(smooth_a, lwd = 2, col = red)
  #tj_m <- taj_m_species_main[which(taj_m_species_main$CHROM == gene_chrom$Chr[1]),]
  #smooth_m <-   smooth.spline(tj_m$BIN_START, tj_m$TajimaD, spar = .1)
  #tj_p <- taj_p_species_main[which(taj_p_species_main$CHROM == gene_chrom$Chr[1]),]
  #smooth_p <-   smooth.spline(tj_p$BIN_START, tj_p$TajimaD, spar = .1)
  #tj_a <- taj_a_species_main[which(taj_a_species_main$CHROM == gene_chrom$Chr[1]),]
  #smooth_a <-   smooth.spline(tj_a$BIN_START, tj_a$TajimaD, spar = .1)
  #lines(smooth_m, lwd = 2, col = gre, lty = 2)
  #lines(smooth_p, lwd = 2, col = blu, lty = 2)
  #lines(smooth_a, lwd = 2, col = red, lty = 2)
  rect(gene_chrom$Start[1], -10, gene_chrom$End[1], 10, border = NA, col = col2alpha(blu,0.2))
  
  sweed_p1$mid <- sweed_p1$V2 + ((sweed_p1$V2[2]-sweed_p1$V2[1])/2)
  sweed_p2$mid <- sweed_p2$V2 + ((sweed_p2$V2[2]-sweed_p2$V2[1])/2)
  sweed_p1$norm_clr = (sweed_p1$V3-min(sweed_p1$V3))/(max(sweed_p1$V3)-min(sweed_p1$V3))
  sweed_p2$norm_clr = (sweed_p2$V3-min(sweed_p2$V3))/(max(sweed_p2$V3)-min(sweed_p2$V3))
  smooth_sweed_p1 <-   smooth.spline(sweed_p1$mid, sweed_p1$norm_clr, spar = 0.1)
  smooth_sweed_p2 <-   smooth.spline(sweed_p2$mid, sweed_p2$norm_clr, spar = 0.1)
  
  plot(sweed_p1$mid, sweed_p1$norm_clr,xlab = "",ylab = "CLR",xlim = xrange, pch = 16,yaxt = "n",
       col = "white", cex.axis = 1.2, xaxt = "n")
  axis(2, at=c(0,0.5,1), labels=c(0,0.5,1),cex.axis=1.2)
  lines(smooth_sweed_p1, lwd = 2, col = p1_col)
  lines(smooth_sweed_p2, lwd = 2, col = p2_col)
  #sweed_p1 <- sweed_p1_species_main[which(sweed_p1_species_main$V1 ==gene_chrom$Chr[1]),]
  #sweed_p2 <- sweed_p2_species_main[which(sweed_p2_species_main$V1 ==gene_chrom$Chr[1]),]
  #sweed_p1$mid <- sweed_p1$V2 + ((sweed_p1$V2[2]-sweed_p1$V2[1])/2)
  #sweed_p2$mid <- sweed_p2$V2 + ((sweed_p2$V2[2]-sweed_p2$V2[1])/2)
  #sweed_p1$norm_clr = (sweed_p1$V3-min(sweed_p1$V3))/(max(sweed_p1$V3)-min(sweed_p1$V3))
  #sweed_p2$norm_clr = (sweed_p2$V3-min(sweed_p2$V3))/(max(sweed_p2$V3)-min(sweed_p2$V3))
  #smooth_sweed_p1 <-   smooth.spline(sweed_p1$mid, sweed_p1$norm_clr, spar = 0.1)
  #smooth_sweed_p2 <-   smooth.spline(sweed_p2$mid, sweed_p2$norm_clr, spar = 0.1)
  #lines(smooth_sweed_p1, lwd = 2, col = p1_col, lty = 2)
  #lines(smooth_sweed_p2, lwd = 2, col = p2_col, lty = 2)
  rect(gene_chrom$Start[1], -10, gene_chrom$End[1], 10, border = NA, col = col2alpha(blu,0.2))
  
  gem <- gem_main[which(gem_main$V1 ==gene_chrom$Chr[1]),]
  smooth_gem <-   smooth.spline(gem$V2, gem$pip, spar = .1)
  options(scipen=999)
  plot(gem$V2, gem$pip,xlab = "",ylab = "PIP",xlim = xrange, 
       col = "white", xaxt = "n",cex.axis=1.2, yaxt = "n", ylim = c(0,0.002))
  lines(smooth_gem, lwd = 2, col = pur)
  axis(2, at=c(0,0.001,0.002), labels=c(0,0.001,0.002),cex.axis=1.2)
  axis(1,cex.axis=1.2)
  rect(gene_chrom$Start[1], -1, gene_chrom$End[1], 1, border = NA, col = col2alpha(blu,0.2))
  #plot(gem$V2, gem$pip,xlab = "",ylab = "PIP", 
  #      col = "white")#, yaxt = "n")
  #lines(smooth_gem, lwd = 2, col = pur)
  
}
#dev.off()




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


