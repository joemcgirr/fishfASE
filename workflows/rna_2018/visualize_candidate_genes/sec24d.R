cool_gene <- "XM_015397247.1"
cool_p1 <- "CRPM"
cool_p2 <- "CRPP"
cool_h <- "CMxCP"
cool_stage <- "8dpf"
cool_cols <- mxp_cols
pop_gen_comp <- "cmxcp"


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