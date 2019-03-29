# first run adaptive_incompatibilities_venn.R

#axm
#ai_genes <- intersect(intersect(caxcm_48_DE,caxcm_48_ME),cmxca_48_ME)
#ai_genes <- intersect(intersect(caxcm_8_DE,caxcm_8_ME),cmxca_8_ME)
#ai_genes <- intersect(intersect(oaxom_48_DE,oaxom_48_ME),omxoa_48_ME)
ai_genes <- intersect(intersect(oaxom_8_DE,oaxom_8_ME),omxoa_8_ME)
cool_p1 <- "OSPA"
cool_p2 <- "OSPM"
cool_h1 <- "OAxOM"
cool_h2 <- "OMxOA"
cool_stage <- "8dpf"
cool_cols <- axm_cols

##### plot counts #####
for (cool_gene in ai_genes)
{
p1_inds <- master[which(master$f1 == cool_p1 & master$stage == cool_stage),]
p2_inds <- master[which(master$f1 == cool_p2 & master$stage == cool_stage),]
hy1_inds <- master[which(master$f1 == cool_h1 & master$stage == cool_stage),]
hy2_inds <- master[which(master$f1 == cool_h2 & master$stage == cool_stage),]
p1_inds <- p1_inds$sample
p2_inds <- p2_inds$sample
hy1_inds <- hy1_inds$sample
hy2_inds <- hy2_inds$sample
cool_genes_out_dir <- "C:/Users/jmcgirr/Documents/all_2018_samples/ase_data/plots/cool_genes/"
norm_cts_gene <- norm_cts[which(norm_cts$related_accession == cool_gene),]
ph_cts <- norm_cts_gene[c(p1_inds,p2_inds,hy1_inds,hy2_inds)]
ph_cts  <- as.numeric(ph_cts[1,])
inds <- c(p1_inds,p2_inds,hy1_inds,hy2_inds)
type_inds <- c(rep(cool_p1,length(p1_inds)),rep(cool_p2,length(p2_inds)),rep(cool_h1,length(hy1_inds)),rep(cool_h2,length(hy2_inds)))
plot_counts <- data.frame(type_ind = type_inds, sample = inds, cts = ph_cts)
plot_counts$type_ind <- as.character(plot_counts$type_ind)
plot_counts$type_ind <- factor(plot_counts$type_ind, levels=unique(plot_counts$type_ind))
plot_counts$type_ind <- factor(plot_counts$type_ind, level = c(cool_p1,cool_h1,cool_h2,cool_p2))
plot_title <- paste(cool_gene,"\n",norm_cts_gene$gene_name[1], sep = "")
cool_gene_out <- paste(cool_genes_out_dir,norm_cts_gene$gene_name[1],"_",cool_h,"_",cool_stage, sep = "")
#tiff(paste(cool_gene_out, "_gene_counts.tiff", sep = ""), width = 5, height = 6, units = 'in', res = 300)
plot(plot_counts$type_ind,plot_counts$cts, ylab = "normalized counts", main = plot_title, border = "white")
p1_cts <- plot_counts[which(plot_counts$type_ind == cool_p1),]
points(rep(1,nrow(p1_cts)), p1_cts$cts, pch = 21,bg =cool_cols[1],  col = "black", cex = 2)
p2_cts <- plot_counts[which(plot_counts$type_ind == cool_p2),]
points(rep(4,nrow(p2_cts)), p2_cts$cts, pch = 21,bg =cool_cols[2],  col = "black", cex = 2)
hy1_cts <- plot_counts[which(plot_counts$type_ind == cool_h1),]
points(rep(2,nrow(hy1_cts)), hy1_cts$cts, pch = 21,bg =cool_cols[3],  col = "black", cex = 2)
hy2_cts <- plot_counts[which(plot_counts$type_ind == cool_h2),]
points(rep(3,nrow(hy2_cts)), hy2_cts$cts, pch = 21,bg =cool_cols[3],  col = "black", cex = 2)
#dev.off()
}

#######################

#axp
#ai_genes <- intersect(caxcp_48_DE,caxcp_48_ME)
ai_genes <- intersect(caxcp_8_DE,caxcp_8_ME)
cool_p1 <- "CRPA"
cool_p2 <- "CRPP"
cool_h <- "CAxCP"
cool_stage <- "8dpf"
cool_cols <- axp_cols

##### plot counts #####
for (cool_gene in ai_genes)
{
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
}
#######################

ai_genes <- intersect(intersect(oaxop_48_DE,oaxop_48_ME),opxoa_48_ME)
cool_p1 <- "OSPA"
cool_p2 <- "OSPP"
cool_h1 <- "OAxOP"
cool_h2 <- "OPxOA"
cool_stage <- "48hpf"
cool_cols <- axp_cols

##### plot counts #####
for (cool_gene in ai_genes)
{
  p1_inds <- master[which(master$f1 == cool_p1 & master$stage == cool_stage),]
  p2_inds <- master[which(master$f1 == cool_p2 & master$stage == cool_stage),]
  hy1_inds <- master[which(master$f1 == cool_h1 & master$stage == cool_stage),]
  hy2_inds <- master[which(master$f1 == cool_h2 & master$stage == cool_stage),]
  p1_inds <- p1_inds$sample
  p2_inds <- p2_inds$sample
  hy1_inds <- hy1_inds$sample
  hy2_inds <- hy2_inds$sample
  cool_genes_out_dir <- "C:/Users/jmcgirr/Documents/all_2018_samples/ase_data/plots/cool_genes/"
  norm_cts_gene <- norm_cts[which(norm_cts$related_accession == cool_gene),]
  ph_cts <- norm_cts_gene[c(p1_inds,p2_inds,hy1_inds,hy2_inds)]
  ph_cts  <- as.numeric(ph_cts[1,])
  inds <- c(p1_inds,p2_inds,hy1_inds,hy2_inds)
  type_inds <- c(rep(cool_p1,length(p1_inds)),rep(cool_p2,length(p2_inds)),rep(cool_h1,length(hy1_inds)),rep(cool_h2,length(hy2_inds)))
  plot_counts <- data.frame(type_ind = type_inds, sample = inds, cts = ph_cts)
  plot_counts$type_ind <- as.character(plot_counts$type_ind)
  plot_counts$type_ind <- factor(plot_counts$type_ind, levels=unique(plot_counts$type_ind))
  plot_counts$type_ind <- factor(plot_counts$type_ind, level = c(cool_p1,cool_h1,cool_h2,cool_p2))
  plot_title <- paste(cool_gene,"\n",norm_cts_gene$gene_name[1], sep = "")
  cool_gene_out <- paste(cool_genes_out_dir,norm_cts_gene$gene_name[1],"_",cool_h,"_",cool_stage, sep = "")
  #tiff(paste(cool_gene_out, "_gene_counts.tiff", sep = ""), width = 5, height = 6, units = 'in', res = 300)
  plot(plot_counts$type_ind,plot_counts$cts, ylab = "normalized counts", main = plot_title, border = "white")
  p1_cts <- plot_counts[which(plot_counts$type_ind == cool_p1),]
  points(rep(1,nrow(p1_cts)), p1_cts$cts, pch = 21,bg =cool_cols[1],  col = "black", cex = 2)
  p2_cts <- plot_counts[which(plot_counts$type_ind == cool_p2),]
  points(rep(4,nrow(p2_cts)), p2_cts$cts, pch = 21,bg =cool_cols[2],  col = "black", cex = 2)
  hy1_cts <- plot_counts[which(plot_counts$type_ind == cool_h1),]
  points(rep(2,nrow(hy1_cts)), hy1_cts$cts, pch = 21,bg =cool_cols[3],  col = "black", cex = 2)
  hy2_cts <- plot_counts[which(plot_counts$type_ind == cool_h2),]
  points(rep(3,nrow(hy2_cts)), hy2_cts$cts, pch = 21,bg =cool_cols[3],  col = "black", cex = 2)
  #dev.off()
}
#######################
ai_genes <- intersect(intersect(oaxop_8_DE,oaxop_8_ME),opxoa_8_ME)
cool_p1 <- "OSPA"
cool_p2 <- "OSPP"
cool_h1 <- "OAxOP"
cool_h2 <- "OPxOA"
cool_stage <- "8dpf"
cool_cols <- axp_cols

##### plot counts #####
for (cool_gene in ai_genes)
{
  p1_inds <- master[which(master$f1 == cool_p1 & master$stage == cool_stage),]
  p2_inds <- master[which(master$f1 == cool_p2 & master$stage == cool_stage),]
  hy1_inds <- master[which(master$f1 == cool_h1 & master$stage == cool_stage),]
  hy2_inds <- master[which(master$f1 == cool_h2 & master$stage == cool_stage),]
  p1_inds <- p1_inds$sample
  p2_inds <- p2_inds$sample
  hy1_inds <- hy1_inds$sample
  hy2_inds <- hy2_inds$sample
  cool_genes_out_dir <- "C:/Users/jmcgirr/Documents/all_2018_samples/ase_data/plots/cool_genes/"
  norm_cts_gene <- norm_cts[which(norm_cts$related_accession == cool_gene),]
  ph_cts <- norm_cts_gene[c(p1_inds,p2_inds,hy1_inds,hy2_inds)]
  ph_cts  <- as.numeric(ph_cts[1,])
  inds <- c(p1_inds,p2_inds,hy1_inds,hy2_inds)
  type_inds <- c(rep(cool_p1,length(p1_inds)),rep(cool_p2,length(p2_inds)),rep(cool_h1,length(hy1_inds)),rep(cool_h2,length(hy2_inds)))
  plot_counts <- data.frame(type_ind = type_inds, sample = inds, cts = ph_cts)
  plot_counts$type_ind <- as.character(plot_counts$type_ind)
  plot_counts$type_ind <- factor(plot_counts$type_ind, levels=unique(plot_counts$type_ind))
  plot_counts$type_ind <- factor(plot_counts$type_ind, level = c(cool_p1,cool_h1,cool_h2,cool_p2))
  plot_title <- paste(cool_gene,"\n",norm_cts_gene$gene_name[1], sep = "")
  cool_gene_out <- paste(cool_genes_out_dir,norm_cts_gene$gene_name[1],"_",cool_h,"_",cool_stage, sep = "")
  #tiff(paste(cool_gene_out, "_gene_counts.tiff", sep = ""), width = 5, height = 6, units = 'in', res = 300)
  plot(plot_counts$type_ind,plot_counts$cts, ylab = "normalized counts", main = plot_title, border = "white")
  p1_cts <- plot_counts[which(plot_counts$type_ind == cool_p1),]
  points(rep(1,nrow(p1_cts)), p1_cts$cts, pch = 21,bg =cool_cols[1],  col = "black", cex = 2)
  p2_cts <- plot_counts[which(plot_counts$type_ind == cool_p2),]
  points(rep(4,nrow(p2_cts)), p2_cts$cts, pch = 21,bg =cool_cols[2],  col = "black", cex = 2)
  hy1_cts <- plot_counts[which(plot_counts$type_ind == cool_h1),]
  points(rep(2,nrow(hy1_cts)), hy1_cts$cts, pch = 21,bg =cool_cols[3],  col = "black", cex = 2)
  hy2_cts <- plot_counts[which(plot_counts$type_ind == cool_h2),]
  points(rep(3,nrow(hy2_cts)), hy2_cts$cts, pch = 21,bg =cool_cols[3],  col = "black", cex = 2)
  #dev.off()
}
#######################

ai_genes <- intersect(cmxcp_48_DE,cmxcp_48_ME)
cool_p1 <- "CRPM"
cool_p2 <- "CRPP"
cool_h <- "CMxCP"
cool_stage <- "48hpf"
cool_cols <- mxp_cols

##### plot counts #####
for (cool_gene in ai_genes)
{
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
}
#######################
ai_genes <- intersect(omxop_48_DE,opxom_48_ME)
ai_genes <- intersect(cmxcp_8_DE,cmxcp_8_ME)
ai_genes <- intersect(omxop_8_DE,opxom_8_ME)

