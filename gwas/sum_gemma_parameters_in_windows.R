
# This script is used after trait_name_gemma_mean_parameters.py
# This script is used to sum GEMMA parameters across 20kb windows
#------------------------------------------------------------------------------------------
# set working directory / trait names / .gff file directory / scaffold sizes file directory
#------------------------------------------------------------------------------------------

wk_dir <- "C:/Users/jmcgirr/Documents/remote_pups/GEMMA/"
traits <- c("jaw","pigment", "nose_height", "nose_length")
sizes <- read.table("D:/Martin Lab/lots_of_pups_project/GEMMA/asm.racon.scaffsizes.txt", header=FALSE, stringsAsFactors = F)


#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("plyranges")

library("plyranges")

colnames(sizes) <- c("seqnames", "end")
sizes$start <- 0
bins <- sizes %>% as_granges() %>% tile_ranges(width = 20000L)

#trait <- "jaw"
for (trait in traits)
{
  snp_params <- read.table(paste(wk_dir,"mean_params_",trait,"_.bed", sep = ""), header = FALSE, row.names = NULL, stringsAsFactors = FALSE)
  colnames(snp_params) <- c("seqnames", "start","end","id","alpha","beta","gamma")
  snp_params$id <- NULL
  #write.table(snp_params, paste(wk_dir,trait,"_mean_10_runs.txt", sep = ""), quote = FALSE, row.names = FALSE, sep = '\t')
  print(trait)
  print('top hits snps')
  print(head(snp_params[order(snp_params$gamma,decreasing = TRUE),]))
  print('  ')
  snps <- snp_params %>% as_granges()  %>% filter_by_overlaps(bins)
  sums <- bins %>%
    join_overlap_inner(snps) %>%
    disjoin_ranges(n = n(), sum_alpha_20kb = sum(alpha), 
                   sum_beta_20kb = sum(beta), 
                   sum_gamma_20kb = sum(gamma))
  
  sums <-  as(sums, "data.frame")
  #write.table(sums, paste(wk_dir,trait,"_mean_10_runs_sum_20kb_windows.txt", sep = ""), quote = FALSE, row.names = FALSE, sep = '\t')
  print('top hits windows')
  print(head(sums[order(sums$sum_gamma_20kb,decreasing = TRUE),]))
}

### overlap with genes ###
colnames(genes) <- c("seqnames", "start","end","gene")
head(genes)
genes <- genes %>% as_granges()

for (trait in traits)
{
  snps <- read.table(paste(wk_dir,trait,"_mean_10_runs.txt", sep = ""), header = TRUE, stringsAsFactors = FALSE)  
  sums <- read.table(paste(wk_dir,trait,"_mean_10_runs_sum_20kb_windows.txt", sep = ""), header = TRUE, stringsAsFactors = FALSE)  
  snps <- snps %>% as_granges()
  sums <- sums %>% as_granges()
  
  snps <- join_overlap_inner(snps, genes)
  snps <- as(snps, "data.frame")
  snps <- snps[order(snps$gamma, decreasing = TRUE),]
  #write.table(snps, paste(wk_dir,trait,"_mean_10_runs_genes.txt", sep = ""), quote = FALSE, row.names = FALSE, sep = '\t')
  sums <- join_overlap_inner(sums, genes)
  sums <- as(sums, "data.frame")
  sums <- sums[order(sums$sum_gamma_20kb, decreasing = TRUE),]
  #write.table(sums, paste(wk_dir,trait,"_mean_10_runs_sum_20kb_windows_genes.txt", sep = ""), quote = FALSE, row.names = FALSE, sep = '\t')
  
  snp_thresh <- quantile(snps$gamma,.999)[[1]]
  snps <- snps[which(snps$gamma >= snp_thresh),]
  #write.table(snps, paste(wk_dir,trait,"_mean_10_runs_genes_99th_sig.txt", sep = ""), quote = FALSE, row.names = FALSE, sep = '\t')
  sum_thresh <- quantile(sums$sum_gamma_20kb,.999)[[1]]
  sums <- sums[which(sums$sum_gamma_20kb >= sum_thresh),]
  #write.table(sums, paste(wk_dir,trait,"_mean_10_runs_sum_20kb_windows_genes_99.9th_sig.txt", sep = ""), quote = FALSE, row.names = FALSE, sep = '\t')
  
}