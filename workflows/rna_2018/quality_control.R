
###### RNA ######
master <- read.table("D:/Martin Lab/rna_2018/all_2018_samples/table_maker_master_outlier_rm.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
map <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/mapping_stats_qc/mapping_stats.txt", stringsAsFactors = FALSE, header = TRUE, sep = "\t")
feat_counts <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/mapping_stats_qc/feature_counts/feature_counts.txt", stringsAsFactors = FALSE, header = TRUE, sep = "\t")
head(feat_counts)

mis_indivs <- c("CVT1","CVT2","CVT3","CWT1","CWT2","CWT3","CXT1","CXT2","CXT3","LFT1","LFT2","LFT3","LIT1","LIT2","LIT3","NAT1","NAT2","NAT3","OVT1","OVT2","OVT3","OYT1","OYT2","OYT3","PAT1","PAT2","PAT3","OYE1","OYE2","OYE3","OYE4","OYE5","CAT1","CAT2","CAT3","CPT1","CPT2","NCT1","NCT2","NCT3","OAT1","OAT2","OMT1","OMT2","OMT3","OPT1","OPT2","OPT3","CAE1","CAE2","CAE3","CAE4","CAE5","CPE1","CPE2","CPE3","CPE4","CPE5","NCE1","NCE2","NCE3","NCE4","NCE5","OAE1","OAE2","OAE3","OAE4","OME1","OME2","OME3","OME4","OME5","OPE1","OPE2","OPE3","OPE4","OPE5","CUT1","CUT2","CUT3","LGT1","LGT2","LGT3","LKT1","LKT2","LKT3","OUT1","OUT2","OUT3","OXT1","OXT2","OXT3","OZT1","OZT2","OZT3","CPU1","CPU3","CPU5","CVE1","CVE2","CVE5","CWE2","CWE3","CWE4","CXE2","CXE3","CXE4","LFE2","LFE3","LFE4","LGE3","LGE4","LIE2","LIE3","LIE5","LKE1","LKE2","LKE3","NAE1","NAE2","NAE4","OUE1","OUE3","OUE4","OVE1","OVE4","OVE5","OXE2","OZE2","OZE4","OZE5","PAE1","PAE2","PAE5","CMT1","CMT2","CMT3","CPT3","CQT1","CQT2","OAT3","CME1","CME2","CME5","CQE1","CQE2","CQE3")
#no lake crosses
#mis_indivs <-c("CAE1","CAE2","CAE3","CAE4","CAE5","CME1","CME2","CME5","CPE1","CPE2","CPE3","CPE4","CPE5","CQE1","CQE2","CQE3","NCE1","NCE2","NCE3","NCE4","NCE5","OAE1","OAE2","OAE3","OAE4","OME1","OME2","OME3","OME4","OME5","OPE1","OPE2","OPE3","OPE4","OPE5","CPU1","CPU3","CPU5","CVE1","CVE2","CVE5","CWE2","CWE3","CWE4","CXE2","CXE3","CXE4","NAE1","NAE2","NAE4","OUE1","OUE3","OUE4","OVE1","OVE4","OVE5","OXE2","OYE1","OYE2","OYE3","OYE4","OYE5","OZE2","OZE4","OZE5","PAE1","PAE2","PAE5","CAT1","CAT2","CAT3","CMT1","CMT2","CMT3","CPT1","CPT2","CPT3","CQT1","CQT2","NCT1","NCT2","NCT3","OAT1","OAT2","OAT3","OMT1","OMT2","OMT3","OPT1","OPT2","OPT3","CUT1","CUT2","CUT3","CVT1","CVT2","CVT3","CWT1","CWT2","CWT3","CXT1","CXT2","CXT3","NAT1","NAT2","NAT3","OUT1","OUT2","OUT3","OVT1","OVT2","OVT3","OXT1","OXT2","OXT3","OYT1","OYT2","OYT3","OZT1","OZT2","OZT3","PAT1","PAT2","PAT3")

map <- map[map$alt_name %in% mis_indivs,]
# number of total reads for misexpression dataset (149 indivs, 147 if PCA outliers removed)
sum(map$total_reads)

# difference in proportion of reads mapped to features?
master$alt_name <- master$sample
m <- merge(feat_counts,master, by = c("alt_name"))
m <- m[m$alt_name %in% mis_indivs,]
nrow(m)
tail(m)

t <- aov(m$percent_assigned~m$species_region)
summary(t)
t <- aov(m$percent_assigned~m$species)
summary(t)

n <- m[which(m$species != "h"),]
t <- aov(n$percent_assigned~n$species_region)
summary(t)
t <- aov(n$percent_assigned~n$species)
summary(t)

o <- n[which(n$species != "h" & n$species != "m" & n$species != "p"),]
t <- aov(o$percent_assigned~o$species_region)
summary(t)

library(ggplot2)
library(ggpubr)
require(gridExtra)
#c('****' = 1e-04, '***' = 0.001, '**' = 0.01, '*' = 0.05,
#install.packages("ggpubr")

info <- m
cols <- c(red,"darkgrey",pur,gre,"black",blu)

box_plts <- info[order(info$species_region),]

p1 <- ggboxplot(box_plts, "species_region", "percent_assigned",color = cols, ylab = "% reads assigned to features\n", xlab="species")+
  stat_compare_means(method = "anova", label.y = max(info$percent_assigned)+(max(info$percent_assigned)*.05),label.x = 2.25)+      # Add global p-value
  stat_compare_means(label = "p.signif", hide.ns = TRUE,method = "t.test",ref.group = ".all.", label.y = (max(info$percent_assigned)+(max(info$percent_assigned)*.01)))+
  scale_x_discrete(labels= c("generalists","new prov", "hybrids","molluscivore", "NC", "scale-eater"))
#tiff("C:/Users/jmcgirr/Documents/all_2018_samples/manuscript_figs/supp/prop_reads_assigned_species.tiff", width = 6, height = 4, units = 'in', res = 1000)
p1
#dev.off()
cols <- c("black","#000000")
box_plts <- info[order(info$stage),] 
p1 <- ggboxplot(box_plts, "stage", "percent_assigned",color = cols, ylab = "% reads assigned to features\n", xlab="stage")+
  #stat_compare_means(method = "t.test", label.y = max(info$percent_assigned)+(max(info$percent_assigned)*.05),label.x = 1.25)+      # Add global p-value
  stat_compare_means(label = "p.signif", method = "t.test",ref.group = ".all.", label.y = (max(info$percent_assigned)+(max(info$percent_assigned)*.01)))+
  scale_x_discrete(labels = c("2dpf", "8dpf"))
#tiff("C:/Users/jmcgirr/Documents/all_2018_samples/manuscript_figs/supp/qc/prop_reads_assigned_stage.tiff", width = 4, height = 4, units = 'in', res = 1000)
p1
#dev.off()

# regressions
cts_table <- read.table("D:/Martin Lab/rna_2018/all_2018_samples/conditions/all_48hpf_and_8dpf_counts_my_gtf_geneid_final", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
norm_cts <- read.table("D:/Martin Lab/rna_2018/all_2018_samples/WGCNA/all_samples_outlier_rm_size_factor_normalized_counts.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
info <- read.table("D:/Martin Lab/rna_2018/all_2018_samples/table_maker_master_outlier_rm.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
map <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/mapping_stats_qc/mapping_qc_stats_region.txt", stringsAsFactors = FALSE, header = TRUE, sep = "\t")
map$sample <- map$alt_name
map <- map[c("sample","percent_assigned","Total","species","region","species_region","species_region_lake_cross_separate","species_region_lake_cross_with_species")]
#mis_indivs <-c("CAE1","CAE2","CAE3","CAE4","CAE5","CME1","CME2","CME5","CPE1","CPE2","CPE3","CPE4","CPE5","CQE1","CQE2","CQE3","NCE1","NCE2","NCE3","NCE4","NCE5","OAE1","OAE2","OAE3","OAE4","OME1","OME2","OME3","OME4","OME5","OPE1","OPE2","OPE3","OPE4","OPE5","CPU1","CPU3","CPU5","CVE1","CVE2","CVE5","CWE2","CWE3","CWE4","CXE2","CXE3","CXE4","NAE1","NAE2","NAE4","OUE1","OUE3","OUE4","OVE1","OVE4","OVE5","OXE2","OYE1","OYE2","OYE3","OYE4","OYE5","OZE2","OZE4","OZE5","PAE1","PAE2","PAE5","CAT1","CAT2","CAT3","CMT1","CMT2","CMT3","CPT1","CPT2","CPT3","CQT1","CQT2","NCT1","NCT2","NCT3","OAT1","OAT2","OAT3","OMT1","OMT2","OMT3","OPT1","OPT2","OPT3","CUT1","CUT2","CUT3","CVT1","CVT2","CVT3","CWT1","CWT2","CWT3","CXT1","CXT2","CXT3","NAT1","NAT2","NAT3","OUT1","OUT2","OUT3","OVT1","OVT2","OVT3","OXT1","OXT2","OXT3","OYT1","OYT2","OYT3","OZT1","OZT2","OZT3","PAT1","PAT2","PAT3")
norm_cts$Geneid <- NULL
cts_table$Geneid <- NULL
cts_table$Chr <- NULL
cts_table$Start <- NULL
cts_table$End <- NULL
cts_table$Strand <- NULL
cts_table$Length <- NULL
cts_table$LGE5 <- NULL
cts_table$OAE5 <- NULL
totals <- data.frame(sample = info$sample,total_raw_cts=colSums(cts_table))
info <- merge(info, totals, by = c("sample"))
totals <- data.frame(sample = info$sample,total_norm_cts=colSums(norm_cts))
info <- merge(info, totals, by = c("sample"))
info <- merge(info,map, by = c("sample"))
info <- info[info$sample %in% mis_indivs,]
head(info)
gc_meds <- c()
dups <- c()
reads <- c()
depths <- c()
tins <- c()


#i <- "PAT3"
for (i in info$sample)
{
  
  gc_file <- paste("C:/Users/jmcgirr/Documents/all_2018_samples/mapping_stats_qc/gc/",i, "_gc.GC.xls", sep = "")
  gc <- read.table(gc_file, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
  gc_meds <- c(gc_meds, median(gc$GC.))
  
  tin_file <- paste("C:/Users/jmcgirr/Documents/all_2018_samples/mapping_stats_qc/tin/",i, ".rna.sort.summary.txt", sep = "")
  tin <- read.table(tin_file, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
  tins <- c(tins, tin$TIN.median.)
  
  #dup_file <- paste("C:/Users/jmcgirr/Documents/all_2018_samples/mapping_stats_qc/dups/",i, "_dups.seq.DupRate.xls", sep = "")
  dup_file <- paste("C:/Users/jmcgirr/Documents/all_2018_samples/mapping_stats_qc/dups/",i, "_dups.pos.DupRate.xls", sep = "")
  dup <- read.table(dup_file, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
  head(dup)
  non_dup <- dup[1,2]
  dup <- dup[-1,]
  non_dup - sum(dup$UniqReadNumber)
  sum(dup$UniqReadNumber) / non_dup
  sum(dup$UniqReadNumber)
  non_dup + sum(dup$UniqReadNumber)
  dups <- c(dups, (sum(dup$UniqReadNumber)))
  
  depth_file <- paste("C:/Users/jmcgirr/Documents/all_2018_samples/mapping_stats_qc/coverage/",i, "_feature_coverage.txt", sep = "")
  depth <- read.table(depth_file, header = FALSE, stringsAsFactors = FALSE, sep = "\t")
  head(depth)
  depth <- mean(depth$V6)
  depths <- c(depths, depth)
  
}
info$median_gc_content <- gc_meds
info$duplicate_reads <- dups
info$prop_duplicate_reads <- info$duplicate_reads/info$Total
info$avg_depth_across_features <- depths
#info$raw_fastq_reads <- (info$raw_fastq_lines / 4) * 2
#info$prop_reads_mapped <- info$reads_mapped / info$raw_fastq_reads
info$tin_median <- tins
head(info)
range(info$avg_depth_across_features)
mean(info$avg_depth_across_features)
median(info$avg_depth_across_features)
#info$type_stage <- paste(info$cross_type, info$stage)
#write.table(info, "C:/Users/jmcgirr/Documents/all_2018_samples/mapping_stats_qc/all_qc_stats.txt", quote = FALSE, row.names = FALSE, sep = "\t")

cols <- c(red,"darkgrey",pur,gre,"black",dkr,blu)
cols <- c(rep("black",7))
unique(info$species_region_lake_cross_with_species)
#info <- info[which(info$stage == "8dpf"),]
box_plts <- info[order(info$species_region_lake_cross_with_species),] 

p1 <- ggboxplot(box_plts, "species_region_lake_cross_with_species", "total_norm_cts",color = cols, ylab = "total normalized counts\n", xlab="")+
  #stat_compare_means(method = "anova", label.y = max(info$total_norm_cts)+(max(info$total_norm_cts)*.05),label.x = 2.25)+      # Add global p-value
  scale_x_discrete(labels= c("generalists","new\nprovidence", "san salvador\nF1 hybrids","molluscivores", "north \ncarolina","outgroup\ngeneralist\nF1 hybrids", "scale-eaters"))+
  theme(axis.text=element_text(size=8.5))
  #stat_compare_means(label = "p.signif", method = "t.test",ref.group = ".all.", label.y = (max(info$total_norm_cts)+(max(info$total_norm_cts)*.01)))
#an <- aov(info$total_norm_cts~info$species_region_lake_cross_with_species)
#summary(an)
#tiff("C:/Users/jmcgirr/Documents/all_2018_samples/manuscript_figs/supp/qc/qc_norm_cts.tiff", width = 6.5, height = 3, units = 'in', res = 1000)
p1
#dev.off()
p1 <- ggboxplot(box_plts, "species_region_lake_cross_with_species", "percent_assigned",color = cols, ylab = "% reads assigned to features\n", xlab="")+
  #stat_compare_means(method = "anova", label.y = max(info$percent_assigned)+(max(info$percent_assigned)*.01),label.x = 2.25)+      # Add global p-value
  scale_x_discrete(labels= c("generalists","new\nprovidence", "san salvador\nF1 hybrids","molluscivores", "north \ncarolina","outgroup\ngeneralist\nF1 hybrids", "scale-eaters"))+
  theme(axis.text=element_text(size=8.5))
#stat_compare_means(label = "p.signif", method = "t.test",ref.group = ".all.", label.y = (max(info$percent_assigned)+(max(info$percent_assigned)*.01)), hide.ns = TRUE)
#an <- aov(info$percent_assigned~info$species_region_lake_cross_with_species)
#summary(an)
#tiff("C:/Users/jmcgirr/Documents/all_2018_samples/manuscript_figs/supp/qc/prop_reads_assigned_species.tiff", width = 6.5, height = 3, units = 'in', res = 1000)
p1
#dev.off()
p1 <- ggboxplot(box_plts, "species_region_lake_cross_with_species", "median_gc_content",color = cols, ylab = "median % GC content\n", xlab="")+
  #stat_compare_means(method = "anova", label.y = max(info$median_gc_content)+(max(info$median_gc_content)*.01),label.x = 2.25)+      # Add global p-value
  scale_x_discrete(labels= c("generalists","new\nprovidence", "san salvador\nF1 hybrids","molluscivores", "north \ncarolina","outgroup\ngeneralist\nF1 hybrids", "scale-eaters"))+
  theme(axis.text=element_text(size=8.5))
#stat_compare_means(label = "p.signif", method = "t.test",ref.group = ".all.", label.y = (max(info$median_gc_content)+(max(info$median_gc_content)*.01)), hide.ns = TRUE)
#an <- aov(info$median_gc_content~info$species_region_lake_cross_with_species)
#summary(an)
#tiff("C:/Users/jmcgirr/Documents/all_2018_samples/manuscript_figs/supp/qc/qc_gc_content.tiff", width =6.5, height = 3, units = 'in', res = 1000)
p1
#dev.off()
p1 <- ggboxplot(box_plts, "species_region_lake_cross_with_species", "avg_depth_across_features",color = cols, ylab = "mean depth across features\n", xlab="")+
  #stat_compare_means(method = "anova", label.y = max(info$avg_depth_across_features)+(max(info$avg_depth_across_features)*.05),label.x = 2.25)+      # Add global p-value
  scale_x_discrete(labels= c("generalists","new\nprovidence", "san salvador\nF1 hybrids","molluscivores", "north \ncarolina","outgroup\ngeneralist\nF1 hybrids", "scale-eaters"))+
  theme(axis.text=element_text(size=8.5))+
  scale_y_continuous(limits = c(0, 2000))
#stat_compare_means(label = "p.signif", method = "t.test",ref.group = ".all.", label.y = (max(info$avg_depth_across_features)+(max(info$avg_depth_across_features)*.01)))
#an <- aov(info$avg_depth_across_features~info$species_region_lake_cross_with_species)
#summary(an)
#tiff("C:/Users/jmcgirr/Documents/all_2018_samples/manuscript_figs/supp/qc/qc_depth.tiff", width = 6.5, height = 3, units = 'in', res = 1000)
p1
#dev.off()
p1 <- ggboxplot(box_plts, "species_region_lake_cross_with_species", "tin_median",color = cols, ylab = "median TIN\n", xlab="")+
  #stat_compare_means(method = "anova", label.y = max(info$tin_median)+(max(info$tin_median)*.05),label.x = 2.25)+      # Add global p-value
  scale_x_discrete(labels= c("generalists","new\nprovidence", "san salvador\nF1 hybrids","molluscivores", "north \ncarolina","outgroup\ngeneralist\nF1 hybrids", "scale-eaters"))+
  theme(axis.text=element_text(size=8.5))+
  stat_compare_means(label = "p.signif", method = "t.test",ref.group = ".all.", label.y = (max(info$tin_median)+(max(info$tin_median)*.01)))
an <- aov(info$tin_median~info$species_region_lake_cross_with_species)
TukeyHSD(an)
#summary(an)
#tiff("C:/Users/jmcgirr/Documents/all_2018_samples/manuscript_figs/supp/qc/qc_tin.tiff", width = 6.5, height = 3, units = 'in', res = 1000)
p1
#dev.off()
p1 <- ggboxplot(box_plts, "species_region_lake_cross_with_species", "prop_duplicate_reads",color = cols, ylab = "prop_duplicate_reads\n", xlab="species")+
  scale_y_continuous(labels = scales::scientific)+
  stat_compare_means(method = "anova", label.y = max(info$prop_duplicate_reads)+(max(info$prop_duplicate_reads)*.05),label.x = 2.25)+      # Add global p-value
  scale_x_discrete(labels= c("generalists","new\nprovidence", "san salvador\nF1 hybrids","molluscivores", "north \ncarolina","outgroup\ngeneralist\nF1 hybrids", "scale-eaters"))+
  theme(axis.text=element_text(size=8.5))+
  stat_compare_means(label = "p.signif", method = "t.test",ref.group = ".all.", label.y = (max(info$prop_duplicate_reads)+(max(info$prop_duplicate_reads)*.01)))
#an <- aov(info$tin_median~info$species_region_lake_cross_with_species)
#summary(an)
#tiff("C:/Users/jmcgirr/Documents/all_2018_samples/manuscript_figs/supp/qc/qc_dups.tiff", width = 6.5, height = 3, units = 'in', res = 1000)
p1
#dev.off()



#####
###### DNA ######
library(reshape2)
dna_inds <- c("CRPA1", "CRPA3", "CRPM1", "CRPM2", "CRPM3", "CRPP2", "CRPP3", "CRPPQ", "CUNA1", "ETAA1", "GREA1", "GREA2", "GREP1", "GREP2", "LILA1", "LILM3", "LILM4", "LILMQ", "LILP3", "LILP4", "LILPQ", "ME2A1", "ME2A2", "ME2P1", "MERA2", "MERA3", "MRKA1", "MRKM1", "MRKM2", "OSPA1", "OSPM1", "OSPM2", "OSPM3", "OSPP1", "OSPP2", "OSPP3", "OYSP1", "OYSP3", "PIGA1", "PIGA3", "MAY1", "SIM1", "DIAB", "OAF1", "CPF1", "NAF1", "CAM1", "OAM1", "CMF1", "CUNP", "OMF1", "OPM1", "CPM1", "OMM1", "OPF1", "CMM1", "CAF1", "CAM2")
depths <- c()
#i <- "CRPA1"
for (i in dna_inds)
{
  
  depth_file <- paste("C:/Users/jmcgirr/Documents/all_2018_samples/mapping_stats_qc/dna/",i, "_avg_coverage.txt", sep = "")
  depth <- read.table(depth_file, header = FALSE, stringsAsFactors = FALSE, sep = "\t")
  head(depth)
  depth <- cbind(depth, colsplit(depth$V1, "  ", c("j", "depth")))
  depths <- c(depths, depth$depth)
  
}
depths <- data.frame(sample = dna_inds ,
           avg_coverage = depths ,
           stringsAsFactors = FALSE)
range(depths$avg_coverage)
mean(depths$avg_coverage)

setwd("C:/Users/jmcgirr/Documents/all_2018_samples/mapping_stats_qc/dna/trim_reports/")
trim_files <- list.files(pattern="*_trimming_report.txt")
trims <- c()
#i <- trim_files[[1]]
for(i in trim_files)
{
  trim <- read.delim(i, header = FALSE, stringsAsFactors = FALSE, sep = "\t")
  trim <- strsplit(trim[20,], ":     ")[[1]][2]
  trim <- strsplit(trim, " \\(")[[1]][1]
  trim <- as.numeric(gsub(",", "", trim))
  trims <- c(trims, trim)
}
sum(trims)
#####


