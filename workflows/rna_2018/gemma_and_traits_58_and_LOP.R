##### gemma result files 58 #####
library(reshape2)

trait_name <- "up_jaw"

al  <- read.table(paste("C:/Users/jmcgirr/Documents/all_2018_samples/GEMMA/output_58/mean_alphas_",trait_name,"_summed_20kb.txt", sep = ""), header=F, stringsAsFactors = F)
be  <- read.table(paste("C:/Users/jmcgirr/Documents/all_2018_samples/GEMMA/output_58/mean_betas_",trait_name,"_summed_20kb.txt", sep = ""), header=F, stringsAsFactors = F)
ga  <- read.table(paste("C:/Users/jmcgirr/Documents/all_2018_samples/GEMMA/output_58/mean_gammas_",trait_name,"_summed_20kb.txt", sep = ""), header=F, stringsAsFactors = F)
head(al)

al <- cbind(al, colsplit(al$V3, "\\|", c("wind_end", "alpha")))
be <- cbind(be, colsplit(be$V3, "\\|", c("wind_end", "beta")))
ga <- cbind(ga, colsplit(ga$V3, "\\|", c("wind_end", "gamma")))
param <- merge(al,be, by =c("V1","V2", "wind_end"))
ga <- ga[c("V1","V2", "wind_end", "gamma")]
param <- merge(param, ga, by = c("V1","V2", "wind_end"))
param <- param[c("V1","V2", "wind_end","alpha","beta", "gamma")]

param <- param[order(param$gamma, decreasing = TRUE),]
head(param)
param <- param[order(param$beta, decreasing = TRUE),]
head(param)
param <- param[order(param$beta, decreasing = FALSE),]
head(param)

quantile(param$gamma,.99)

candidates_99 <- param[which(param$gamma >= 0.00175),]
#candidates_95 <- param[which(param$gamma >= 0.00118),]
candidates_99 <- candidates_99[c("V1","V2","wind_end", "gamma")]
#write.table(candidates_99,"C:/Users/jmcgirr/Documents/all_2018_samples/GEMMA/candidates_low_jaw_99percentile.bed", col.names = FALSE, row.names = FALSE, quote = FALSE)
#bedtools intersect -a candidates_low_jaw_99percentile.sort.bed  -b Cyprinodon_NW_noheader_mrna.sort.bed -wa -wb > candidates_low_jaw_99percentile.genes.txt


param <- param[order(param$V1, param$wind_end),]
head(param)
param$x_ax <- c(1:nrow(param))
param$effect_size <- param$beta * param$gamma
col.topo <- c("black","darkgrey")
palette(col.topo)
family <- as.factor(param$V1)
png("C:/Users/jmcgirr/Documents/all_2018_samples/manuscript_figs/supp/gemma_pip.png", width = 9, height = 4, units = 'in', res = 1000)
plot(param$x_ax, param$gamma, col = family, pch = 16, 
     xaxt = 'n', ylab="PIP", xlab = "scaffolds", yaxt = "n", ylim = c(0,0.11))
axis(side=2, at=c(0,0.05,0.1), labels = TRUE)
abline(h=0.00175, lty = 2, col = red)
dev.off()


candidates_95_up <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/GEMMA/candidates_95_percentile/candidates_up_jaw_95percentile.genes.txt", header=F, stringsAsFactors = F)
candidates_95_low <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/GEMMA/candidates_95_percentile/candidates_low_jaw_95percentile.genes.txt", header=F, stringsAsFactors = F)
nrow(candidates_95_up)
nrow(candidates_95_low)
length(intersect(candidates_95_up$V4, candidates_95_low$V4))
length(setdiff(candidates_95_up$V4, candidates_95_low$V4))
length(setdiff(candidates_95_low$V4, candidates_95_up$V4))


iters <- c(2:10)
hyp1 <- read.table(paste("C:/Users/jmcgirr/Documents/all_2018_samples/GEMMA/output_58/gemma_",trait_name,
                        "_1.hyp.txt",sep = ""), header=T, stringsAsFactors = F)
png("C:/Users/jmcgirr/Documents/all_2018_samples/manuscript_figs/supp/gemma_pve.png", width = 4, height = 4, units = 'in', res = 1000)
plot(density(hyp1$pve), xlab="PVE", ylab="density", main=NA, lty = 1, lwd = 2, 
     col = blu, cex.axis = 1.5, cex.lab = 1.5, ylim = c(0,80))
for(i in iters)
{
  hyp <- read.table(paste("C:/Users/jmcgirr/Documents/all_2018_samples/GEMMA/output_58/gemma_",trait_name,
                          "_", i,".hyp.txt",sep = ""), header=T, stringsAsFactors = F) 
  lines(density(hyp$pve), col = blu)
  
}
dev.off()
png("C:/Users/jmcgirr/Documents/all_2018_samples/manuscript_figs/supp/gemma_pge.png", width = 4, height = 4, units = 'in', res = 1000)
plot(density(hyp1$pge), xlab="PGE", ylab="density", main=NA, lty = 1, lwd = 2, col = blu, cex.axis = 1.5, cex.lab = 1.5)
for(i in iters)
{
  hyp <- read.table(paste("C:/Users/jmcgirr/Documents/all_2018_samples/GEMMA/output_58/gemma_",trait_name,
                          "_", i,".hyp.txt",sep = ""), header=T, stringsAsFactors = F) 
  lines(density(hyp$pge), col = blu)
  
}
dev.off()
png("C:/Users/jmcgirr/Documents/all_2018_samples/manuscript_figs/supp/gemma_nsnp.png", width = 4, height = 4, units = 'in', res = 1000)
plot(density(hyp1$n_gamma), xlab="n-SNP", ylab="density", main=NA, lty = 1, lwd = 2, col = blu, cex.axis = 1.5, cex.lab = 1.5)
for(i in iters)
{
  hyp <- read.table(paste("C:/Users/jmcgirr/Documents/all_2018_samples/GEMMA/output_58/gemma_",trait_name,
                          "_", i,".hyp.txt",sep = ""), header=T, stringsAsFactors = F) 
  lines(density(hyp$n_gamma), col = blu)
  
}
dev.off()

head(hyp)
plot(density(hyp$pve), xlab="PVE", ylab="density", main=NA, lty = 1, lwd = 2, col = "blue4", cex.axis = 1.5, cex.lab = 1.5)
plot(density(hyp$pge), xlab="PGE", ylab="density", main=NA, lty = 1, lwd = 2, col = "blue4", cex.axis = 1.5, cex.lab = 1.5)
plot(density(hyp$n_gamma), xlab="n-SNP", ylab="density", main=NA, lty = 1, lwd = 2, col = "blue4", cex.axis = 1.5, cex.lab = 1.5)

##### overlap gemma and tajimas d #####
require(data.table)

gem <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/GEMMA/candidates_95_percentile/candidates_up_jaw_95percentile.genes.txt", header = FALSE, stringsAsFactors = FALSE)
head(gem)
colnames(gem) <- c("CHROM","START","END","mrna","gene","zeb_gene")
b <- gem[c("CHROM", "START", "END")]

comps <- c("")
for (i in  c(1:length(comps)))
{

taj <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/fst_dna/a_taj_d_20kb.txt.Tajima.D", header = TRUE, stringsAsFactors = FALSE)
head(taj)
taj$END <- taj$BIN_START + 19999
taj$START <- taj$BIN_START
a <- taj[c("CHROM", "START", "END")]#,"INDEX", "TajimaD")]

setDT(taj)
c <- foverlaps(taj, gem,by.x = c("CHROM", "START", "END"),by.y = c("CHROM", "START", "END"),type="within", nomatch=0L)
c <- as.data.frame(c)
rownames(c) <- NULL
thresh <- quantile(c$TajimaD, 0.1)
c <- c[which(c$TajimaD <= thresh),]
head(c)
}


#####

#####
##### gemma result files 94 #####
library(reshape2)

trait_name <- "low_jaw"

al  <- read.table(paste("D:/Martin Lab/lots_of_pups_project/GEMMA/output/mean_alphas_",trait_name,"_summed_20kb.txt", sep = ""), header=F, stringsAsFactors = F)
be  <- read.table(paste("D:/Martin Lab/lots_of_pups_project/GEMMA/output/mean_betas_",trait_name,"_summed_20kb.txt", sep = ""), header=F, stringsAsFactors = F)
ga  <- read.table(paste("D:/Martin Lab/lots_of_pups_project/GEMMA/output/mean_gammas_",trait_name,"_summed_20kb.txt", sep = ""), header=F, stringsAsFactors = F)
head(al)

al <- cbind(al, colsplit(al$V3, "\\|", c("wind_end", "alpha")))
be <- cbind(be, colsplit(be$V3, "\\|", c("wind_end", "beta")))
ga <- cbind(ga, colsplit(ga$V3, "\\|", c("wind_end", "gamma")))
param <- merge(al,be, by =c("V1","V2", "wind_end"))
ga <- ga[c("V1","V2", "wind_end", "gamma")]
param <- merge(param, ga, by = c("V1","V2", "wind_end"))
param <- param[c("V1","V2", "wind_end","alpha","beta", "gamma")]

param <- param[order(param$gamma, decreasing = TRUE),]
head(param)
param <- param[order(param$beta, decreasing = TRUE),]
head(param)
param <- param[order(param$beta, decreasing = FALSE),]
head(param)
quantile(param$gamma,.95)

candidates_95 <- param[which(param$gamma >= 0.00018),]
candidates_95 <- candidates_95[c("V1","V2","wind_end", "gamma")]
#write.table(candidates_95,"D:/Martin Lab/lots_of_pups_project/GEMMA/candidates_low_jaw_95percentile.bed", col.names = FALSE, row.names = FALSE, quote = FALSE)

trait_name <- "low_jaw"
iters <- c(2:10)
hyp1 <- read.table(paste("D:/Martin Lab/lots_of_pups_project/GEMMA/output/gemma_",trait_name,
                         "_94_1.hyp.txt",sep = ""), header=T, stringsAsFactors = F)
plot(density(hyp1$pve), xlab="PVE", ylab="density", main=NA, lty = 1, lwd = 2, col = "blue4", cex.axis = 1.5, cex.lab = 1.5)
for(i in iters)
{
  hyp <- read.table(paste("D:/Martin Lab/lots_of_pups_project/GEMMA/output/gemma_",trait_name,
                          "_94_", i,".hyp.txt",sep = ""), header=T, stringsAsFactors = F) 
  lines(density(hyp$pve), col = "blue4")
  
}
plot(density(hyp1$pge), xlab="PGE", ylab="density", main=NA, lty = 1, lwd = 2, col = "blue4", cex.axis = 1.5, cex.lab = 1.5)
for(i in iters)
{
  hyp <- read.table(paste("D:/Martin Lab/lots_of_pups_project/GEMMA/output/gemma_",trait_name,
                          "_94_", i,".hyp.txt",sep = ""), header=T, stringsAsFactors = F) 
  lines(density(hyp$pge), col = "blue4")
  
}
plot(density(hyp1$n_gamma), xlab="n-SNP", ylab="density", main=NA, lty = 1, lwd = 2, col = "blue4", cex.axis = 1.5, cex.lab = 1.5)
for(i in iters)
{
  hyp <- read.table(paste("D:/Martin Lab/lots_of_pups_project/GEMMA/output/gemma_",trait_name,
                          "_94_", i,".hyp.txt",sep = ""), header=T, stringsAsFactors = F) 
  lines(density(hyp$n_gamma), col = "blue4")
  
}



#####
####### traits for 58 pups project #####
traits <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/GEMMA/pupfish_measurements.txt", header=T, stringsAsFactors = F)

traits$log_std_len <- log(traits$standard_length)
traits$log_up_jaw <- log(traits$upperjawlength)
traits$log_low_jaw <- log(traits$lowerjawlength)
up_jaw_lm <- lm(traits$log_up_jaw~traits$log_std_len)
low_jaw_lm <- lm(traits$log_low_jaw~traits$log_std_len)
traits$resid_up_jaw <- data.frame(up_jaw_lm$residuals)[,1]
traits$resid_low_jaw <- data.frame(low_jaw_lm$residuals)[,1]
head(traits)

plot(traits$log_std_len, traits$log_low_jaw, col= c("red", "forestgreen", "blue","lightblue")[as.numeric(traits$speciesID)], 
     xlab = "log body length", ylab = "log lower jaw length", main = "Standardized Jaw Measurements", cex = 3, pch = 16)
abline(low_jaw_lm)
plot(traits$log_std_len, traits$log_up_jaw, col= c("red", "forestgreen", "blue","lightblue")[as.numeric(traits$speciesID)], 
     xlab = "log body length", ylab = "log upper jaw length", main = "Standardized Jaw Measurements", cex = 3, pch = 16)
abline(up_jaw_lm)

#write.table(traits,"C:/Users/jmcgirr/Documents/all_2018_samples/GEMMA/jaw_residuals.txt", row.names = FALSE, quote = FALSE)

fam <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/GEMMA/58_dna_maf_0.5_maxmiss_0.9_snps.recode.vcf.plink.fam", header=F, stringsAsFactors = F)
fam$sample <- fam$V1
fam <- merge(fam, traits,all = TRUE, by = c("sample"))
fam <- fam[c("sample","V1","V3","V4","V5","V6","resid_up_jaw","resid_low_jaw")]
fam
#write.table(fam,"C:/Users/jmcgirr/Documents/all_2018_samples/GEMMA/58_dna_maf_0.5_maxmiss_0.9_snps.recode.vcf.plink.fam",col.names = FALSE, row.names = FALSE, quote = FALSE)
#####
##### traits for lots of pups project ####

h <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/GEMMA/lop/mar_2019/all_pupfish_SSI_filtered.Q20.MAF0.05.MAXMISS0.9.Carib_filtered.Q20.MAF0.05.MAXMISS0.8.fam.useable", header=F, stringsAsFactors = F)
h$vcf_label <- h$V1
k <- read.table("D:/Martin Lab/lots_of_pups_project/measures_matched_to_vcf.txt", header=T, stringsAsFactors = F, sep = "\t")
c <- read.table("D:/Martin Lab/lots_of_pups_project/match_measures_caudal_lightness.txt", header=T, stringsAsFactors = F, sep = "\t")
setdiff(c$vcf_label, k$vcf_label)
l <- merge(c,k,all = TRUE, by = c("vcf_label"))
l <- merge(h,c,all.x = TRUE, by = c("vcf_label"))
huh <- c("GEOA11","GEOA10","GEOA7","GEOA6","OSPM10","OSPM1","OSPM2","MRKM5","GREA1","OSPM8","OSPM3","OSPM9","CRPP5","GREA2","OSPM7","OSPP3","CRPM3","CRPP9","CRPP2","OSPM6","OSPP5","LILM5","CRPP4","OYSM2","OSPP1","GEO2A9","GEO2A6","CRPP8","RHPA1","GEO2A8","LILM3","OSPP10","OSPP6","LILM4","OSPP11","GEO2A7","LILP4","GEO2A10","OSPP7","LILP5","GEO2A4","GEO2A1","OSPA9","OYSM8","CRPM10","GEO2A5","OYSM4","CRPM6","GNYA1","OSPA13","OYSP6","CRPM8","MRKM1","OSPA8","LILP3","CRPM9","NLLA1","MERA3","CRPM2","OYSM6","CRPM11","GEO2A3","CRPM5","CRPP3","GEO2A2","CRPM7","OSPA12","OYSM3","LILP-QTL","CRPA1","CRPA3","LILM-QTL","LILA1","BAVA11","BAVA12","BAVA13","BAVA14","BAVA2","BAVA4","BAVA5","BAVA6","BAVA7","BAVA8","CLRA1","CRPA1000","CRPA1003","CRPM1","CRPM1000","CRPM1001","CRPP-QTL","CRPP1000","CRPP1001","CRPP7","CUNA1","CUNA10","CUNA2","CUNA3","CUNA4","CUNA5","CUNA6","CUNA8","CUNA9","CUNP10","CUNP3","CUNP4","CUNP5","CUNP6","CUNP7","GEOA1","GEOA2","GEOA4","KILA1","LGIA1","LILM4-2","LILP4-2","ME2A1","ME2A2","ME2P1","ME2P1-1","MERA2","MRKA1","MRKM2","MRKM2-2","MRKM4","MRKM4-2","NCCA1","NCCA10","NCCA1000","NCCA12","NCCA15","NCCA2","NCCA3","NCCA5","NCCA6","NCCA7","NCCA8","NCCA9","OSPA1","OSPA10","OSPA1000","OSPA1001","OSPA11","OSPA4","OSPA5","OSPA6","OSPM1000","OSPM1001","OSPM11","OSPM4","OSPM5","OSPP1000","OSPP1001","OSPP9","OYSA1","OYSA2","OYSM1","OYSM5","OYSM7","PAIA1","VENA1","VENA10","VENA11-2","VENA12","VENA13","VENA2","VENA3","VENA4","VENA5","VENA8","VENA9","WDPA1","ACKA1","APHFAS","ARTA1","ARTA2","BAVA10","BAVA9","BONA1","CAIA1","CATA1","CRPA1001","CUNA7","CUNP11","CURA1","CURA2","CURA21","CYPALBI","CYPERE","CYPEXI","CYPFONTCARB","CYPLONGI","CYPMACCOCH17.4","CYPMACRO","CYPRAD18.1","CYPVER","DEAA1","EPLA1","ETA1","EXUA1","EXUA2","FCTA1","FLSA1","GREP1","GREP2","GRES3","GRES4","LILH1","LILH2","LILS1","LILS2","LILS3","LILS4","MAFA1","MAY1","MEG-Q1","MERP1","MERP2","MRKM3","NBIA1","NCCA11","NCCA4","OSPA7","OSPH1","OSPP2","OSPP4","OSPP8","OSPS1","OSPS10","OSPS11","OSPS2","OSPS3","OSPS4","OSPS5","OSPS6","OSPS7","OSPS8","OSPS9","OYSP1","OYSP3","OYSP4","OYSP7","OYSS3","OYSS4","OYSS5","OYSS6","OYSS9","PIGA1","PIGA3","PWLA1","SALA1","SCLA1","SIM1","SPPA1","VENA12.2","VENA7","WDPM2","CUATESS1","VENA1-2")
duplicated(huh)
#write.table(l, "C:/Users/jmcgirr/Documents/all_2018_samples/GEMMA/lop/mar_2019/test.txt", quote=FALSE, row.names = FALSE, sep = "\t")

l <- l[which(l$V7 == "YES" & is.na(l$lowerjawlength)),]
#write.table(l, "C:/Users/jmcgirr/Documents/all_2018_samples/GEMMA/lop/mar_2019/need_measures.txt", quote=FALSE, row.names = FALSE, sep = "\t")


traits <- read.table("D:/Martin Lab/lots_of_pups_project/measures_matched_to_vcf_nose_pigment_jaw.txt", header=T, stringsAsFactors = F, sep = "\t")
nrow(traits)
traits$log_std_len <- log(traits$standard_length)
traits$log_low_jaw <- log(traits$lowerjawlength)
traits$log_nose_length <- log((traits$nose_length +1))
traits$log_nose_height <- log(traits$nose_height)
low_jaw_lm <- lm(traits$log_low_jaw~traits$log_std_len)
nose_length_lm  <- lm(traits$log_nose_length~traits$log_std_len)
nose_height_lm  <- lm(traits$log_nose_height~traits$log_std_len)
traits$resid_low_jaw <- data.frame(low_jaw_lm$residuals)[,1]
traits[names(nose_length_lm$residuals),"resid_nose_length"]<-nose_length_lm$residuals
traits[names(nose_height_lm$residuals),"resid_nose_height"]<-nose_height_lm$residuals

head(traits)

pdf("D:/Martin Lab/lots_of_pups_project/transformed_traits.pdf",width=6,height=6)
plot(traits$log_std_len, traits$log_low_jaw, col= c("red", "forestgreen", "blue")[as.factor(traits$species)], 
     xlab = "log body length", ylab = "log lower jaw length", main = "Standardized Jaw Measurements", cex = 3, pch = 16)
abline(low_jaw_lm)
plot(traits$log_std_len, traits$log_nose_length, col= c("red", "forestgreen", "blue","lightblue")[as.factor(traits$species)], 
     xlab = "log body length", ylab = "log nose length", main = "Standardized Nose Measurements", cex = 3, pch = 16)
abline(nose_length_lm)
plot(traits$log_std_len, traits$log_nose_height, col= c("red", "forestgreen", "blue","lightblue")[as.factor(traits$species)], 
     xlab = "log body length", ylab = "log nose height", main = "Standardized Nose Measurements", cex = 3, pch = 16)
abline(nose_height_lm)
dev.off()

fam <- read.table("D:/Martin Lab/lots_of_pups_project/GEMMA/all_pupfish_SSI_filtered.Q20.MAF0.05.MAXMISS0.9.Carib_filtered.Q20.MAF0.05.MAXMISS0.5.fam.original", header=FALSE, stringsAsFactors = F)
fam$vcf_label <- fam$V1
fam <- merge(fam, traits,all = TRUE, by = c("vcf_label"))
fam <- fam[c("vcf_label","V1","V3","V4","V5","V6","caudal_lightness_average","resid_low_jaw", "resid_nose_length", "resid_nose_height")]
# -n 2 = pigment 3 = low jaw 4 = nose length 5 = nose height
#write.table(fam,"D:/Martin Lab/lots_of_pups_project/GEMMA/all_pupfish_SSI_filtered.Q20.MAF0.05.MAXMISS0.9.Carib_filtered.Q20.MAF0.05.MAXMISS0.5.fam",col.names = FALSE, row.names = FALSE, quote = FALSE)

traits[which(traits$log_nose_length < -4),]




#####
##### overlap with MBE results #####

mbe <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/GEMMA/mp_gemma_mbe.txt", header=T, stringsAsFactors = F)
scaff_key <- read.table("D:/Cyprinodon/scaffold_conversion.txt", header = TRUE, stringsAsFactors = FALSE)
head(mbe)
head(scaff_key)
mbe <- merge(mbe, scaff_key, by = c("scaffold"))
mbe <- mbe[c("chr","strt","stp")]
#write.table(mbe,"C:/Users/jmcgirr/Documents/all_2018_samples/GEMMA/candidates_31_mbe.bed", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

candidates_95_up <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/GEMMA/candidates_95_percentile/candidates_up_jaw_95percentile.genes.txt", header=F, stringsAsFactors = F)
candidates_95_low <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/GEMMA/candidates_95_percentile/candidates_low_jaw_95percentile.genes.txt", header=F, stringsAsFactors = F)
nrow(candidates_95_up)
nrow(candidates_95_low)
length(intersect(candidates_95_up$V4, candidates_95_low$V4))
length(setdiff(candidates_95_up$V4, candidates_95_low$V4))
length(setdiff(candidates_95_low$V4, candidates_95_up$V4))

#####