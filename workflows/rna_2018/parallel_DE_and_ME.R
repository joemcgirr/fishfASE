#load DE/ME
{
  DE_genes_dir <- "C:/Users/jmcgirr/Documents/all_2018_samples/conditions/"
  venn_dir <- "C:/Users/jmcgirr/Documents/all_2018_samples/mse_venns/"
  stage <- "8dpf"
  # load misexpressed
  {
    # am crosses
    caxcm <- read.csv(paste(DE_genes_dir,"DE_CRPA_and_CRPM_vs_CAxCM_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    cmxca <- read.csv(paste(DE_genes_dir,"DE_CRPM_and_CRPA_vs_CMxCA_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    oaxom <- read.csv(paste(DE_genes_dir,"DE_OSPA_and_OSPM_vs_OAxOM_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    omxoa <- read.csv(paste(DE_genes_dir,"DE_OSPM_and_OSPA_vs_OMxOA_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    # ap crosses
    caxcp <- read.csv(paste(DE_genes_dir,"DE_CRPA_and_CRPP_vs_CAxCP_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    oaxop <- read.csv(paste(DE_genes_dir,"DE_OSPA_and_OSPP_vs_OAxOP_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    opxoa <- read.csv(paste(DE_genes_dir,"DE_OSPP_and_OSPA_vs_OPxOA_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    # specialist crosses
    cmxcp <- read.csv(paste(DE_genes_dir,"DE_CRPM_and_CRPP_vs_CMxCP_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    opxom <- read.csv(paste(DE_genes_dir,"DE_OSPM_and_OSPP_vs_OPxOM_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    # outgroup crosses
    naxca <-  read.csv(paste(DE_genes_dir,"DE_CRPA_and_NCA_vs_NAxCA_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    upxca <-  read.csv(paste(DE_genes_dir,"DE_CRPA_and_UPxUA_vs_UPxCA_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    # lake crosses
    oaxca <- read.csv(paste(DE_genes_dir,"DE_CRPA_and_OSPA_vs_OAxCA_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    caxoa <- read.csv(paste(DE_genes_dir,"DE_CRPA_and_OSPA_vs_CAxOA_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    cmxom <- read.csv(paste(DE_genes_dir,"DE_CRPM_and_OSPM_vs_CMxOM_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    opxcp <- read.csv(paste(DE_genes_dir,"DE_CRPP_and_OSPP_vs_OPxCP_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    
    # sig misexpressed
    caxcm <- caxcm[which(caxcm$padj <= 0.05),] 
    cmxca <- cmxca[which(cmxca$padj <= 0.05),] 
    oaxom <- oaxom[which(oaxom$padj <= 0.05),] 
    omxoa <- omxoa[which(omxoa$padj <= 0.05),] 
    caxcp <- caxcp[which(caxcp$padj <= 0.05),] 
    oaxop <- oaxop[which(oaxop$padj <= 0.05),] 
    opxoa <- opxoa[which(opxoa$padj <= 0.05),] 
    cmxcp <- cmxcp[which(cmxcp$padj <= 0.05),] 
    opxom <- opxom[which(opxom$padj <= 0.05),] 
    naxca <- naxca[which(naxca$padj <= 0.05),] 
    upxca <- upxca[which(upxca$padj <= 0.05),] 
    oaxca <- oaxca[which(oaxca$padj <= 0.05),] 
    caxoa <- caxoa[which(caxoa$padj <= 0.05),] 
    cmxom <- cmxom[which(cmxom$padj <= 0.05),] 
    opxcp <- opxcp[which(opxcp$padj <= 0.05),] 
    
    caxcm_8_ME <- caxcm$related_accession
    cmxca_8_ME <- cmxca$related_accession
    oaxom_8_ME <- oaxom$related_accession
    omxoa_8_ME <- omxoa$related_accession
    caxcp_8_ME <- caxcp$related_accession
    oaxop_8_ME <- oaxop$related_accession
    opxoa_8_ME <- opxoa$related_accession
    cmxcp_8_ME <- cmxcp$related_accession
    opxom_8_ME <- opxom$related_accession
    naxca_8_ME <- naxca$related_accession
    upxca_8_ME <- upxca$related_accession
    oaxca_8_ME <- oaxca$related_accession
    caxoa_8_ME <- caxoa$related_accession
    cmxom_8_ME <- cmxom$related_accession
    opxcp_8_ME <- opxcp$related_accession
  }
  # load DE
  {
    # am crosses
    caxcm <- read.csv(paste(DE_genes_dir,"DE_CRPA_vs_CRPM_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    oaxom <- read.csv(paste(DE_genes_dir,"DE_OSPA_vs_OSPM_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    # ap crosses
    caxcp <- read.csv(paste(DE_genes_dir,"DE_CRPA_vs_CRPP_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    oaxop <- read.csv(paste(DE_genes_dir,"DE_OSPA_vs_OSPP_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    # specialist crosses
    cmxcp <- read.csv(paste(DE_genes_dir,"DE_CRPM_vs_CRPP_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    omxop <- read.csv(paste(DE_genes_dir,"DE_OSPM_vs_OSPP_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    # outgroup crosses
    caxna <-  read.csv(paste(DE_genes_dir,"DE_CRPA_vs_NCA_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    oaxna <-  read.csv(paste(DE_genes_dir,"DE_OSPA_vs_NCA_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    caxup <-  read.csv(paste(DE_genes_dir,"DE_CRPA_vs_UPxUA_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    oaxup <-  read.csv(paste(DE_genes_dir,"DE_OSPA_vs_UPxUA_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    
    # lake crosses
    caxoa <- read.csv(paste(DE_genes_dir,"DE_CRPA_vs_OSPA_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    cmxom <- read.csv(paste(DE_genes_dir,"DE_CRPM_vs_OSPM_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    cpxop <- read.csv(paste(DE_genes_dir,"DE_CRPP_vs_OSPP_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    
    # number of informitive transcripts
    caxcm_n_8_DE <- nrow(caxcm)
    oaxom_n_8_DE <- nrow(oaxom)
    caxcp_n_8_DE <- nrow(caxcp)
    oaxop_n_8_DE <- nrow(oaxop)
    cmxcp_n_8_DE <- nrow(cmxcp)
    omxop_n_8_DE <- nrow(omxop)
    caxna_n_8_DE <- nrow(caxna)
    oaxna_n_8_DE <- nrow(oaxna)
    caxup_n_8_DE <- nrow(caxup)
    oaxup_n_8_DE <- nrow(oaxup)
    caxoa_n_8_DE <- nrow(caxoa)
    cmxom_n_8_DE <- nrow(cmxom)
    cpxop_n_8_DE <- nrow(cpxop)
    
    # sig DE
    caxcm <- caxcm[which(caxcm$padj <= 0.05),] 
    oaxom <- oaxom[which(oaxom$padj <= 0.05),]
    caxcp <- caxcp[which(caxcp$padj <= 0.05),]
    oaxop <- oaxop[which(oaxop$padj <= 0.05),]
    cmxcp <- cmxcp[which(cmxcp$padj <= 0.05),]
    omxop <- omxop[which(omxop$padj <= 0.05),]
    caxna <- caxna[which(caxna$padj <= 0.05),]
    oaxna <- oaxna[which(oaxna$padj <= 0.05),]
    caxup <- caxup[which(caxup$padj <= 0.05),]
    oaxup <- oaxup[which(oaxup$padj <= 0.05),]
    caxoa <- caxoa[which(caxoa$padj <= 0.05),]
    cmxom <- cmxom[which(cmxom$padj <= 0.05),]
    cpxop <- cpxop[which(cpxop$padj <= 0.05),]
    
    caxcm_8_DE <- caxcm$related_accession
    oaxom_8_DE <- oaxom$related_accession
    caxcp_8_DE <- caxcp$related_accession
    oaxop_8_DE <- oaxop$related_accession
    cmxcp_8_DE <- cmxcp$related_accession
    omxop_8_DE <- omxop$related_accession
    caxna_8_DE <- caxna$related_accession
    oaxna_8_DE <- oaxna$related_accession
    caxup_8_DE <- caxup$related_accession
    oaxup_8_DE <- oaxup$related_accession
    caxoa_8_DE <- caxoa$related_accession
    cmxom_8_DE <- cmxom$related_accession
    cpxop_8_DE <- cpxop$related_accession
  }
  stage <- "48hpf"
  # load misexpressed
  {
    # am crosses
    caxcm <- read.csv(paste(DE_genes_dir,"DE_CRPA_and_CRPM_vs_CAxCM_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    cmxca <- read.csv(paste(DE_genes_dir,"DE_CRPM_and_CRPA_vs_CMxCA_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    oaxom <- read.csv(paste(DE_genes_dir,"DE_OSPA_and_OSPM_vs_OAxOM_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    omxoa <- read.csv(paste(DE_genes_dir,"DE_OSPM_and_OSPA_vs_OMxOA_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    # ap crosses
    caxcp <- read.csv(paste(DE_genes_dir,"DE_CRPA_and_CRPP_vs_CAxCP_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    oaxop <- read.csv(paste(DE_genes_dir,"DE_OSPA_and_OSPP_vs_OAxOP_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    opxoa <- read.csv(paste(DE_genes_dir,"DE_OSPP_and_OSPA_vs_OPxOA_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    # specialist crosses
    cmxcp <- read.csv(paste(DE_genes_dir,"DE_CRPM_and_CRPP_vs_CMxCP_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    opxom <- read.csv(paste(DE_genes_dir,"DE_OSPM_and_OSPP_vs_OPxOM_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    # outgroup crosses
    naxca <-  read.csv(paste(DE_genes_dir,"DE_CRPA_and_NCA_vs_NAxCA_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    upxca <-  read.csv(paste(DE_genes_dir,"DE_CRPA_and_UPxUA_vs_UPxCA_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    # lake crosses
    oaxca <- read.csv(paste(DE_genes_dir,"DE_CRPA_and_OSPA_vs_OAxCA_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    caxoa <- read.csv(paste(DE_genes_dir,"DE_CRPA_and_OSPA_vs_CAxOA_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    cmxom <- read.csv(paste(DE_genes_dir,"DE_CRPM_and_OSPM_vs_CMxOM_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    opxcp <- read.csv(paste(DE_genes_dir,"DE_CRPP_and_OSPP_vs_OPxCP_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    
    # sig misexpressed
    caxcm <- caxcm[which(caxcm$padj <= 0.05),] 
    cmxca <- cmxca[which(cmxca$padj <= 0.05),] 
    oaxom <- oaxom[which(oaxom$padj <= 0.05),] 
    omxoa <- omxoa[which(omxoa$padj <= 0.05),] 
    caxcp <- caxcp[which(caxcp$padj <= 0.05),] 
    oaxop <- oaxop[which(oaxop$padj <= 0.05),] 
    opxoa <- opxoa[which(opxoa$padj <= 0.05),] 
    cmxcp <- cmxcp[which(cmxcp$padj <= 0.05),] 
    opxom <- opxom[which(opxom$padj <= 0.05),] 
    naxca <- naxca[which(naxca$padj <= 0.05),] 
    upxca <- upxca[which(upxca$padj <= 0.05),] 
    oaxca <- oaxca[which(oaxca$padj <= 0.05),] 
    caxoa <- caxoa[which(caxoa$padj <= 0.05),] 
    cmxom <- cmxom[which(cmxom$padj <= 0.05),] 
    opxcp <- opxcp[which(opxcp$padj <= 0.05),] 
    
    caxcm_48_ME <- caxcm$related_accession
    cmxca_48_ME <- cmxca$related_accession
    oaxom_48_ME <- oaxom$related_accession
    omxoa_48_ME <- omxoa$related_accession
    caxcp_48_ME <- caxcp$related_accession
    oaxop_48_ME <- oaxop$related_accession
    opxoa_48_ME <- opxoa$related_accession
    cmxcp_48_ME <- cmxcp$related_accession
    opxom_48_ME <- opxom$related_accession
    naxca_48_ME <- naxca$related_accession
    upxca_48_ME <- upxca$related_accession
    oaxca_48_ME <- oaxca$related_accession
    caxoa_48_ME <- caxoa$related_accession
    cmxom_48_ME <- cmxom$related_accession
    opxcp_48_ME <- opxcp$related_accession
  }
  # load DE
  {
    # am crosses
    caxcm <- read.csv(paste(DE_genes_dir,"DE_CRPA_vs_CRPM_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    oaxom <- read.csv(paste(DE_genes_dir,"DE_OSPA_vs_OSPM_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    # ap crosses
    caxcp <- read.csv(paste(DE_genes_dir,"DE_CRPA_vs_CRPP_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    oaxop <- read.csv(paste(DE_genes_dir,"DE_OSPA_vs_OSPP_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    # specialist crosses
    cmxcp <- read.csv(paste(DE_genes_dir,"DE_CRPM_vs_CRPP_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    omxop <- read.csv(paste(DE_genes_dir,"DE_OSPM_vs_OSPP_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    # outgroup crosses
    caxna <-  read.csv(paste(DE_genes_dir,"DE_CRPA_vs_NCA_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    oaxna <-  read.csv(paste(DE_genes_dir,"DE_OSPA_vs_NCA_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    caxup <-  read.csv(paste(DE_genes_dir,"DE_CRPA_vs_UPxUA_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    oaxup <-  read.csv(paste(DE_genes_dir,"DE_OSPA_vs_UPxUA_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    
    # lake crosses
    caxoa <- read.csv(paste(DE_genes_dir,"DE_CRPA_vs_OSPA_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    cmxom <- read.csv(paste(DE_genes_dir,"DE_CRPM_vs_OSPM_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    cpxop <- read.csv(paste(DE_genes_dir,"DE_CRPP_vs_OSPP_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    
    # number of informitive transcripts
    caxcm_n_48_DE <- nrow(caxcm)
    oaxom_n_48_DE <- nrow(oaxom)
    caxcp_n_48_DE <- nrow(caxcp)
    oaxop_n_48_DE <- nrow(oaxop)
    cmxcp_n_48_DE <- nrow(cmxcp)
    omxop_n_48_DE <- nrow(omxop)
    caxna_n_48_DE <- nrow(caxna)
    oaxna_n_48_DE <- nrow(oaxna)
    caxup_n_48_DE <- nrow(caxup)
    oaxup_n_48_DE <- nrow(oaxup)
    caxoa_n_48_DE <- nrow(caxoa)
    cmxom_n_48_DE <- nrow(cmxom)
    cpxop_n_48_DE <- nrow(cpxop)
    
    # sig DE
    caxcm <- caxcm[which(caxcm$padj <= 0.05),] 
    oaxom <- oaxom[which(oaxom$padj <= 0.05),]
    caxcp <- caxcp[which(caxcp$padj <= 0.05),]
    oaxop <- oaxop[which(oaxop$padj <= 0.05),]
    cmxcp <- cmxcp[which(cmxcp$padj <= 0.05),]
    omxop <- omxop[which(omxop$padj <= 0.05),]
    caxna <- caxna[which(caxna$padj <= 0.05),]
    oaxna <- oaxna[which(oaxna$padj <= 0.05),]
    caxup <- caxup[which(caxup$padj <= 0.05),]
    oaxup <- oaxup[which(oaxup$padj <= 0.05),]
    caxoa <- caxoa[which(caxoa$padj <= 0.05),]
    cmxom <- cmxom[which(cmxom$padj <= 0.05),]
    cpxop <- cpxop[which(cpxop$padj <= 0.05),]
    
    caxcm_48_DE <- caxcm$related_accession
    oaxom_48_DE <- oaxom$related_accession
    caxcp_48_DE <- caxcp$related_accession
    oaxop_48_DE <- oaxop$related_accession
    cmxcp_48_DE <- cmxcp$related_accession
    omxop_48_DE <- omxop$related_accession
    caxna_48_DE <- caxna$related_accession
    oaxna_48_DE <- oaxna$related_accession
    caxup_48_DE <- caxup$related_accession
    oaxup_48_DE <- oaxup$related_accession
    caxoa_48_DE <- caxoa$related_accession
    cmxom_48_DE <- cmxom$related_accession
    cpxop_48_DE <- cpxop$related_accession
  }
}


# parallel incompatibilities
venn(list(caxcm_48_ME=caxcm_48_ME,cmxca_48_ME=cmxca_48_ME,caxcp_48_ME=caxcp_48_ME), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
intersect(caxcp_48_ME,intersect(caxcm_48_ME,cmxca_48_ME))
venn(list(caxcm_8_ME=caxcm_8_ME,cmxca_8_ME=cmxca_8_ME,caxcp_8_ME=caxcp_8_ME), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
intersect(caxcp_8_ME,intersect(caxcm_8_ME,cmxca_8_ME))

venn(list(oaxom_48_ME=oaxom_48_ME,omxoa_48_ME=omxoa_48_ME,oaxop_48_ME=oaxop_48_ME,opxoa_48_ME=opxoa_48_ME), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
venn(list(oaxom_8_ME=oaxom_8_ME,omxoa_8_ME=omxoa_8_ME,oaxop_8_ME=oaxop_8_ME,opxoa_8_ME=opxoa_8_ME), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
#intersect(caxcp_48_ME,intersect(caxcm_48_ME,cmxca_48_ME))
venn(list(caxcm_8_ME=caxcm_8_ME,cmxca_8_ME=cmxca_8_ME,caxcp_8_ME=caxcp_8_ME), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
#intersect(caxcp_8_ME,intersect(caxcm_8_ME,cmxca_8_ME))

# parallel DE http://nemates.org/MA/progs/overlap_stats.html
# crescent pond
venn(list(caxcm_48_DE=caxcm_48_DE,caxcp_48_DE=caxcp_48_DE), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
#intersect(caxcp_48_DE,caxcm_48_DE)
test <- matrix(c(length(intersect(caxcp_48_DE,caxcm_48_DE)),
                 length(setdiff(caxcp_48_DE,caxcm_48_DE)),
                 length(setdiff(caxcm_48_DE,caxcp_48_DE)),
                 (caxcp_n_48_DE - length(unique(c(caxcm_48_DE,caxcp_48_DE))))),ncol=2)
fisher.test(test, alternative = 'g')
venn(list(caxcm_8_DE=caxcm_8_DE,caxcp_8_DE=caxcp_8_DE), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
#intersect(caxcp_8_DE,intersect(caxcm_8_DE,cmxca_8_DE))
#intersect(intersect(caxcm_48_DE,caxcp_48_DE),intersect(caxcm_8_DE,caxcp_8_DE))
test <- matrix(c(length(intersect(caxcp_8_DE,caxcm_8_DE)),
                 length(setdiff(caxcp_8_DE,caxcm_8_DE)),
                 length(setdiff(caxcm_8_DE,caxcp_8_DE)),
                 (caxcp_n_8_DE - length(unique(c(caxcm_8_DE,caxcp_8_DE))))),ncol=2)
fisher.test(test, alternative = 'g')

# osprey
venn(list(oaxom_48_DE=oaxom_48_DE,oaxop_48_DE=oaxop_48_DE), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
#intersect(oaxom_48_DE,oaxop_48_DE)
test <- matrix(c(length(intersect(oaxop_48_DE,oaxom_48_DE)),
                 length(setdiff(oaxop_48_DE,oaxom_48_DE)),
                 length(setdiff(oaxom_48_DE,oaxop_48_DE)),
                 (oaxop_n_48_DE - length(unique(c(oaxom_48_DE,oaxop_48_DE))))),ncol=2)
fisher.test(test, alternative = 'g')
venn(list(oaxom_8_DE=oaxom_8_DE,oaxop_8_DE=oaxop_8_DE), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
#intersect(oaxom_8_DE,oaxop_8_DE)
#intersect(intersect(oaxom_48_DE,oaxop_48_DE),intersect(oaxom_8_DE,oaxop_8_DE))
test <- matrix(c(length(intersect(oaxop_8_DE,oaxom_8_DE)),
                 length(setdiff(oaxop_8_DE,oaxom_8_DE)),
                 length(setdiff(oaxom_8_DE,oaxop_8_DE)),
                 (oaxop_n_8_DE - length(unique(c(oaxom_8_DE,oaxop_8_DE))))),ncol=2)
fisher.test(test, alternative = 'g')


# scale venns
library(gridExtra)
library(eulerr)

de_am <- list(caxcm_48_DE,
               caxcm_8_DE,
               oaxom_48_DE,
               oaxom_8_DE)
de_ap <- list(caxcp_48_DE,
              caxcp_8_DE,
              oaxop_48_DE,
              oaxop_8_DE)


comp_names <- c("CRP_48hpf", 
                "CRP_8dpf", 
                "OSP_48hpf",
                "OSP_8dpf") 

bo_geness <- c()
comps <- c()
same_dirs <- c()
opp_dirs <-  c()
#same_dirs_genes <- c()
 
for (i in c(1:4)) 
#i <- 1
{
#dev.off()
  lake <- strsplit(comp_names[i],"_")[[1]][1]
  stage <- strsplit(comp_names[i],"_")[[1]][2]
  ai_genes <- de_am[[i]]
  can_genes <- de_ap[[i]]
  bo_genes <- intersect(ai_genes,can_genes)
  #bo_geness <- c(bo_geness,paste(bo_genes, sep = ";",collapse=";"))
  #bo_geness <- c(bo_geness,bo_genes)
  #comps <- c(comps,comp_names[i])
  
  ai <- length(ai_genes)
  cans <- length(can_genes)
  bo <- length(bo_genes)

  am_file <- paste(DE_genes_dir, "DE_",lake,"A_vs_",lake,"M_",stage,"_genes",".csv",sep = "")
  ap_file <- paste(DE_genes_dir, "DE_",lake,"A_vs_",lake,"P_",stage,"_genes",".csv",sep = "")
  am    <- read.csv(am_file  , header = TRUE, stringsAsFactors = FALSE)
  ap    <- read.csv(ap_file  , header = TRUE, stringsAsFactors = FALSE)
  head(am)
  am_bo <- am[am$related_accession %in% bo_genes, ]
  ap_bo <- ap[ap$related_accession %in% bo_genes, ]
  am_dn <- am_bo[which(am_bo$log2FoldChange < 0),]
  ap_dn <- ap_bo[which(ap_bo$log2FoldChange < 0),]
  am_up <- am_bo[which(am_bo$log2FoldChange > 0),]
  ap_up <- ap_bo[which(ap_bo$log2FoldChange > 0),]
  nrow(am_dn)
  nrow(am_up)
  nrow(ap_dn)
  nrow(ap_up)
  same_dir_up <- merge(am_up,ap_up, by = c("related_accession"))
  same_dir_dn <- merge(am_dn,ap_dn, by = c("related_accession"))
  opp_dir_m_up <- merge(am_up,ap_dn, by = c("related_accession"))
  opp_dir_m_dn <- merge(am_dn,ap_up, by = c("related_accession"))
  nrow(same_dir_up)
  nrow(same_dir_dn)
  nrow(opp_dir_m_up)
  nrow(opp_dir_m_dn)
  intersect(am_dn$related_accession, ap_dn$related_accession)
  
  tiff(paste("C:/Users/jmcgirr/Documents/all_2018_samples/manuscript_figs/parallel_",comp_names[i],".tiff", sep = ""), width = 4, height = 4, units = 'in', res = 1000)
  VennDiag <- euler(c("A" = ai, "B" = cans, "A&B" = bo), counts=TRUE)
  plt <- plot(VennDiag, alpha=1, lty = 1, col = "black",counts = TRUE, quantities = list(fontsize = 20, col = "white"),
              fill=c(gre,blu), labels = c("",""))#, main = comp_names[i])
  print(plt)
  dev.off()
  bars <- c(nrow(same_dir_up)+nrow(same_dir_dn),nrow(opp_dir_m_up)+nrow(opp_dir_m_dn))
  #png(paste("C:/Users/jmcgirr/Documents/all_2018_samples/manuscript_figs/parallel_bars_",comp_names[i],".png", sep = ""), width = 4, height = 4, units = 'in', res = 1000)
  
  #barplot(bars,names.arg =c("same","opposite"), col = "lightgrey", 
  #        ylim = c(0,max(bars)+(max(bars)*0.1)),cex.axis = 1.5,cex.names = 1.5)
  #dev.off()
  same_dirs <- c(same_dirs,same_dir_up$related_accession,same_dir_dn$related_accession)
  opp_dirs <- c(opp_dirs, opp_dir_m_up$related_accession,opp_dir_m_dn$related_accession)
  parallel_mp <- rbind(same_dir_up, same_dir_dn)
  parallel_mp$name.x <- NULL
  parallel_mp$name.y <- NULL
  #write.table(parallel_mp,paste("C:/Users/jmcgirr/Documents/all_2018_samples/parallel/parallel_same_dir_",comp_names[i],".csv",sep = ""), sep = ",",quote = FALSE, row.names = FALSE)
  
}

total_parallel <- length(c(unique(same_dirs),unique(opp_dirs)))
length(unique(same_dirs)) / total_parallel
length(unique(opp_dirs))
#tiff("C:/Users/jmcgirr/Documents/all_2018_samples/manuscript_figs/parallel_pie.tiff", width = 4, height = 4, units = 'in', res = 1000)
pie(c(0.9655725,0.0344275), labels = c("",""), col = c("#B2C4D2",lir))
#dev.off()


# generic read counts for parallel expression
tiff("C:/Users/jmcgirr/Documents/all_2018_samples/manuscript_figs/parallel_bars.tiff", width = 2.8, height = 2.8, units = 'in', res = 1000)
x <- c(10,30,30)
par(mfrow=c(2,2))
par(mai=c(0.2,0.4,0.2,0))
barplot(x, ylab = "", col = c(red,gre,blu), yaxt = "n")
abline(v=0.06, lwd = 1.5)
title(ylab="read counts", line=0.5, cex.lab=1)
par(mai=c(0.2,0.2,0.2,0.2))
x <- c(20,10,30)
par(mai=c(0.2,0.4,0.2,0))
barplot(x, ylab = "", col = c(red,gre,blu), yaxt = "n")
abline(v=0.06, lwd = 1.5)
par(mai=c(0.4,0.4,0,0))
x <- c(30,10,10)
par(mai=c(0.2,0.4,0.2,0))
barplot(x, ylab = "", col = c(red,gre,blu), yaxt = "n")
abline(v=0.06, lwd = 1.5)
title(ylab="read counts", line=0.5, cex.lab=1)
par(mai=c(0.4,0.2,0,0.2))
x <- c(20,30,10)
par(mai=c(0.2,0.4,0.2,0))
barplot(x, ylab = "", col = c(red,gre,blu), yaxt = "n")
abline(v=0.06, lwd = 1.5)
dev.off()





pmp_c_48 <- read.csv("C:/Users/jmcgirr/Documents/all_2018_samples/parallel/parallel_same_dir_mp_crp_48.csv", header = TRUE, stringsAsFactors = FALSE)
pmp_c_8 <-  read.csv("C:/Users/jmcgirr/Documents/all_2018_samples/parallel/parallel_same_dir_mp_crp_8.csv", header = TRUE, stringsAsFactors = FALSE)
pmp_o_48 <- read.csv("C:/Users/jmcgirr/Documents/all_2018_samples/parallel/parallel_same_dir_mp_osp_48.csv", header = TRUE, stringsAsFactors = FALSE)
pmp_o_8 <-  read.csv("C:/Users/jmcgirr/Documents/all_2018_samples/parallel/parallel_same_dir_mp_osp_8.csv", header = TRUE, stringsAsFactors = FALSE)
pmp_c_48 <- pmp_c_48$related_accession
pmp_c_8 <-  pmp_c_8$related_accession
pmp_o_48 <- pmp_o_48$related_accession
pmp_o_8 <-  pmp_o_8$related_accession

all_pmp <- c(pmp_c_48,
             pmp_c_8, 
             pmp_o_48,
             pmp_o_8, 
             pmp_c_48,
             pmp_c_8, 
             pmp_o_48,
             pmp_o_8 )

oaxom_8_ai <- intersect(intersect(oaxom_8_DE,oaxom_8_ME),omxoa_8_ME)
caxcp_8_ai <- intersect(caxcp_8_DE,caxcp_8_ME)
oaxop_48_ai <- intersect(intersect(oaxop_48_DE,oaxop_48_ME),opxoa_48_ME)
oaxop_8_ai  <- intersect(intersect(oaxop_8_DE,oaxop_8_ME),opxoa_8_ME)
cmxcp_48_ai <- intersect(cmxcp_48_DE,cmxcp_48_ME)
cmxcp_8_ai <-intersect(cmxcp_8_DE,cmxcp_8_ME)
omxop_48_ai <- intersect(omxop_48_DE,opxom_48_ME)
omxop_8_ai <- intersect(omxop_8_DE,opxom_8_ME)

# AI and parallel DE MP
intersect(caxcp_8_ai,pmp_c_8)
intersect(oaxop_48_ai,pmp_o_48)
intersect(oaxop_8_ai,pmp_o_8)

# parallel AI
intersect(caxcp_48_ai,caxcm_48_ai)
intersect(caxcp_8_ai,caxcm_8_ai)
intersect(oaxop_48_ai,oaxom_48_ai)
intersect(oaxop_8_ai,oaxom_8_ai)

# parallel DE MP and ME MP
intersect(cmxcp_48_ME,pmp_c_48)
intersect(cmxcp_8_ME,pmp_c_8)
intersect(opxom_48_ME,pmp_o_48)
intersect(opxom_8_ME,pmp_o_8)

# parallel DE MP and super ME MP
intersect(cmxcp_48_SME,pmp_c_48)
intersect(cmxcp_8_SME,pmp_c_8)
intersect(opxom_48_SME,pmp_o_48)
intersect(opxom_8_SME,pmp_o_8)


all_ai <- c(oaxom_8_ai, 
            caxcp_8_ai, 
            oaxop_48_ai,
            oaxop_8_ai, 
            cmxcp_48_ai,
            cmxcp_8_ai, 
            omxop_48_ai,
            omxop_8_ai )

intersect(all_pmp,all_ai)



