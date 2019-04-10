######### load all ai, parallel, DE, ME, super ME ######
library(venn)
xm_to_gene <- function(xms){
  
  xm_transcripts <- data.frame(related_accession = xms)
  xm_transcripts <- merge(xm_transcripts,final_features, by = c("related_accession"))
  xm_transcripts<- merge(xm_transcripts, blast_key, by = c("product_accession"))
  print("cyprinodon gene annotations")
  print(xm_transcripts$symbol)
  print("danio ortholog annotations")
  print(xm_transcripts$zeb_gene_symbol_one_way)
  print("danio orthologs annotated for skeletal effects")
  xm_cranial <- merge(xm_transcripts, cranial_genes, by = c("zeb_gene_symbol_one_way"))
  print(xm_cranial$zeb_gene_symbol_one_way)
  print("cyprinodon genes annotated for skeletal effects")
  xm_cranial_cyp <- merge(xm_transcripts, cranial_genes, by = c("symbol"))
  print(xm_cranial_cyp$symbol)
}

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
#load parallel mp DE
pmp_c_48 <- read.csv("C:/Users/jmcgirr/Documents/all_2018_samples/parallel/parallel_same_dir_CRP_48hpf.csv", header = TRUE, stringsAsFactors = FALSE)
pmp_c_8 <-  read.csv("C:/Users/jmcgirr/Documents/all_2018_samples/parallel/parallel_same_dir_CRP_8dpf.csv", header = TRUE, stringsAsFactors = FALSE)
pmp_o_48 <- read.csv("C:/Users/jmcgirr/Documents/all_2018_samples/parallel/parallel_same_dir_OSP_48hpf.csv", header = TRUE, stringsAsFactors = FALSE)
pmp_o_8 <-  read.csv("C:/Users/jmcgirr/Documents/all_2018_samples/parallel/parallel_same_dir_OSP_8dpf.csv", header = TRUE, stringsAsFactors = FALSE)
pmp_c_48 <- pmp_c_48$related_accession
pmp_c_8 <-  pmp_c_8$related_accession
pmp_o_48 <- pmp_o_48$related_accession
pmp_o_8 <-  pmp_o_8$related_accession
#load super misexpressed
{
  DE_genes_dir <- "C:/Users/jmcgirr/Documents/all_2018_samples/super_misexpression/conditions/"
  stage <- "8dpf"
  # load misexpressed
  {
    # am crosses
    caxcm <- read.csv(paste(DE_genes_dir,"DE_CAxCM_vs_all_pure_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    cmxca <- read.csv(paste(DE_genes_dir,"DE_CMxCA_vs_all_pure_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    oaxom <- read.csv(paste(DE_genes_dir,"DE_OAxOM_vs_all_pure_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    omxoa <- read.csv(paste(DE_genes_dir,"DE_OMxOA_vs_all_pure_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    # ap crosses
    caxcp <- read.csv(paste(DE_genes_dir,"DE_CAxCP_vs_all_pure_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    oaxop <- read.csv(paste(DE_genes_dir,"DE_OAxOP_vs_all_pure_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    opxoa <- read.csv(paste(DE_genes_dir,"DE_OPxOA_vs_all_pure_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    # specialist crosses
    cmxcp <- read.csv(paste(DE_genes_dir,"DE_CMxCP_vs_all_pure_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    opxom <- read.csv(paste(DE_genes_dir,"DE_OPxOM_vs_all_pure_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    
    # sig super misexpressed
    caxcm <- caxcm[which(caxcm$padj <= 0.05),] 
    cmxca <- cmxca[which(cmxca$padj <= 0.05),] 
    oaxom <- oaxom[which(oaxom$padj <= 0.05),] 
    omxoa <- omxoa[which(omxoa$padj <= 0.05),] 
    caxcp <- caxcp[which(caxcp$padj <= 0.05),] 
    oaxop <- oaxop[which(oaxop$padj <= 0.05),] 
    opxoa <- opxoa[which(opxoa$padj <= 0.05),] 
    cmxcp <- cmxcp[which(cmxcp$padj <= 0.05),] 
    opxom <- opxom[which(opxom$padj <= 0.05),] 
    
    
    caxcm_8_SME <- caxcm$related_accession
    cmxca_8_SME <- cmxca$related_accession
    oaxom_8_SME <- oaxom$related_accession
    omxoa_8_SME <- omxoa$related_accession
    caxcp_8_SME <- caxcp$related_accession
    oaxop_8_SME <- oaxop$related_accession
    opxoa_8_SME <- opxoa$related_accession
    cmxcp_8_SME <- cmxcp$related_accession
    opxom_8_SME <- opxom$related_accession
    
  }
  stage <- "48hpf"
  # load misexpressed
  {
    # am crosses
    caxcm <- read.csv(paste(DE_genes_dir,"DE_CAxCM_vs_all_pure_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    cmxca <- read.csv(paste(DE_genes_dir,"DE_CMxCA_vs_all_pure_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    oaxom <- read.csv(paste(DE_genes_dir,"DE_OAxOM_vs_all_pure_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    omxoa <- read.csv(paste(DE_genes_dir,"DE_OMxOA_vs_all_pure_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    # ap crosses
    caxcp <- read.csv(paste(DE_genes_dir,"DE_CAxCP_vs_all_pure_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    oaxop <- read.csv(paste(DE_genes_dir,"DE_OAxOP_vs_all_pure_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    opxoa <- read.csv(paste(DE_genes_dir,"DE_OPxOA_vs_all_pure_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    # specialist crosses
    cmxcp <- read.csv(paste(DE_genes_dir,"DE_CMxCP_vs_all_pure_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    opxom <- read.csv(paste(DE_genes_dir,"DE_OPxOM_vs_all_pure_",stage,"_genes.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    
    # sig super misexpressed
    caxcm <- caxcm[which(caxcm$padj <= 0.05),] 
    cmxca <- cmxca[which(cmxca$padj <= 0.05),] 
    oaxom <- oaxom[which(oaxom$padj <= 0.05),] 
    omxoa <- omxoa[which(omxoa$padj <= 0.05),] 
    caxcp <- caxcp[which(caxcp$padj <= 0.05),] 
    oaxop <- oaxop[which(oaxop$padj <= 0.05),] 
    opxoa <- opxoa[which(opxoa$padj <= 0.05),] 
    cmxcp <- cmxcp[which(cmxcp$padj <= 0.05),] 
    opxom <- opxom[which(opxom$padj <= 0.05),] 
    
    
    caxcm_48_SME <- caxcm$related_accession
    cmxca_48_SME <- cmxca$related_accession
    oaxom_48_SME <- oaxom$related_accession
    omxoa_48_SME <- omxoa$related_accession
    caxcp_48_SME <- caxcp$related_accession
    oaxop_48_SME <- oaxop$related_accession
    opxoa_48_SME <- opxoa$related_accession
    cmxcp_48_SME <- cmxcp$related_accession
    opxom_48_SME <- opxom$related_accession
    
  }
}
#load ai
oaxom_8_ai <- intersect(intersect(oaxom_8_DE,oaxom_8_ME),omxoa_8_ME)
caxcp_8_ai <- intersect(caxcp_8_DE,caxcp_8_ME)
oaxop_48_ai <- intersect(intersect(oaxop_48_DE,oaxop_48_ME),opxoa_48_ME)
oaxop_8_ai  <- intersect(intersect(oaxop_8_DE,oaxop_8_ME),opxoa_8_ME)
cmxcp_48_ai <- intersect(cmxcp_48_DE,cmxcp_48_ME)
cmxcp_8_ai <-intersect(cmxcp_8_DE,cmxcp_8_ME)
omxop_48_ai <- intersect(omxop_48_DE,opxom_48_ME)
omxop_8_ai <- intersect(omxop_8_DE,opxom_8_ME)
#combine ME for reciprocal crosses
caxcm_48_ME <- intersect(caxcm_48_ME,cmxca_48_ME)
caxcm_8_ME <- intersect(caxcm_8_ME,cmxca_8_ME)
oaxom_48_ME <- intersect(oaxom_48_ME,omxoa_48_ME)
oaxom_8_ME <- intersect(oaxom_8_ME,omxoa_8_ME)
oaxop_48_ME <- intersect(oaxop_48_ME,opxoa_48_ME)
oaxop_8_ME <- intersect(oaxop_8_ME,opxoa_8_ME)
#load genes associated with jaw size
require(data.table)
library(reshape2)

final_features <- read.table("D:/Martin Lab/rna_2018/all_2018_samples/features_gff.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t", quote = "")
blast_key <- read.table("D:/Martin Lab/rna_2018/all_2018_samples/blast/cyprinodon_to_danio_one_way_best_hit_symbols.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
mrna <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/ase_data/mrna.saf", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
mrna <- cbind(mrna, colsplit(mrna$GeneID, ";", c("related_accession", "gene_name")))
mrna_table <- mrna
colnames(mrna_table) <- c("GeneID","CHROM","START","END","Strand","related_accession", "gene_name")
mrna_table$GeneID <- NULL
mrna_table$START <- mrna_table$START -10000
mrna_table$END <- mrna_table$END +10000
setDT(mrna_table)
gem <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/GEMMA/candidates_up_jaw_99percentile.bed", header = FALSE, stringsAsFactors = FALSE)
colnames(gem) <- c("CHROM","START","END","PIP")
setDT(gem)
setkey(mrna_table)
c <- foverlaps(gem, mrna_table,by.x = c("CHROM", "START", "END"),by.y = c("CHROM", "START", "END"),type="any", nomatch=0L)
c <- as.data.frame(c)
c <- c[order(c$PIP, decreasing = TRUE),]
rownames(c) <- NULL
gem_cans <- unique(c$related_accession)

#load genes with fixed snps
caxcm_fixed <- read.csv("C:/Users/jmcgirr/Documents/all_2018_samples/interesting_genes_tables/caxcm_fixed_snp_genes.csv", header = TRUE, row.names = NULL,stringsAsFactors = FALSE)
caxcp_fixed <- read.csv("C:/Users/jmcgirr/Documents/all_2018_samples/interesting_genes_tables/caxcp_fixed_snp_genes.csv", header = TRUE, row.names = NULL,stringsAsFactors = FALSE)
oaxom_fixed <- read.csv("C:/Users/jmcgirr/Documents/all_2018_samples/interesting_genes_tables/oaxom_fixed_snp_genes.csv", header = TRUE, row.names = NULL,stringsAsFactors = FALSE)
oaxop_fixed <- read.csv("C:/Users/jmcgirr/Documents/all_2018_samples/interesting_genes_tables/oaxop_fixed_snp_genes.csv", header = TRUE, row.names = NULL,stringsAsFactors = FALSE)
cmxcp_fixed <- read.csv("C:/Users/jmcgirr/Documents/all_2018_samples/interesting_genes_tables/cmxcp_fixed_snp_genes.csv", header = TRUE, row.names = NULL,stringsAsFactors = FALSE)
omxop_fixed <- read.csv("C:/Users/jmcgirr/Documents/all_2018_samples/interesting_genes_tables/omxop_fixed_snp_genes.csv", header = TRUE, row.names = NULL,stringsAsFactors = FALSE)
caxcm_fixed <- unique(caxcm_fixed$related_accession)  
caxcp_fixed <- unique(caxcp_fixed$related_accession)
oaxom_fixed <- unique(oaxom_fixed$related_accession)
oaxop_fixed <- unique(oaxop_fixed$related_accession)
cmxcp_fixed <- unique(cmxcp_fixed$related_accession)
omxop_fixed <- unique(omxop_fixed$related_accession)

#load genes showing selection sweed
ca_sweed <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/sweed/correct_seq/ca_pop_bottle_58_sweeps_90_genes.txt", header = TRUE, row.names = NULL, stringsAsFactors = FALSE)
cm_sweed <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/sweed/correct_seq/cm_pop_bottle_58_sweeps_90_genes.txt", header = TRUE, row.names = NULL, stringsAsFactors = FALSE)
cp_sweed <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/sweed/correct_seq/cp_pop_bottle_58_sweeps_90_genes.txt", header = TRUE, row.names = NULL, stringsAsFactors = FALSE)
oa_sweed <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/sweed/correct_seq/oa_pop_bottle_58_sweeps_90_genes.txt", header = TRUE, row.names = NULL, stringsAsFactors = FALSE)
om_sweed <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/sweed/correct_seq/om_pop_bottle_58_sweeps_90_genes.txt", header = TRUE, row.names = NULL, stringsAsFactors = FALSE)
op_sweed <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/sweed/correct_seq/op_pop_bottle_58_sweeps_90_genes.txt", header = TRUE, row.names = NULL, stringsAsFactors = FALSE)
ca_sweed <- unique(ca_sweed$related_accession)
cm_sweed <- unique(cm_sweed$related_accession)
cp_sweed <- unique(cp_sweed$related_accession)
oa_sweed <- unique(oa_sweed$related_accession)
om_sweed <- unique(om_sweed$related_accession)
op_sweed <- unique(op_sweed$related_accession)
caxcm_sweed <- unique(c(ca_sweed,cm_sweed))
caxcp_sweed<-  unique(c(ca_sweed,cp_sweed))
cmxcp_sweed<-  unique(c(cm_sweed,cp_sweed))
oaxom_sweed<-  unique(c(oa_sweed,om_sweed))
oaxop_sweed<-  unique(c(oa_sweed,op_sweed))
omxop_sweed<-  unique(c(om_sweed,op_sweed))
#load genes showing selection taj
ca_taj <- read.csv("C:/Users/jmcgirr/Documents/all_2018_samples/interesting_genes_tables/ca_taj_90.csv", header = TRUE,row.names = NULL, stringsAsFactors = FALSE)
cm_taj <- read.csv("C:/Users/jmcgirr/Documents/all_2018_samples/interesting_genes_tables/cm_taj_90.csv", header = TRUE,row.names = NULL, stringsAsFactors = FALSE)
cp_taj <- read.csv("C:/Users/jmcgirr/Documents/all_2018_samples/interesting_genes_tables/cp_taj_90.csv", header = TRUE,row.names = NULL, stringsAsFactors = FALSE)
oa_taj <- read.csv("C:/Users/jmcgirr/Documents/all_2018_samples/interesting_genes_tables/oa_taj_90.csv", header = TRUE,row.names = NULL, stringsAsFactors = FALSE)
om_taj <- read.csv("C:/Users/jmcgirr/Documents/all_2018_samples/interesting_genes_tables/om_taj_90.csv", header = TRUE,row.names = NULL, stringsAsFactors = FALSE)
op_taj <- read.csv("C:/Users/jmcgirr/Documents/all_2018_samples/interesting_genes_tables/op_taj_90.csv", header = TRUE,row.names = NULL, stringsAsFactors = FALSE)
ca_taj <- unique(ca_taj$related_accession)
cm_taj <- unique(cm_taj$related_accession)
cp_taj <- unique(cp_taj$related_accession)
oa_taj <- unique(oa_taj$related_accession)
om_taj <- unique(om_taj$related_accession)
op_taj <- unique(op_taj$related_accession)
caxcm_taj <- unique(c(ca_taj,cm_taj))
caxcp_taj<-  unique(c(ca_taj,cp_taj))
cmxcp_taj<-  unique(c(cm_taj,cp_taj))
oaxom_taj<-  unique(c(oa_taj,om_taj))
oaxop_taj<-  unique(c(oa_taj,op_taj))
omxop_taj<-  unique(c(om_taj,op_taj))
#load dxy
caxcm_dxy <- read.csv("C:/Users/jmcgirr/Documents/all_2018_samples/interesting_genes_tables/caxcm_corr_dxy_90.csv", header = TRUE, row.names = NULL,stringsAsFactors = FALSE)
caxcp_dxy <- read.csv("C:/Users/jmcgirr/Documents/all_2018_samples/interesting_genes_tables/caxcp_corr_dxy_90.csv", header = TRUE, row.names = NULL,stringsAsFactors = FALSE)
oaxom_dxy <- read.csv("C:/Users/jmcgirr/Documents/all_2018_samples/interesting_genes_tables/oaxom_corr_dxy_90.csv", header = TRUE, row.names = NULL,stringsAsFactors = FALSE)
oaxop_dxy <- read.csv("C:/Users/jmcgirr/Documents/all_2018_samples/interesting_genes_tables/oaxop_corr_dxy_90.csv", header = TRUE, row.names = NULL,stringsAsFactors = FALSE)
cmxcp_dxy <- read.csv("C:/Users/jmcgirr/Documents/all_2018_samples/interesting_genes_tables/cmxcp_corr_dxy_90.csv", header = TRUE, row.names = NULL,stringsAsFactors = FALSE)
omxop_dxy <- read.csv("C:/Users/jmcgirr/Documents/all_2018_samples/interesting_genes_tables/omxop_corr_dxy_90.csv", header = TRUE, row.names = NULL,stringsAsFactors = FALSE)
caxcm_dxy <- unique(caxcm_dxy$related_accession)  
caxcp_dxy <- unique(caxcp_dxy$related_accession)
oaxom_dxy <- unique(oaxom_dxy$related_accession)
oaxop_dxy <- unique(oaxop_dxy$related_accession)
cmxcp_dxy <- unique(cmxcp_dxy$related_accession)
omxop_dxy <- unique(omxop_dxy$related_accession)






# number of genes filtered through pipline
# ai 
all_ai <- unique(c(oaxom_8_ai, 
                caxcp_8_ai ,
                oaxop_48_ai,
                oaxop_8_ai ,
                cmxcp_48_ai,
                cmxcp_8_ai ,
                omxop_48_ai,
                omxop_8_ai ))
# ai_fst_dxy
ai_fst_dxy <- unique(c(intersect(oaxom_8_ai ,intersect(oaxom_dxy,oaxom_fixed)),
                intersect(caxcp_8_ai ,intersect(caxcp_dxy,caxcp_fixed)),
                intersect(oaxop_48_ai,intersect(oaxop_dxy,oaxop_fixed)),
                intersect(oaxop_8_ai ,intersect(oaxop_dxy,oaxop_fixed)),
                intersect(cmxcp_48_ai,intersect(cmxcp_dxy,cmxcp_fixed)),
                intersect(cmxcp_8_ai ,intersect(cmxcp_dxy,cmxcp_fixed)),
                intersect(omxop_48_ai,intersect(omxop_dxy,omxop_fixed)),
                intersect(omxop_8_ai ,intersect(omxop_dxy,omxop_fixed))))
ai_fst_dxy_sweed <- unique(c(intersect(intersect(oaxom_8_ai ,intersect(oaxom_dxy,oaxom_fixed)),om_taj),
                             intersect(intersect(caxcp_8_ai ,intersect(caxcp_dxy,caxcp_fixed)),cp_taj),
                             intersect(intersect(oaxop_48_ai,intersect(oaxop_dxy,oaxop_fixed)),op_taj),
                             intersect(intersect(oaxop_8_ai ,intersect(oaxop_dxy,oaxop_fixed)),op_taj),
                             intersect(intersect(cmxcp_48_ai,intersect(cmxcp_dxy,cmxcp_fixed)),cmxcp_taj),
                             intersect(intersect(cmxcp_8_ai ,intersect(cmxcp_dxy,cmxcp_fixed)),cmxcp_taj),
                             intersect(intersect(omxop_48_ai,intersect(omxop_dxy,omxop_fixed)),omxop_taj),
                             intersect(intersect(omxop_8_ai ,intersect(omxop_dxy,omxop_fixed)),omxop_taj)))
xm_to_gene(intersect(oaxom_8_ai ,intersect(oaxom_dxy,oaxom_fixed)))
xm_to_gene(intersect(caxcp_8_ai ,intersect(caxcp_dxy,caxcp_fixed)))
xm_to_gene(intersect(oaxop_48_ai,intersect(oaxop_dxy,oaxop_fixed)))
xm_to_gene(intersect(oaxop_8_ai ,intersect(oaxop_dxy,oaxop_fixed)))
xm_to_gene(intersect(cmxcp_48_ai,intersect(cmxcp_dxy,cmxcp_fixed)))
xm_to_gene(intersect(cmxcp_8_ai ,intersect(cmxcp_dxy,cmxcp_fixed)))
xm_to_gene(intersect(omxop_48_ai,intersect(omxop_dxy,omxop_fixed)))
xm_to_gene(intersect(omxop_8_ai ,intersect(omxop_dxy,omxop_fixed)))

# ai_fst_gemma
ai_fst_dxy_gemma <- unique(c(intersect(gem_cans,intersect(oaxom_8_ai ,intersect(oaxom_dxy,oaxom_fixed))),
                intersect(gem_cans,intersect(caxcp_8_ai ,intersect(caxcp_dxy,caxcp_fixed))),
                intersect(gem_cans,intersect(oaxop_48_ai,intersect(oaxop_dxy,oaxop_fixed))),
                intersect(gem_cans,intersect(oaxop_8_ai ,intersect(oaxop_dxy,oaxop_fixed))),
                intersect(gem_cans,intersect(cmxcp_48_ai,intersect(cmxcp_dxy,cmxcp_fixed))),
                intersect(gem_cans,intersect(cmxcp_8_ai ,intersect(cmxcp_dxy,cmxcp_fixed))),
                intersect(gem_cans,intersect(omxop_48_ai,intersect(omxop_dxy,omxop_fixed))),
                intersect(gem_cans,intersect(omxop_8_ai ,intersect(omxop_dxy,omxop_fixed)))))

intersect(gem_cans,intersect(oaxom_8_ai ,intersect(oaxom_dxy,oaxom_fixed)))
intersect(gem_cans,intersect(caxcp_8_ai ,intersect(caxcp_dxy,caxcp_fixed)))
intersect(gem_cans,intersect(oaxop_48_ai,intersect(oaxop_dxy,oaxop_fixed)))
intersect(gem_cans,intersect(oaxop_8_ai ,intersect(oaxop_dxy,oaxop_fixed)))
intersect(gem_cans,intersect(cmxcp_48_ai,intersect(cmxcp_dxy,cmxcp_fixed)))
intersect(gem_cans,intersect(cmxcp_8_ai ,intersect(cmxcp_dxy,cmxcp_fixed)))
intersect(gem_cans,intersect(omxop_48_ai,intersect(omxop_dxy,omxop_fixed)))
intersect(gem_cans,intersect(omxop_8_ai ,intersect(omxop_dxy,omxop_fixed)))

# pmp ME
pmp_me <- unique(c(intersect(pmp_c_8 ,cmxcp_8_ME),
                intersect(pmp_o_8 ,opxom_8_ME),
                intersect(pmp_c_48,cmxcp_48_ME),
                intersect(pmp_o_48,opxom_48_ME)))

intersect(pmp_c_8 ,cmxcp_8_ME)
intersect(pmp_o_8 ,opxom_8_ME)
intersect(pmp_c_48,cmxcp_48_ME)
intersect(pmp_o_48,opxom_48_ME)
intersect(intersect(pmp_c_8 ,cmxcp_8_ME),intersect(pmp_o_8 ,opxom_8_ME))
intersect(intersect(pmp_c_48 ,cmxcp_48_ME),intersect(pmp_o_48 ,opxom_48_ME))

# pmp ME fst_dxy
pmp_me_fst_dxy <- unique(c(intersect(intersect(pmp_c_8 ,cmxcp_8_ME) ,intersect(cmxcp_dxy,cmxcp_fixed)), 
                           intersect(intersect(pmp_o_8 ,opxom_8_ME) ,intersect(omxop_dxy,omxop_fixed)), 
                           intersect(intersect(pmp_c_48,cmxcp_48_ME),intersect(cmxcp_dxy,cmxcp_fixed)), 
                           intersect(intersect(pmp_o_48,opxom_48_ME),intersect(omxop_dxy,omxop_fixed))))
pmp_me_fst_dxy_sweed <- unique(c(intersect(intersect(intersect(pmp_c_8 ,cmxcp_8_ME) ,intersect(cmxcp_dxy,cmxcp_fixed)),cmxcp_taj), 
                                 intersect(intersect(intersect(pmp_o_8 ,opxom_8_ME) ,intersect(omxop_dxy,omxop_fixed)),omxop_taj), 
                                 intersect(intersect(intersect(pmp_c_48,cmxcp_48_ME),intersect(cmxcp_dxy,cmxcp_fixed)),cmxcp_taj), 
                                 intersect(intersect(intersect(pmp_o_48,opxom_48_ME),intersect(omxop_dxy,omxop_fixed)),omxop_taj)))
intersect(intersect(pmp_c_8 ,cmxcp_8_ME) ,intersect(cmxcp_dxy,cmxcp_fixed)) 
intersect(intersect(pmp_o_8 ,opxom_8_ME) ,intersect(omxop_dxy,omxop_fixed)) 
intersect(intersect(pmp_c_48,cmxcp_48_ME),intersect(cmxcp_dxy,cmxcp_fixed)) 
intersect(intersect(pmp_o_48,opxom_48_ME),intersect(omxop_dxy,omxop_fixed)) 
# pmp ME fst_dxy_gemma
pmp_fst_dxy_gemma <- unique(c(intersect(gem_cans,intersect(intersect(pmp_c_8 ,cmxcp_8_ME) ,intersect(cmxcp_dxy,cmxcp_fixed))), 
                              intersect(gem_cans,intersect(intersect(pmp_o_8 ,opxom_8_ME) ,intersect(omxop_dxy,omxop_fixed))), 
                              intersect(gem_cans,intersect(intersect(pmp_c_48,cmxcp_48_ME),intersect(cmxcp_dxy,cmxcp_fixed))), 
                              intersect(gem_cans,intersect(intersect(pmp_o_48,opxom_48_ME),intersect(omxop_dxy,omxop_fixed)))))

intersect(gem_cans,intersect(intersect(pmp_c_8 ,cmxcp_8_ME) ,intersect(cmxcp_dxy,cmxcp_fixed))) 
intersect(gem_cans,intersect(intersect(pmp_o_8 ,opxom_8_ME) ,intersect(omxop_dxy,omxop_fixed))) 
intersect(gem_cans,intersect(intersect(pmp_c_48,cmxcp_48_ME),intersect(cmxcp_dxy,cmxcp_fixed))) 
intersect(gem_cans,intersect(intersect(pmp_o_48,opxom_48_ME),intersect(omxop_dxy,omxop_fixed))) 

out_table <- mrna[mrna$related_accession %in% pmp_me_fst_dxy, ]
#write.table(out_table,"C:/Users/jmcgirr/Documents/all_2018_samples/GO/pmp_me_fst_dxy_genes.txt", row.names = FALSE, quote = FALSE, sep = "\t")


"XM_015398227.1" %in% intersect(oaxom_8_ai ,intersect(oaxom_dxy,oaxom_fixed))
"XM_015398227.1" %in% intersect(caxcp_8_ai ,intersect(caxcp_dxy,caxcp_fixed))
"XM_015398227.1" %in% intersect(oaxop_48_ai,intersect(oaxop_dxy,oaxop_fixed))
"XM_015398227.1" %in% intersect(oaxop_8_ai ,intersect(oaxop_dxy,oaxop_fixed))
"XM_015398227.1" %in% intersect(cmxcp_48_ai,intersect(cmxcp_dxy,cmxcp_fixed))
"XM_015398227.1" %in% intersect(cmxcp_8_ai ,intersect(cmxcp_dxy,cmxcp_fixed))
"XM_015398227.1" %in% intersect(omxop_48_ai,intersect(omxop_dxy,omxop_fixed))
"XM_015398227.1" %in% intersect(omxop_8_ai ,intersect(omxop_dxy,omxop_fixed))


"XM_015369861.1" %in% intersect(intersect(pmp_c_8 ,cmxcp_8_ME) ,intersect(cmxcp_dxy,cmxcp_fixed)) 
"XM_015369861.1" %in% intersect(intersect(pmp_o_8 ,opxom_8_ME) ,intersect(omxop_dxy,omxop_fixed)) 
"XM_015369861.1" %in% intersect(intersect(pmp_c_48,cmxcp_48_ME),intersect(cmxcp_dxy,cmxcp_fixed)) 
"XM_015369861.1" %in% intersect(intersect(pmp_o_48,opxom_48_ME),intersect(omxop_dxy,omxop_fixed))

"XM_015374108.1" %in% intersect(pmp_c_8 ,cmxcp_8_ME) 
"XM_015374108.1" %in% intersect(pmp_o_8 ,opxom_8_ME) 
"XM_015374108.1" %in% intersect(pmp_c_48,cmxcp_48_ME)
"XM_015374108.1" %in% intersect(pmp_o_48,opxom_48_ME)




par(mfrow=c(1,1))
# AI and selection and gemma and fixed
venn(list(cmxcp_8_ai=cmxcp_8_ai,cmxcp_fixed=cmxcp_fixed,cmxcp_sweed=cmxcp_sweed, gem_cans= gem_cans,cmxcp_taj=cmxcp_taj), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
venn(list(caxcp_8_ai=caxcp_8_ai,caxcp_fixed=caxcp_fixed,caxcp_sweed=caxcp_sweed, gem_cans= gem_cans,caxcp_taj=caxcp_taj), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
venn(list(omxop_8_ai=omxop_8_ai,omxop_fixed=omxop_fixed,omxop_sweed=omxop_sweed, gem_cans= gem_cans,omxop_taj=omxop_taj), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
venn(list(oaxop_8_ai=oaxop_8_ai,oaxop_fixed=oaxop_fixed,oaxop_sweed=oaxop_sweed, gem_cans= gem_cans,oaxop_taj=oaxop_taj), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)

# AI and fixed and gemma
venn(list(cmxcp_8_ai=cmxcp_8_ai,cmxcp_fixed=cmxcp_fixed, gem_cans= gem_cans,cmxcp_dxy=cmxcp_dxy), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
venn(list(caxcp_8_ai=caxcp_8_ai,caxcp_fixed=caxcp_fixed, gem_cans= gem_cans,caxcp_dxy=caxcp_dxy), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
venn(list(omxop_8_ai=omxop_8_ai,omxop_fixed=omxop_fixed, gem_cans= gem_cans), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
venn(list(oaxop_8_ai=oaxop_8_ai,oaxop_fixed=oaxop_fixed, gem_cans= gem_cans,oaxop_dxy=oaxop_dxy), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
venn(list(oaxom_8_ai=oaxom_8_ai,oaxom_fixed=oaxom_fixed, gem_cans= gem_cans,oaxop_dxy=oaxop_dxy), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
venn(list(cmxcp_48_ai=cmxcp_48_ai,cmxcp_fixed=cmxcp_fixed, gem_cans= gem_cans,oaxop_dxy=oaxop_dxy), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
venn(list(omxop_48_ai=omxop_48_ai,omxop_fixed=omxop_fixed, gem_cans= gem_cans,omxop_dxy=omxop_dxy), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
venn(list(oaxop_48_ai=oaxop_48_ai,oaxop_fixed=oaxop_fixed, gem_cans= gem_cans,oaxop_dxy=oaxop_dxy), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
length(unique(c(cmxcp_8_ai,caxcp_8_ai,omxop_8_ai,oaxom_8_ai,cmxcp_48_ai,omxop_48_ai,oaxop_48_ai)))

intersect(intersect(oaxom_8_ai ,gem_cans),intersect(oaxom_dxy,oaxom_fixed))
intersect(intersect(caxcp_8_ai ,gem_cans),intersect(caxcp_dxy,caxcp_fixed))
intersect(intersect(oaxop_48_ai,gem_cans),intersect(oaxop_dxy,oaxop_fixed))
intersect(intersect(oaxop_8_ai ,gem_cans),intersect(oaxop_dxy,oaxop_fixed))
intersect(intersect(cmxcp_48_ai,gem_cans),intersect(cmxcp_dxy,cmxcp_fixed))
intersect(intersect(cmxcp_8_ai ,gem_cans),intersect(cmxcp_dxy,cmxcp_fixed))
intersect(intersect(omxop_48_ai,gem_cans),intersect(omxop_dxy,omxop_fixed))
intersect(intersect(omxop_8_ai ,gem_cans),intersect(omxop_dxy,omxop_fixed))

length(intersect(oaxom_8_ai ,intersect(oaxom_dxy,oaxom_fixed)))
length(intersect(caxcp_8_ai ,intersect(caxcp_dxy,caxcp_fixed)))
length(intersect(oaxop_48_ai,intersect(oaxop_dxy,oaxop_fixed)))
length(intersect(oaxop_8_ai ,intersect(oaxop_dxy,oaxop_fixed)))
length(intersect(cmxcp_48_ai,intersect(cmxcp_dxy,cmxcp_fixed)))
length(intersect(cmxcp_8_ai ,intersect(cmxcp_dxy,cmxcp_fixed)))
length(intersect(omxop_48_ai,intersect(omxop_dxy,omxop_fixed)))
length(intersect(omxop_8_ai ,intersect(omxop_dxy,omxop_fixed)))
xm_to_gene(intersect(intersect(oaxom_8_ai ,gem_cans),intersect(oaxom_dxy,oaxom_fixed)))
xm_to_gene(intersect(intersect(caxcp_8_ai ,gem_cans),intersect(caxcp_dxy,caxcp_fixed)))
xm_to_gene(intersect(intersect(oaxop_48_ai,gem_cans),intersect(oaxop_dxy,oaxop_fixed)))
xm_to_gene(intersect(intersect(oaxop_8_ai ,gem_cans),intersect(oaxop_dxy,oaxop_fixed)))
xm_to_gene(intersect(intersect(cmxcp_48_ai,gem_cans),intersect(cmxcp_dxy,cmxcp_fixed)))
xm_to_gene(intersect(intersect(cmxcp_8_ai ,gem_cans),intersect(cmxcp_dxy,cmxcp_fixed)))
xm_to_gene(intersect(intersect(omxop_48_ai,gem_cans),intersect(omxop_dxy,omxop_fixed)))
xm_to_gene(intersect(intersect(omxop_8_ai ,gem_cans),intersect(omxop_dxy,omxop_fixed)))
all <- c("XM_015371309.1","XM_015372614.1 ","XM_015372485.1 ","XM_015369602.1 ","XM_015379462.1","XM_015394812.1","XM_015386614.1 ","XM_015392887.1 ","XM_015402993.1 ","XM_015394812.1 ","XM_015384645.1 ","XM_015384780.1","XM_015377809.1 ","XM_015397247.1 ","XM_015399685.1")
all <- trimws(all, which = c("both", "left", "right"))
xm_to_gene(all)
# parallel DE MP and fixed and gemma
venn(list(pmp_c_8=pmp_c_8,cmxcp_fixed=cmxcp_fixed, gem_cans= gem_cans, cmxcp_dxy=cmxcp_dxy), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
venn(list(pmp_o_8=pmp_o_8,omxop_fixed=omxop_fixed, gem_cans= gem_cans,omxop_dxy=omxop_dxy), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
venn(list(pmp_c_48=pmp_c_48,cmxcp_fixed=cmxcp_fixed, gem_cans= gem_cans, cmxcp_dxy=cmxcp_dxy), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
venn(list(pmp_o_48=pmp_o_48,omxop_fixed=omxop_fixed, gem_cans= gem_cans,omxop_dxy=omxop_dxy), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
intersect(intersect(pmp_c_8 ,gem_cans),intersect(cmxcp_dxy,cmxcp_fixed))
intersect(intersect(pmp_o_8 ,gem_cans),intersect(omxop_dxy,omxop_fixed))
intersect(intersect(pmp_c_48,gem_cans),intersect(cmxcp_dxy,cmxcp_fixed))
intersect(intersect(pmp_o_48,gem_cans),intersect(omxop_dxy,omxop_fixed))
xm_to_gene(intersect(intersect(pmp_c_8,gem_cans), intersect(cmxcp_dxy,cmxcp_fixed)))
xm_to_gene(intersect(intersect(pmp_o_8,gem_cans), intersect(omxop_dxy,omxop_fixed)))
xm_to_gene(intersect(intersect(pmp_c_48,gem_cans),intersect(cmxcp_dxy,cmxcp_fixed)))
xm_to_gene(intersect(intersect(pmp_o_48,gem_cans),intersect(omxop_dxy,omxop_fixed)))
length(unique(c(pmp_c_8,pmp_o_8,pmp_c_48,pmp_o_48)))
length(unique(c(intersect(pmp_c_8,intersect(cmxcp_dxy,cmxcp_fixed)))))
length(unique(c(intersect(pmp_c_48,intersect(cmxcp_dxy,cmxcp_fixed)))))
length(unique(c(intersect(pmp_o_8,intersect(omxop_dxy,omxop_fixed)))))
length(unique(c(intersect(pmp_o_48,intersect(omxop_dxy,omxop_fixed)))))
all1 <- c("XM_015404651.1","XM_015385242.1","XM_015399645.1","XM_015370331.1","XM_015370604.1","XM_015370640.1","XM_015370652.1","XM_015371134.1","XM_015375030.1","XM_015375570.1","XM_015376119.1","XM_015376708.1","XM_015378696.1","XM_015379768.1","XM_015380008.1","XM_015381079.1","XM_015382890.1","XM_015385225.1","XM_015387512.1","XM_015387941.1","XM_015389801.1","XM_015392646.1","XM_015395321.1","XM_015395713.1","XM_015398199.1","XM_015402116.1","XM_015404131.1","XM_015372096.1","XM_015373727.1","XM_015376991.1","XM_015377614.1","XM_015378775.1","XM_015379276.1","XM_015380414.1","XM_015382197.1","XM_015385349.1","XM_015386614.1","XM_015386615.1","XM_015393592.1","XM_015398785.1","XM_015396291.1")
all1 <- trimws(all1, which = c("both", "left", "right"))
xm_to_gene(all1)




# AI gene table
ai_genes <- c(oaxom_8_ai ,
              caxcp_8_ai ,
              oaxop_48_ai,
              oaxop_8_ai ,
              cmxcp_48_ai,
              cmxcp_8_ai ,
              omxop_48_ai,
              omxop_8_ai )
ai_genes <- c(intersect(oaxom_8_ai ,intersect(oaxom_dxy,oaxom_fixed)),
              intersect(caxcp_8_ai ,intersect(caxcp_dxy,caxcp_fixed)),
              intersect(oaxop_48_ai,intersect(oaxop_dxy,oaxop_fixed)),
              intersect(oaxop_8_ai ,intersect(oaxop_dxy,oaxop_fixed)),
              intersect(cmxcp_48_ai,intersect(cmxcp_dxy,cmxcp_fixed)),
              intersect(cmxcp_8_ai ,intersect(cmxcp_dxy,cmxcp_fixed)),
              intersect(omxop_48_ai,intersect(omxop_dxy,omxop_fixed)),
              intersect(omxop_8_ai ,intersect(omxop_dxy,omxop_fixed)))
ai_genes <- intersect(ai_genes, gem_cans)
ai_genes <- mrna[mrna$related_accession %in% ai_genes, ]
#write.table(ai_genes,"C:/Users/jmcgirr/Documents/all_2018_samples/interesting_genes_tables/adaptive_misexpression_genes_fst_dxy_gemma.txt", row.names = FALSE, quote = FALSE, sep = "\t")
ai_genes <- merge(ai_genes, final_features, all.x = TRUE, by = c("related_accession"))
ai_genes <- merge(ai_genes, blast_key, all.x = TRUE, by = c("product_accession"))
nrow(ai_genes)
head(ai_genes)
#write.table(ai_genes,"C:/Users/jmcgirr/Documents/all_2018_samples/interesting_genes_tables/adaptive_misexpression_genes_zeb.txt", row.names = FALSE, quote = FALSE,sep = "\t")
pmp_genes <- c(pmp_c_8,
               pmp_o_8,
               pmp_c_48,
               pmp_o_48)
pmp_genes <- c(intersect(pmp_c_8,intersect(cmxcp_dxy,cmxcp_fixed)),
               intersect(pmp_o_8,intersect(omxop_dxy,omxop_fixed)),
               intersect(pmp_c_48,intersect(cmxcp_dxy,cmxcp_fixed)),
               intersect(pmp_o_48,intersect(omxop_dxy,omxop_fixed)))
pmp_genes <- intersect(pmp_genes,gem_cans)

pmp_genes <- mrna[mrna$related_accession %in% pmp_genes, ]
#write.table(pmp_genes,"C:/Users/jmcgirr/Documents/all_2018_samples/interesting_genes_tables/parallel_mp_genes_fst_dxy_gemma.txt", row.names = FALSE, quote = FALSE, sep = "\t")
pmp_genes <- merge(pmp_genes, final_features, all.x = TRUE, by = c("related_accession"))
pmp_genes <- merge(pmp_genes, blast_key, all.x = TRUE, by = c("product_accession"))
nrow(pmp_genes)
head(pmp_genes)
#write.table(pmp_genes,"C:/Users/jmcgirr/Documents/all_2018_samples/interesting_genes_tables/parallel_mp_genes_zeb.txt", row.names = FALSE, quote = FALSE, sep = "\t")

# overlapping adaptive misexpressed genes
ai_genes <- c(intersect(cmxcp_8_ai,omxop_8_ai),intersect(caxcp_8_ai, oaxop_8_ai))
#write.table(ai_genes,"C:/Users/jmcgirr/Documents/all_2018_samples/interesting_genes_tables/adaptive_misexpression_genes_crp_and_osp.txt", row.names = FALSE, quote = FALSE, sep = "\t")


# parallel DE MP and AI and fixed and gemma
venn(list(pmp_o_8=pmp_o_8,omxop_fixed=omxop_fixed, gem_cans= gem_cans,omxop_8_ai=omxop_8_ai), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
venn(list(pmp_c_8=pmp_c_8,cmxcp_fixed=cmxcp_fixed, gem_cans= gem_cans,cmxcp_8_ai=cmxcp_8_ai), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)

# parallel DE MP and super ME
venn(list(pmp_o_8=pmp_o_8,omxop_fixed=omxop_fixed, gem_cans= gem_cans,opxom_8_SME=opxom_8_SME), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
venn(list(pmp_c_8=pmp_c_8,cmxcp_fixed=cmxcp_fixed, gem_cans= gem_cans,cmxcp_8_SME=cmxcp_8_SME), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
intersect(intersect(pmp_c_8,cmxcp_fixed),intersect(gem_cans,cmxcp_8_SME))
intersect(pmp_c_8,cmxcp_8_SME)

# parallel DE MP and ME MP
venn(list(pmp_o_8=pmp_o_8,opxom_8_ME=opxom_8_ME), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
venn(list(pmp_c_8=pmp_c_8,cmxcp_8_ME=cmxcp_8_ME), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
venn(list(pmp_o_48=pmp_o_48,opxom_48_ME=opxom_48_ME), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
venn(list(pmp_c_48=pmp_c_48,cmxcp_48_ME=cmxcp_48_ME), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
length(unique(c(intersect(pmp_o_8,opxom_8_ME),intersect(pmp_c_8,cmxcp_8_ME),intersect(pmp_o_48,opxom_48_ME),
                intersect(pmp_c_48,cmxcp_48_ME))))
#tiff("C:/Users/jmcgirr/Documents/all_2018_samples/manuscript_figs/parallel_mse_pie.tiff", width = 4, height = 4, units = 'in', res = 1000)
pie(c(0.9626866,0.0373134), labels = c("",""), col = c(grb,"#BC1F8A"))
#dev.off()
#tiff("C:/Users/jmcgirr/Documents/all_2018_samples/manuscript_figs/parallel_pie.tiff", width = 4, height = 4, units = 'in', res = 1000)
pie(c(43/1249,45/1249,1161/1249), labels = c("",""), col = c(grb,"#BC1F8A",yel))
#dev.off()

#####

library(venn)
library(reshape2)
st <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/ase_data/ase_and_mse/super_lenient_ase_unphased_mbased_transtest_by_individual_no_ase_all_parents.txt", stringsAsFactors = FALSE, header = TRUE, sep = "\t")
st <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/ase_data/ase_and_mse/THE_ONE_lenient_ase_unphased_mbased_transtest_by_individual_no_ase_all_parents.txt", stringsAsFactors = FALSE, header = TRUE, sep = "\t")
genes <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/ase_data/mrna.saf", stringsAsFactors = FALSE, header = TRUE, sep = "\t")
genes <- cbind(genes, colsplit(genes$GeneID, ";", c("related_accession", "symbol")))
final_features <- read.table("D:/Martin Lab/rna_2018/all_2018_samples/features_gff.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t", quote = "")
blast_key <- read.table("D:/Martin Lab/rna_2018/all_2018_samples/blast/cyprinodon_to_danio_one_way_best_hit_symbols.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
#cranial_genes <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/GO/cranial_skeletal_system_development_GO_1904888.txt", stringsAsFactors = FALSE, header = FALSE, sep = "\t")
cranial_genes <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/GO/skeletal_system_morphogenesis_GO_0048705.txt", stringsAsFactors = FALSE, header = FALSE, sep = "\t")
cranial_genes$zeb_gene_symbol_one_way <- tolower(cranial_genes$V1)
cranial_genes$symbol <- tolower(cranial_genes$V1)


cis         <- c()  
comp        <- c()  
trans       <- c()  
mis_ase_comp<- c()  
mis_ase_cis <- c()  
mis_trans <- c() 
row_nums_specialists <- c(7:15,22:29)

#for (row_num in c(1:nrow(st)))
for (row_num in row_nums_specialists)
{
  cis                      <- c(cis, strsplit(st$cis_genes[row_num],";")[[1]])
  comp                    <- c(comp, strsplit(st$comp_genes[row_num],";")[[1]])
  trans                  <- c(trans, strsplit(st$trans_genes[row_num],";")[[1]])
  mis_ase_comp    <- c(mis_ase_comp, strsplit(st$mis_ase_comp_genes[row_num],";")[[1]])
  mis_ase_cis     <- c(mis_ase_cis,  strsplit(st$mis_ase_cis_genes[row_num],";")[[1]])
  mis_trans     <- c(mis_trans,  strsplit(st$mis_trans_genes[row_num],";")[[1]])
  
  length(mis_ase_comp)
  length(mis_ase_cis)
  length(mis_trans  )
}
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

## # am crosses
## caxcm 
## cmxca 
## oaxom 
## omxoa 
## # ap crosses
## caxcp 
## oaxop 
## opxoa 
## # specialist
## cmxcp
## opxom 
## # outgroup c
## naxca 
## upxca 
## # lake cross
## oaxca 
## caxoa 
## cmxom 
## opxcp 

venn(list(caxcm_48_DE=caxcm_48_DE,caxcm_48_ME=caxcm_48_ME,cmxca_48_ME=cmxca_48_ME), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
venn(list(caxcm_8_DE=caxcm_8_DE,caxcm_8_ME=caxcm_8_ME,cmxca_8_ME=cmxca_8_ME), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
venn(list(oaxom_48_DE=oaxom_48_DE,oaxom_48_ME=oaxom_48_ME,omxoa_48_ME=omxoa_48_ME), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
venn(list(oaxom_8_DE=oaxom_8_DE,oaxom_8_ME=oaxom_8_ME,omxoa_8_ME=omxoa_8_ME), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)

venn(list(caxcp_48_DE=caxcp_48_DE,caxcp_48_ME=caxcp_48_ME), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
venn(list(caxcp_8_DE=caxcp_8_DE,caxcp_8_ME=caxcp_8_ME), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
venn(list(oaxop_48_DE=oaxop_48_DE,oaxop_48_ME=oaxop_48_ME,opxoa_48_ME=opxoa_48_ME), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
venn(list(oaxop_8_DE=oaxop_8_DE,oaxop_8_ME=oaxop_8_ME,opxoa_8_ME=opxoa_8_ME), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
venn(list(oaxop_48_DE=oaxop_48_DE,oaxop_48_ME=oaxop_48_ME,opxoa_48_ME=opxoa_48_ME,caxcp_48_DE=caxcp_48_DE,caxcp_48_ME=caxcp_48_ME), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
venn(list(oaxop_8_DE=oaxop_8_DE,oaxop_8_ME=oaxop_8_ME,opxoa_8_ME=opxoa_8_ME,caxcp_8_DE=caxcp_8_DE,caxcp_8_ME=caxcp_8_ME), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)


venn(list(cmxcp_48_DE=cmxcp_48_DE,cmxcp_48_ME=cmxcp_48_ME), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
venn(list(omxop_48_DE=omxop_48_DE,opxom_48_ME=opxom_48_ME), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
venn(list(cmxcp_8_DE=cmxcp_8_DE ,cmxcp_8_ME=cmxcp_8_ME), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
venn(list(omxop_8_DE=omxop_8_DE,opxom_8_ME=opxom_8_ME), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
venn(list(omxop_48_DE=omxop_48_DE,opxom_48_ME=opxom_48_ME,cmxcp_48_DE=cmxcp_48_DE,cmxcp_48_ME=cmxcp_48_ME), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)
venn(list(omxop_8_DE=omxop_8_DE,opxom_8_ME=opxom_8_ME,cmxcp_8_DE=cmxcp_8_DE ,cmxcp_8_ME=cmxcp_8_ME), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)

adaptive_incompatibilities_am <- c(intersect(intersect(caxcm_48_DE,caxcm_48_ME),cmxca_48_ME),
                                   intersect(intersect(caxcm_8_DE,caxcm_8_ME),cmxca_8_ME),
                                   intersect(intersect(oaxom_48_DE,oaxom_48_ME),omxoa_48_ME),
                                   intersect(intersect(oaxom_8_DE,oaxom_8_ME),omxoa_8_ME))

adaptive_incompatibilities_ap <- c(intersect(caxcp_48_DE,caxcp_48_ME),intersect(caxcp_8_DE,caxcp_8_ME),
                                   intersect(intersect(oaxop_48_DE,oaxop_48_ME),opxoa_48_ME),
                                   intersect(intersect(oaxop_8_DE,oaxop_8_ME),opxoa_8_ME))

adaptive_incompatibilities_mp <- c(intersect(cmxcp_48_DE,cmxcp_48_ME),intersect(omxop_48_DE,opxom_48_ME),
                                intersect(cmxcp_8_DE,cmxcp_8_ME),intersect(omxop_8_DE,opxom_8_ME))

adaptive_incompatibilities_48 <- c(intersect(intersect(caxcm_48_DE,caxcm_48_ME),cmxca_48_ME),
                                   intersect(intersect(oaxom_48_DE,oaxom_48_ME),omxoa_48_ME),
                                   intersect(intersect(oaxop_48_DE,oaxop_48_ME),opxoa_48_ME),
                                   intersect(caxcp_48_DE,caxcp_48_ME),
                                   intersect(cmxcp_48_DE,cmxcp_48_ME),
                                   intersect(omxop_48_DE,opxom_48_ME))

adaptive_incompatibilities_8 <- c(intersect(intersect(caxcm_8_DE,caxcm_8_ME),cmxca_8_ME),
                                   intersect(intersect(oaxom_8_DE,oaxom_8_ME),omxoa_8_ME),
                                   intersect(intersect(oaxop_8_DE,oaxop_8_ME),opxoa_8_ME),
                                   intersect(caxcp_8_DE,caxcp_8_ME),
                                   intersect(cmxcp_8_DE,cmxcp_8_ME),
                                   intersect(omxop_8_DE,opxom_8_ME))



venn(list(adaptive_incompatibilities_am=adaptive_incompatibilities_am,
          adaptive_incompatibilities_ap=adaptive_incompatibilities_ap,
          adaptive_incompatibilities_mp=adaptive_incompatibilities_mp), 
     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85)

adaptive_incompatibilities <- c(adaptive_incompatibilities_am,
                                adaptive_incompatibilities_ap,
                                adaptive_incompatibilities_mp)
xm_to_gene <- function(xms){
  
  xm_transcripts <- data.frame(related_accession = xms)
  xm_transcripts <- merge(xm_transcripts,final_features, by = c("related_accession"))
  xm_transcripts<- merge(xm_transcripts, blast_key, by = c("product_accession"))
  print("cyprinodon gene annotations")
  print(xm_transcripts$symbol)
  print("danio ortholog annotations")
  print(xm_transcripts$zeb_gene_symbol_one_way)
  print("danio orthologs annotated for skeletal effects")
  xm_cranial <- merge(xm_transcripts, cranial_genes, by = c("zeb_gene_symbol_one_way"))
  print(xm_cranial$zeb_gene_symbol_one_way)
  print("cyprinodon genes annotated for skeletal effects")
  xm_cranial_cyp <- merge(xm_transcripts, cranial_genes, by = c("symbol"))
  print(xm_cranial_cyp$symbol)
}

xm_to_gene(adaptive_incompatibilities_8)
intersect(adaptive_incompatibilities,mis_ase_cis)
intersect(adaptive_incompatibilities,mis_trans)
xms <- fixed_snps_neg_taj_ai_genes_all
xms <- ai_and_jaw_and_selection



library(gridExtra)
library(eulerr)

caxcm_48_ME <- intersect(caxcm_48_ME,cmxca_48_ME)
caxcm_8_ME <- intersect(caxcm_8_ME,cmxca_8_ME)
oaxom_48_ME <- intersect(oaxom_48_ME,omxoa_48_ME)
oaxom_8_ME <- intersect(oaxom_8_ME,omxoa_8_ME)
oaxop_48_ME <- intersect(oaxop_48_ME,opxoa_48_ME)
oaxop_8_ME <- intersect(oaxop_8_ME,opxoa_8_ME)


de_all <- list(caxcm_48_DE,
       caxcm_8_DE,
       caxcp_48_DE,
       caxcp_8_DE,
       oaxom_48_DE,
       oaxom_8_DE,
       oaxop_48_DE,
       oaxop_8_DE,
       cmxcp_48_DE,
       cmxcp_8_DE,
       omxop_48_DE,
       omxop_8_DE)
me_all <- list(caxcm_48_ME,
               caxcm_8_ME,
               caxcp_48_ME,
               caxcp_8_ME,
               oaxom_48_ME,
               oaxom_8_ME,
               oaxop_48_ME,
               oaxop_8_ME,
               cmxcp_48_ME,
               cmxcp_8_ME,
               opxom_48_ME,
               opxom_8_ME)
comp_names <- c("caxcm_48",
               "caxcm_8",
               "caxcp_48",
               "caxcp_8",
               "oaxom_48",
               "oaxom_8",
               "oaxop_48",
               "oaxop_8",
               "cmxcp_48",
               "cmxcp_8",
               "opxom_48",
               "opxom_8")
i <- 1
for (i in c(1:12)) 
{
  #png(paste("C:/Users/jmcgirr/Desktop/test_plots/",comp_names[i],"test.png", sep = ""), width = 5, height = 5, units = 'in', res = 500,bg = "transparent")
  #png(paste("C:/Users/jmcgirr/Documents/all_2018_samples/manuscript_figs/de_me_venns/",comp_names[i],"de_me.png", sep = ""), width = 5, height = 5, units = 'in', res = 500,bg = "transparent")
  
  de_genes <- de_all[[i]]
  me_genes <- me_all[[i]]
  bo_genes <- intersect(de_genes,me_genes)
  
  de <- length(de_genes)
  me <- length(me_genes)
  bo <- length(bo_genes)
  
  VennDiag <- euler(c("A" = de, "B" = me, "A&B" = bo), counts=TRUE)
  plt <- plot(VennDiag, alpha=1, lty = 1, col = "black",counts = TRUE,quantities = list(fontsize = 15),
       fill=c(red,blu), labels = c("",""))#, main = comp_names[i])
  print(plt)
  #dev.off()
  
  
}

#png("C:/Users/jmcgirr/Documents/all_2018_samples/manuscript_figs/de_me_venns/ap_both_de_me.png", width = 5, height = 5, units = 'in', res = 500,bg = "transparent")

VennDiag <- euler(c("A" = 84, "B" = 72, "A&B" = 3), counts=TRUE)
plt <- plot(VennDiag, alpha=1, lty = 1, col = "black",counts = TRUE,quantities = FALSE,
            fill=c("#BC1F8A","#BC1F8A"), labels = c("",""))#, main = comp_names[i])
plt
#dev.off()


###############################################################################################
############# overlap adaptive incompatibility genes with fixed snps and selective sweeps #####
###############################################################################################

##### overlap gemma and tajimas d #####

require(data.table)
library(reshape2)

mrna <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/ase_data/mrna.saf", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
mrna <- cbind(mrna, colsplit(mrna$GeneID, ";", c("related_accession", "gene_name")))
mrna_table <- mrna
colnames(mrna_table) <- c("GeneID","CHROM","START","END","Strand","related_accession", "gene_name")
mrna_table$GeneID <- NULL
mrna_table$START <- mrna_table$START -10000
mrna_table$END <- mrna_table$END +10000
setDT(mrna_table)
gem <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/GEMMA/candidates_up_jaw_95percentile.bed", header = FALSE, stringsAsFactors = FALSE)
colnames(gem) <- c("CHROM","START","END","PIP")
#gem <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/GEMMA/candidates_95_percentile/candidates_up_jaw_95percentile.genes.txt", header = FALSE, stringsAsFactors = FALSE)
#colnames(gem) <- c("CHROM","START","END","mrna","gene","zeb_gene")
setDT(gem)
setkey(mrna_table)
c <- foverlaps(gem, mrna_table,by.x = c("CHROM", "START", "END"),by.y = c("CHROM", "START", "END"),type="any", nomatch=0L)
c <- as.data.frame(c)
rownames(c) <- NULL
head(c)
gem_cans <- c$related_accession


crosses <- c("caxcp","oaxom", "oaxop", "cmxcp", "omxop")
crosses <- c("caxcm","caxcp","oaxom", "oaxop", "cmxcp", "omxop")

oaxom <- c(intersect(intersect(oaxom_8_DE,oaxom_8_ME),omxoa_8_ME))
caxcp <- intersect(caxcp_8_DE,caxcp_8_ME)
oaxop <- c(intersect(intersect(oaxop_48_DE,oaxop_48_ME),opxoa_48_ME),
           intersect(intersect(oaxop_8_DE,oaxop_8_ME),opxoa_8_ME))
cmxcp <- c(intersect(cmxcp_48_DE,cmxcp_48_ME),intersect(cmxcp_8_DE,cmxcp_8_ME))
omxop <- c(intersect(omxop_48_DE,opxom_48_ME),intersect(omxop_8_DE,opxom_8_ME))
ai_genes <- list(caxcp=caxcp,oaxom=oaxom, oaxop=oaxop, cmxcp=cmxcp, omxop=omxop)
#i <- 1
##### fixed snps and negative taj D

for (i in c(1:length(crosses)))
{
  cross <- crosses[i]
  p1 <- strsplit(cross, "x")[[1]][1]
  p2 <- strsplit(cross, "x")[[1]][2]
  
  taj_p1 <-    read.table(paste("C:/Users/jmcgirr/Documents/all_2018_samples/fst_dna/",p1,"_taj_d_20kb.txt.Tajima.D", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  taj_p2 <-    read.table(paste("C:/Users/jmcgirr/Documents/all_2018_samples/fst_dna/",p2,"_taj_d_20kb.txt.Tajima.D", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  
  taj_p1$END <- taj_p1$BIN_START + 19999
  taj_p1$START <- taj_p1$BIN_START
  taj_p2$END <- taj_p2$BIN_START + 19999
  taj_p2$START <- taj_p2$BIN_START
  
  
  # find genes in outlier negative tajima's d (10th %ile)  
  setDT(taj_p1)
  setkey(mrna_table)
  c <- foverlaps(taj_p1, mrna_table,by.x = c("CHROM", "START", "END"),by.y = c("CHROM", "START", "END"),type="any", nomatch=0L)
  c <- as.data.frame(c)
  rownames(c) <- NULL
  #thresh <- mean(taj_p1$TajimaD, na.rm = TRUE) - sd(taj_p1$TajimaD, na.rm = TRUE)
  thresh <- quantile(na.omit(taj_p1$TajimaD), 0.1)
  #thresh <- -0.5
  c <- c[which(c$TajimaD <= thresh),]
  neg_taj_p1 <- c$related_accession
  c <- merge(c,final_features,all.x = TRUE, by = c("related_accession"))
  c<- merge(c, blast_key,all.x = TRUE, by = c("product_accession"))
  c$name <- NULL
  c$tag <- NULL
  
  write.csv(c,paste("C:/Users/jmcgirr/Documents/all_2018_samples/interesting_genes_tables/",p1,"_taj_90.csv",sep=""), row.names = FALSE, quote = FALSE)
  
  setDT(taj_p2)
  setkey(mrna_table)
  c <- foverlaps(taj_p2, mrna_table,by.x = c("CHROM", "START", "END"),by.y = c("CHROM", "START", "END"),type="any", nomatch=0L)
  c <- as.data.frame(c)
  rownames(c) <- NULL
  thresh <- quantile(na.omit(taj_p2$TajimaD), 0.1)
  c <- c[which(c$TajimaD <= thresh),]
  neg_taj_p2 <- c$related_accession
  c <- merge(c,final_features,all.x = TRUE, by = c("related_accession"))
  c<- merge(c, blast_key,all.x = TRUE, by = c("product_accession"))
  c$name <- NULL
  c$tag <- NULL
  
  write.csv(c,paste("C:/Users/jmcgirr/Documents/all_2018_samples/interesting_genes_tables/",p2,"_taj_90.csv",sep=""), row.names = FALSE, quote = FALSE)
  
  neg_taj_p1_or_p2 <- unique(c(neg_taj_p1,neg_taj_p2))
  
  
  fst <- read.table(paste("C:/Users/jmcgirr/Documents/all_2018_samples/fst_dna/",cross,"_fst",sep = ""), header = TRUE, stringsAsFactors = FALSE)
  fst$START <- fst$POS 
  fst$END <- fst$POS + 1
  fst <- fst[which(fst$WEIR_AND_COCKERHAM_FST == 1),]
  
  setDT(fst)
  setkey(mrna_table)
  c <- foverlaps(fst, mrna_table,by.x = c("CHROM", "START", "END"),by.y = c("CHROM", "START", "END"),type="any", nomatch=0L)
  c <- as.data.frame(c)
  rownames(c) <- NULL
  fixed <- c$related_accession
  c <- merge(c,final_features,all.x = TRUE, by = c("related_accession"))
  c<- merge(c, blast_key,all.x = TRUE, by = c("product_accession"))
  c$name <- NULL
  c$tag <- NULL
  
  write.csv(c,paste("C:/Users/jmcgirr/Documents/all_2018_samples/interesting_genes_tables/",cross,"_fixed_snp_genes.csv",sep=""), row.names = FALSE, quote = FALSE)

  
  #ai <- data.frame(ai_genes[cross])
  
  #venn(list(ai_genes=ai[,1],near_fixed_snp=fixed,gemma=gem_cans,neg_taj_p1_or_p2=neg_taj_p1_or_p2), 
  #     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85, title = cross)
  #venn(list(ai_genes=ai[,1],near_fixed_snp=fixed,gemma=gem_cans,sig_sweed_p1_or_p2=sig_sweed_p1_or_p2), 
  #     ilabels = TRUE, cexsn = 0.70, zcolor = "style", cexil = 0.85, title = cross)
  
  dxy <- read.csv(paste("D:/Martin Lab/rna_2018/fst/dna/",cross,"_popgen_dna_stats_corr_dxy.csv",sep=""), header = TRUE, stringsAsFactors = FALSE)
  colnames(dxy)[colnames(dxy )=="start"] <- "START"
  colnames(dxy)[colnames(dxy )=="end"] <- "END"
  colnames(dxy)[colnames(dxy )=="scaffold"] <- "CHROM"

  setDT(dxy)
  setkey(mrna_table)
  c <- foverlaps(dxy, mrna_table,by.x = c("CHROM", "START", "END"),by.y = c("CHROM", "START", "END"),type="any", nomatch=0L)
  c <- as.data.frame(c)
  rownames(c) <- NULL
  thresh <- quantile(na.omit(dxy$corr_dxy), 0.90)
  c <- c[which(c$corr_dxy >= thresh),]
  c <- merge(c,final_features,all.x = TRUE, by = c("related_accession"))
  c<- merge(c, blast_key,all.x = TRUE, by = c("product_accession"))
  c$name <- NULL
  c$tag <- NULL
  
  write.csv(c,paste("C:/Users/jmcgirr/Documents/all_2018_samples/interesting_genes_tables/",cross,"_corr_dxy_90.csv",sep=""), row.names = FALSE, quote = FALSE)
  

}


fixed_snps_neg_taj_ai_genes_table <- data.frame(comparison = crosses,
                                        fixed_snps_neg_taj_ai_genes = fixed_snps_neg_taj_ai_genes,
                                        stringsAsFactors = FALSE)


fixed_snps_ai_genes_table <- data.frame(comparison = crosses,
                                        fixed_snps_ai_genes = fixed_snps_ai_genes,
                                        stringsAsFactors = FALSE)       
neg_taj_ai_genes_table    <-  data.frame(comparison = crosses,
                                         neg_taj_ai_genes = neg_taj_ai_genes,
                                         stringsAsFactors = FALSE) 

all_tables <- merge(fixed_snps_neg_taj_ai_genes_table, fixed_snps_ai_genes_table, by = c("comparison"))
all_tables <- merge(all_tables, neg_taj_ai_genes_table, by = c("comparison"))
#write.table(all_tables,        "C:/Users/jmcgirr/Documents/all_2018_samples/ase_data/ase_and_mse/fixed_snps_neg_taj_ai_genes_table_separated.txt", row.names = FALSE, quote = FALSE, sep ="\t")

##### fixed snps and negative taj D and sweed sig
#fixed_snps_neg_taj_ai_genes <- c()
#crosses <- c("caxcp","oaxom", "oaxop", "cmxcp", "omxop")
#
#for (i in c(1:length(crosses)))
#{
  cross <- crosses[i]
  p1 <- strsplit(cross, "x")[[1]][1]
  p2 <- strsplit(cross, "x")[[1]][2]
  
  fst <- read.table(paste("C:/Users/jmcgirr/Documents/all_2018_samples/fst_dna/",cross,"_fst",sep = ""), header = TRUE)
  taj_p1 <-    read.table(paste("C:/Users/jmcgirr/Documents/all_2018_samples/fst_dna/",p1,"_taj_d.txt.Tajima.D", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  taj_p2 <-    read.table(paste("C:/Users/jmcgirr/Documents/all_2018_samples/fst_dna/",p2,"_taj_d.txt.Tajima.D", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  sweed_p1 <-    read.table(paste("C:/Users/jmcgirr/Documents/all_2018_samples/sweed/",p1,"_pop_bottle_58_sweeps.txt", sep = ""), header = FALSE, stringsAsFactors = FALSE)
  sweed_p2 <-    read.table(paste("C:/Users/jmcgirr/Documents/all_2018_samples/sweed/",p2,"_pop_bottle_58_sweeps.txt", sep = ""), header = FALSE, stringsAsFactors = FALSE)
  
  # set fst and tajimas D threshold
  # fixed snps and negative D
  
  fst <- fst[which(fst$WEIR_AND_COCKERHAM_FST == 1),]
  taj_p1 <- taj_p1[which(taj_p1$TajimaD < 0),]
  taj_p2 <- taj_p2[which(taj_p2$TajimaD < 0),]
  
  ai <- data.frame(ai_genes[cross])
  can_genes <- c()
  #j <- 1
  for (j in c(1:nrow(ai)))
  {
    
    mrna <- genes[which(genes$related_accession == ai[j,1]),]
    snp_region <- c(mrna$Start[1]-10000,mrna$End[1]+10000)
    fst_mrna <- fst[which(fst$CHROM == mrna$Chr[1] & fst$POS > snp_region[1] & fst$POS < snp_region[2]),]
    taj_p1_mrna <- taj_p1[which(taj_p1$BIN_START > (snp_region[1]-10000) & taj_p1$BIN_START < (snp_region[1]+10000)),]
    taj_p2_mrna <- taj_p2[which(taj_p2$BIN_START > (snp_region[1]-10000) & taj_p2$BIN_START < (snp_region[1]+10000)),]
    sweed_p1_mrna <- sweed_p1[which(sweed_p1$V1 == mrna$Chr[1]),]
    sweed_p2_mrna <- sweed_p2[which(sweed_p2$V1 == mrna$Chr[1]),]
    sweed_p1_sig <- quantile(sweed_p1_mrna$V3, .90)[[1]]
    sweed_p2_sig <- quantile(sweed_p2_mrna$V3, .90)[[1]]
    sweed_p1_mrna <- sweed_p1_mrna[which(sweed_p1_mrna$V2 > (snp_region[1]-10000) & sweed_p1_mrna$V2 < (snp_region[1]+10000) & sweed_p1_mrna$V4 >= sweed_p1_sig),]
    sweed_p2_mrna <- sweed_p2_mrna[which(sweed_p2_mrna$V2 > (snp_region[1]-10000) & sweed_p2_mrna$V2 < (snp_region[1]+10000) & sweed_p2_mrna$V4 >= sweed_p2_sig),]
    
    if(nrow(fst_mrna)>0 & (nrow(taj_p1_mrna)>0|nrow(taj_p2_mrna)>0) & (nrow(sweed_p1_mrna)>0|nrow(sweed_p2_mrna)>0))
    {
      can_genes <- c(can_genes,mrna$related_accession[1])
    }
  }
  
  fixed_snps_neg_taj_ai_genes<- c(fixed_snps_neg_taj_ai_genes,paste(can_genes, sep = ";",collapse=";"))
  
  
}
#fixed_snps_neg_taj_ai_genes
##fixed_snps_neg_taj_ai_genes_table <- data.frame(comparison = crosses,
##                                        fixed_snps_neg_taj_ai_genes = fixed_snps_neg_taj_ai_genes,
##                                        stringsAsFactors = FALSE)
#
##write.table(fixed_snps_neg_taj_ai_genes_table,        "C:/Users/jmcgirr/Documents/all_2018_samples/ase_data/ase_and_mse/fixed_snps_neg_taj_ai_genes_table.txt", row.names = FALSE, quote = FALSE, sep ="\t")
#
#
#fixed_snps_neg_taj_ai_genes_all <- c("XM_015371309.1","XM_015401310.1","XM_015400530.1","XM_015378973.1","XM_015402122.1","XM_015404671.1","XM_015398227.1","XM_015402903.1","XM_015388091.1","XM_015403848.1","XM_015389560.1","XM_015377048.1","XM_015372161.1","XM_015371166.1","XM_015377824.1","XM_015370806.1","XM_015402115.1","XM_015394568.1","XM_015374217.1","XM_015396828.1","XM_015388332.1","XM_015379795.1","XM_015389317.1","XM_015378526.1","XM_015380427.1","XM_015403156.1","XM_015384312.1","XM_015388601.1","XM_015398326.1","XM_015381657.1","XM_015394934.1","XM_015377629.1","XM_015371423.1","XM_015393320.1","XM_015375370.1","XM_015378968.1","XM_015388668.1","XM_015387654.1","XM_015403298.1","XM_015387876.1","XM_015392140.1","XM_015399576.1","XM_015386511.1","XM_015369927.1","XM_015386306.1","XM_015384972.1","XM_015389807.1","XM_015372823.1","XM_015400295.1","XM_015378331.1","XM_015372614.1","XM_015394388.1","XM_015396414.1","XM_015372485.1","XM_015374552.1","XM_015390249.1","XM_015378342.1","XM_015386804.1","XM_015393506.1","XM_015379919.1","XM_015373765.1","XM_015389438.1","XM_015403216.1","XM_015369602.1","XM_015381070.1","XM_015383682.1","XM_015390740.1","XM_015401302.1","XM_015395051.1","XM_015398647.1","XM_015384652.1","XM_015390234.1","XM_015384373.1","XM_015372979.1","XM_015393272.1","XM_015380884.1","XM_015391433.1","XM_015379462.1","XM_015372595.1","XM_015373262.1","XM_015377864.1","XM_015392918.1","XM_015377701.1","XM_015388906.1","XM_015393460.1","XM_015397741.1","XM_015369916.1","XM_015394622.1","XM_015376842.1","XM_015384851.1","XM_015394567.1","XM_015380941.1","XM_015372477.1","XM_015391083.1","XM_015379690.1","XM_015392458.1","XM_015381839.1","XM_015402488.1",
#                                     "XM_015392920.1","XM_015392909.1","XM_015375205.1","XM_015373543.1","XM_015400182.1","XM_015370560.1","XM_015370867.1","XM_015381826.1","XM_015378276.1","XM_015401963.1","XM_015371438.1","XM_015381007.1","XM_015374547.1","XM_015394812.1","XM_015387154.1","XM_015393317.1","XM_015389364.1","XM_015398568.1","XM_015372001.1","XM_015373000.1","XM_015402320.1","XM_015404364.1","XM_015370214.1","XM_015392416.1","XM_015370650.1","XM_015392909.1","XM_015381015.1","XM_015394487.1","XM_015400705.1","XM_015404671.1","XM_015400851.1","XM_015404394.1","XM_015380597.1","XM_015388364.1","XM_015402096.1","XM_015381663.1","XM_015405336.1","XM_015386614.1","XM_015378312.1","XM_015392887.1","XM_015403147.1","XM_015395339.1","XM_015380986.1","XM_015376203.1","XM_015382551.1","XM_015402993.1","XM_015398425.1","XM_015395512.1","XM_015384989.1","XM_015375968.1","XM_015378212.1","XM_015384345.1","XM_015399552.1","XM_015397648.1","XM_015395626.1","XM_015391768.1","XM_015370278.1","XM_015379798.1","XM_015370896.1","XM_015394645.1","XM_015404053.1","XM_015389787.1","XM_015381072.1","XM_015379246.1","XM_015391944.1","XM_015380636.1","XM_015397082.1","XM_015371310.1","XM_015372278.1","XM_015403625.1","XM_015383707.1","XM_015369058.1","XM_015378004.1","XM_015403788.1","XM_015383401.1","XM_015401791.1","XM_015403848.1","XM_015381159.1","XM_015400835.1","XM_015372123.1","XM_015370249.1","XM_015381215.1","XM_015396313.1","XM_015401231.1","XM_015404250.1","XM_015394620.1","XM_015369894.1","XM_015394812.1","XM_015398330.1","XM_015369960.1","XM_015386468.1","XM_015388212.1","XM_015385070.1","XM_015381323.1","XM_015384304.1","XM_015390008.1","XM_015391639.1","XM_015385092.1",
#                                     "XM_015376003.1","XM_015371199.1","XM_015401907.1","XM_015403491.1","XM_015395742.1","XM_015384645.1","XM_015380698.1","XM_015384865.1","XM_015375262.1","XM_015403322.1","XM_015390673.1","XM_015385630.1","XM_015401458.1","XM_015394246.1","XM_015383469.1","XM_015400637.1","XM_015394513.1","XM_015374248.1","XM_015404213.1","XM_015374119.1","XM_015369954.1","XM_015379278.1","XM_015373295.1","XM_015401227.1","XM_015381007.1","XM_015384780.1","XM_015373412.1","XM_015376342.1","XM_015370883.1","XM_015405113.1","XM_015380040.1","XM_015384753.1","XM_015388344.1","XM_015402474.1","XM_015404303.1","XM_015390244.1","XM_015397277.1","XM_015377809.1","XM_015402355.1","XM_015394482.1","XM_015375372.1","XM_015401607.1","XM_015401944.1","XM_015391900.1","XM_015383508.1","XM_015371807.1","XM_015376523.1","XM_015381436.1","XM_015377065.1","XM_015370913.1","XM_015379783.1","XM_015369741.1","XM_015377296.1","XM_015399020.1","XM_015373683.1","XM_015399122.1","XM_015391631.1","XM_015405177.1","XM_015396687.1","XM_015397247.1","XM_015394175.1","XM_015369415.1","XM_015379514.1","XM_015395327.1","XM_015374387.1","XM_015376194.1","XM_015393235.1","XM_015376245.1","XM_015392582.1","XM_015398874.1","XM_015401797.1","XM_015404369.1","XM_015385606.1","XM_015397675.1","XM_015383425.1","XM_015388814.1","XM_015399685.1","XM_015383709.1","XM_015376150.1","XM_015376168.1","XM_015396544.1","XM_015394526.1","XM_015388362.1","XM_015401135.1","XM_015380662.1","XM_015373181.1","XM_015375713.1","XM_015381839.1","XM_015385063.1","XM_015369861.1","XM_015377792.1","XM_015386141.1","XM_015377529.1")

intersect(fixed_snps_neg_taj_ai_genes_all,mis_ase_cis)
intersect(fixed_snps_neg_taj_ai_genes_all,mis_trans)

length(fixed_snps_neg_taj_ai_genes_all)
xm_to_gene(fixed_snps_neg_taj_ai_genes_all)
xm_to_gene(intersect(fixed_snps_neg_taj_ai_genes_all,mis_ase_cis))

ai_and_selection <- intersect(fixed_snps_neg_taj_ai_genes_all,adaptive_incompatibilities)
ai_and_jaw_and_reg_mechanism <- c(intersect(candidates_95$V4,intersect(adaptive_incompatibilities,mis_ase_cis)),
                                  intersect(candidates_95$V4,intersect(adaptive_incompatibilities,mis_trans)))
ai_and_jaw_and_selection <- intersect(candidates_95$V4,intersect(fixed_snps_neg_taj_ai_genes_all,adaptive_incompatibilities))
ai_and_jaw_and_selection_and_reg_mechanism <- intersect(ai_and_jaw_and_selection,ai_and_jaw_and_reg_mechanism)

xm_to_gene(ai_and_jaw_and_selection)


library(gridExtra)
library(eulerr)
xm_to_gene <- function(xms){
  
  xm_transcripts <- data.frame(related_accession = xms)
  xm_transcripts <- merge(xm_transcripts,final_features, by = c("related_accession"))
  xm_transcripts<- merge(xm_transcripts, blast_key, by = c("product_accession"))
  print("cyprinodon gene annotations")
  print(xm_transcripts$symbol)
  print("danio ortholog annotations")
  print(xm_transcripts$zeb_gene_symbol_one_way)
  print("danio orthologs annotated for skeletal effects")
  xm_cranial <- merge(xm_transcripts, cranial_genes, by = c("zeb_gene_symbol_one_way"))
  print(xm_cranial$zeb_gene_symbol_one_way)
  print("cyprinodon genes annotated for skeletal effects")
  xm_cranial_cyp <- merge(xm_transcripts, cranial_genes, by = c("symbol"))
  print(xm_cranial_cyp$symbol)
}

fixed_snps_neg_taj_ai_genes_table <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/ase_data/ase_and_mse/fixed_snps_neg_taj_ai_genes_table_separated.txt", header = TRUE, stringsAsFactors = FALSE)
candidates_95_up <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/GEMMA/candidates_95_percentile/candidates_up_jaw_95percentile.genes.txt", header=F, stringsAsFactors = F)
candidates_95_low <- read.table("C:/Users/jmcgirr/Documents/all_2018_samples/GEMMA/candidates_95_percentile/candidates_low_jaw_95percentile.genes.txt", header=F, stringsAsFactors = F)
candidates_99_up <- read.delim("C:/Users/jmcgirr/Documents/all_2018_samples/GEMMA/candidates_up_jaw_99percentile.genes.txt", header=F, stringsAsFactors = F, sep = "\t")
candidates_99_low <- read.delim("C:/Users/jmcgirr/Documents/all_2018_samples/GEMMA/candidates_low_jaw_99percentile.genes.txt", header=F, stringsAsFactors = F)

caxcp_cans <- intersect(unlist(strsplit(fixed_snps_neg_taj_ai_genes_table[1,2], split=";")),candidates_95_up$V4)
oaxom_cans <- intersect(unlist(strsplit(fixed_snps_neg_taj_ai_genes_table[2,2], split=";")),candidates_95_up$V4)
oaxop_cans <- intersect(unlist(strsplit(fixed_snps_neg_taj_ai_genes_table[3,2], split=";")),candidates_95_up$V4)
cmxcp_cans <- intersect(unlist(strsplit(fixed_snps_neg_taj_ai_genes_table[4,2], split=";")),candidates_95_up$V4)
omxop_cans <- intersect(unlist(strsplit(fixed_snps_neg_taj_ai_genes_table[5,2], split=";")),candidates_95_up$V4)

#caxcp_cans <- intersect(unlist(strsplit(fixed_snps_neg_taj_ai_genes_table[1,2], split=";")),candidates_99_up$V5)
#oaxom_cans <- intersect(unlist(strsplit(fixed_snps_neg_taj_ai_genes_table[2,2], split=";")),candidates_99_up$V5)
#oaxop_cans <- intersect(unlist(strsplit(fixed_snps_neg_taj_ai_genes_table[3,2], split=";")),candidates_99_up$V5)
#cmxcp_cans <- intersect(unlist(strsplit(fixed_snps_neg_taj_ai_genes_table[4,2], split=";")),candidates_99_up$V5)
#omxop_cans <- intersect(unlist(strsplit(fixed_snps_neg_taj_ai_genes_table[5,2], split=";")),candidates_99_up$V5)


#caxcp_cans <- unlist(strsplit(fixed_snps_neg_taj_ai_genes_table[1,2], split=";"))
#oaxom_cans <- unlist(strsplit(fixed_snps_neg_taj_ai_genes_table[2,2], split=";"))
#oaxop_cans <- unlist(strsplit(fixed_snps_neg_taj_ai_genes_table[3,2], split=";"))
#cmxcp_cans <- unlist(strsplit(fixed_snps_neg_taj_ai_genes_table[4,2], split=";"))
#omxop_cans <- unlist(strsplit(fixed_snps_neg_taj_ai_genes_table[5,2], split=";"))

oaxom_8_ai <- intersect(intersect(oaxom_8_DE,oaxom_8_ME),omxoa_8_ME)
caxcp_8_ai <- intersect(caxcp_8_DE,caxcp_8_ME)
oaxop_48_ai <- intersect(intersect(oaxop_48_DE,oaxop_48_ME),opxoa_48_ME)
oaxop_8_ai  <- intersect(intersect(oaxop_8_DE,oaxop_8_ME),opxoa_8_ME)
cmxcp_48_ai <- intersect(cmxcp_48_DE,cmxcp_48_ME)
cmxcp_8_ai <-intersect(cmxcp_8_DE,cmxcp_8_ME)
omxop_48_ai <- intersect(omxop_48_DE,opxom_48_ME)
omxop_8_ai <- intersect(omxop_8_DE,opxom_8_ME)



cans_all <- list(oaxom_cans,
                 caxcp_cans,
                 oaxop_cans,
                 oaxop_cans,
                 cmxcp_cans,
                 cmxcp_cans,
                 omxop_cans,
                 omxop_cans)

ai_all <- list(oaxom_8_ai, 
               caxcp_8_ai, 
               oaxop_48_ai,
               oaxop_8_ai, 
               cmxcp_48_ai,
               cmxcp_8_ai, 
               omxop_48_ai,
               omxop_8_ai)

comp_names <- c("oaxom_8", 
                "caxcp_8", 
                "oaxop_48",
                "oaxop_8", 
                "cmxcp_48",
                "cmxcp_8", 
                "omxop_48",
                "omxop_8")

bo_geness <- c()
comps <- c()
#for (i in c(1:8)) 
i <- 2
{
  #png(paste("C:/Users/jmcgirr/Desktop/test_plots/",comp_names[i],"test_ai.png", sep = ""), width = 4, height = 4, units = 'in', res = 1000,bg = "transparent")
  png(paste("C:/Users/jmcgirr/Documents/all_2018_samples/manuscript_figs/",comp_names[i],"de_me.png", sep = ""), width = 4, height = 4, units = 'in', res = 1000,bg = "transparent")
  
  ai_genes <- ai_all[[i]]
  can_genes <- cans_all[[i]]
  bo_genes <- intersect(ai_genes,can_genes)
  bo_geness <- c(bo_geness,paste(bo_genes, sep = ";",collapse=";"))
  #bo_geness <- c(bo_geness,bo_genes)
  comps <- c(comps,comp_names[i])
  
  ai <- length(ai_genes)
  cans <- length(can_genes)
  bo <- length(bo_genes)
  
  VennDiag <- euler(c("A" = ai, "B" = cans, "A&B" = bo), counts=TRUE)
  plt <- plot(VennDiag, alpha=1, lty = 1, col = "black",counts = TRUE, quantities = list(fontsize = 15),
              fill=c(blu,blu), labels = c("",""), main = comp_names[i])
  print(plt)
  dev.off()
  
  
}
best_ai_genes    <-  data.frame(comparison = comps,
                                         ai_genes = bo_geness,
                                         stringsAsFactors = FALSE)
#write.table(best_ai_genes, "C:/Users/jmcgirr/Documents/all_2018_samples/ase_data/ase_and_mse/fixed_snps_neg_taj_ai_genes_best_candidates.txt", row.names = FALSE, quote = FALSE, sep ="\t")

# caxcp
jaw <- candidates_95_up$V4
#jaw <- candidates_99_up$V5
ai <- intersect(caxcp_8_DE,caxcp_8_ME)
# use 'separated' fixed_snps_neg_taj_ai
taj <- unlist(strsplit(fixed_snps_neg_taj_ai_genes_table[1,4], split=";"))
fst <- unlist(strsplit(fixed_snps_neg_taj_ai_genes_table[1,3], split=";"))
cols <- c(pur,yel,lir,grb)
venn(list(ai=ai,jaw=jaw,fst=fst,taj=taj), 
     ilabels = TRUE, cexsn = 0, zcolor = cols, cexil = 1.8,opacity = .7)
