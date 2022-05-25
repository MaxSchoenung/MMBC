###################################################################
#
# Script: 14_CREs_expression_2022-03.R
# Correlation transcriptional programs which are annotated to CREs
# Date: 2021-05-14
# Author: Maximilian Schönung
#
###################################################################

# Load Libraries ----------------------------------------------------------
#remotes::install_github("epigen/RnBeads",ref="feature/MouseMethylationBeadChip")
library(RnBeads)
library(RnBeads.mm10)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(DESeq2)
library(biomaRt)
library(RColorBrewer)
library(VennDiagram)
library(UpSetR)

# Set Directories ---------------------------------------------------------
setwd("/icgc/analysis/OE0565_projects/mmbc_hierarchy/")
data.dir <- "/icgc/analysis/OE0565_projects/mmbc_hierarchy/data/"
save.dir <- "/icgc/analysis/OE0565_projects/mmbc_hierarchy/RDS/differential/"
plot.dir <- "/icgc/analysis/OE0565_projects/mmbc_hierarchy/plots/"
table.dir <- "/icgc/analysis/OE0565_projects/mmbc_hierarchy/tables/"


# Load the annotation from stephen ----------------------------------------
annot <- read.delim(paste0(table.dir,"gene-annos_primary_one-row.bed"),stringsAsFactors = F)
colnames(annot)[1] <- "Chromosome"
rownames(annot) <- annot$name

# Load the CREs -------------------------------------------------------
mdCRE <- readRDS(paste0(data.dir,"2022-03-08_overlapping_CRE.RDS"))
chromCRE <- readRDS(paste0(data.dir,"2022-03-08_previous_CRE.RDS"))
novel_mdCRE <- readRDS(paste0(data.dir,"2022-03-08_novel_mdCRE.RDS"))
cre.list <- list("univ_mdCRE"=mdCRE,"chromCRE"=chromCRE,"novel_mdCRE"=novel_mdCRE)

# Load the RNA Seq data ---------------------------------------------------
dds <- readRDS(paste0(data.dir,"haemopedia_dds_mmbc.RDS"))
rld <- vst(dds,blind = F)
mat <- assay(rld)
colnames(mat) <- dds@colData$maxi_mmbc

# Biomart for gene names --------------------------------------------------
ensembl <- useMart('ensembl', dataset = 'mmusculus_gene_ensembl',host="http://nov2020.archive.ensembl.org")
genes <- getBM(
  attributes = c('external_gene_name','ensembl_gene_id'),
  filters = 'ensembl_gene_id',
  values = rownames(mat),
  mart = ensembl)%>%
  na_if("")%>% 
  na.omit()

# Function for CRE set to gene names --------------------------------------
cor_cre <- function(x){
gnames <- unique(unlist(strsplit(annot[x,"gene_name"], ",")))
gnames_ens <- genes[genes$external_gene_name%in%gnames,"ensembl_gene_id"]
return(mat[gnames_ens,])
}

# calculate the expression correlation for the CREs
cor.list <- lapply(cre.list,cor_cre)

# Plot the data -----------------------------------------------------------
colours <- c("LSK"="gray50",
             "CMP"="darkorange",
             "GMP"="darkgoldenrod2",
             "MEP"="orangered4",
             "Monocyte"="gold1",
             "Neutrophil"="darksalmon",
             "B-cell"="navy",
             "CD8_t-cell"="lightsteelblue1",
             "CD4_t-cell"="lightskyblue")

anno.df <- dds@colData[,"maxi_mmbc",drop=F]
rownames(anno.df) <- dds@colData$maxi_mmbc

order.vec <- c(9:10,15:16,23:26,11:14,17:22,1:8,27:30)
pdf(paste0(plot.dir,Sys.Date(),"_cCRE_corr_exp.pdf"))
for(i in 1:length(cor.list)){
  print(
    pheatmap(
      cor(cor.list[[i]])[order.vec,order.vec],
      color = colorRampPalette(c("white","lavender","navy"))(length(seq(0.7, 1, by = 0.001))),
      breaks=seq(0.7, 1, by = 0.001),
      cluster_rows = F,
      cluster_cols = F,
      border_color = NA,
      show_colnames = F,
      show_rownames=F,
      annotation_col = anno.df,
      annotation_row = anno.df,
      annotation_colors = list("maxi_mmbc"=colours),
      main=names(cor.list)[i]))
  print(
    pheatmap(
      cor(cor.list[[i]])[order.vec,order.vec],
      color = colorRampPalette(c("white","lavender","navy"))(length(seq(0.0, 1, by = 0.001))),
      breaks=seq(0.0, 1, by = 0.001),
      cluster_rows = F,
      cluster_cols = F,
      border_color = NA,
      show_colnames = F,
      show_rownames=F,
      annotation_col = anno.df,
      annotation_row = anno.df,
      annotation_colors = list("maxi_mmbc"=colours),
      main=names(cor.list)[i]))
}
dev.off()

#what is the overlap of these genes?
subject.venn <- lapply(cor.list,rownames)

pdf(paste0(plot.dir,Sys.Date(),"_CRE_gnames_Venn.pdf"))
grid.draw(venn.diagram(subject.venn,filename = NULL,
                       fill = c("#999999", "#E69F00", "#56B4E9")))
dev.off()

pdf(paste0(plot.dir,Sys.Date(),"_CRE_gnames_upset.pdf"))
upset(fromList(subject.venn), order.by = "freq")
dev.off()
