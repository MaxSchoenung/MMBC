###################################################################
#
# Script: 02_exploratory.R
# Exploratory analysis of mouse methylation array data
# Date: 2021-03-15
# Author: Maximilian Schönung
#
###################################################################

# Load Libraries ----------------------------------------------------------
#remotes::install_github("epigen/RnBeads",ref="feature/MouseMethylationBeadChip")
library(RnBeads)
library(RnBeads.mm10)
library(ggplot2)
library(ggpointdensity)
library(viridis)
library(ggpubr)


# Set Directories ---------------------------------------------------------
setwd("/icgc/analysis/OE0565_projects/mmbc_hierarchy/")
data.dir <- "/icgc/analysis/OE0565_projects/mmbc_hierarchy/data/"
plot.dir <- "/icgc/analysis/OE0565_projects/mmbc_hierarchy/plots/"
table.dir <- "/icgc/analysis/OE0565_projects/mmbc_hierarchy/tables/"


# Set the colours ---------------------------------------------------------

colours <- c("LSK"="gray50",
             "CMP"="darkorange",
             "GMP"="darkgoldenrod2",
             "MEP"="orangered4",
             "Monocyte"="gold1",
             "Neutrophil"="darksalmon",
             "B-cell"="navy",
             "CD8_t-cell"="lightsteelblue1",
             "CD4_t-cell"="lightskyblue")

# Load the dataset --------------------------------------------------------
rnb <- readRDS("/icgc/analysis/OE0565_projects/mmbc_hierarchy/data/2021-03-15_rnbeadsset_betas_scaling.RDS")

## access the meth and pheno data
pheno.rnb <- pheno(rnb)
pheno.rnb$Sample_Name
meth.rnb <- meth(rnb,row.names=T)
colnames(meth.rnb) <- pheno.rnb$Sample_Name

## focus on the cell types
meth.ct <- meth.rnb[,1:27]
pheno.ct <- pheno.rnb[1:27,]

# PCA ---------------------------------------------------------------------
meth.t <- t(meth.ct)
pca <- prcomp(meth.t,center = T,scale. = T)

pdf(paste0(plot.dir,Sys.Date(),"_all_cpg_pca.pdf"))
ggplot(as.data.frame(pca$x),aes(PC1,PC2,colour=pheno.ct$CellType))+
  geom_point()+
  theme_classic()
ggplot(as.data.frame(pca$x),aes(PC1,PC2,label=pheno.ct$Sample_Name))+
  geom_text()+
  theme_classic()
dev.off()


# 5k mv CpGs --------------------------------------------------------------
meth.sd <- apply(meth.ct,1,sd)
sorted <- order(meth.sd,decreasing = T)[1:5000]
meth.var <- meth.ct[sorted,]

meth.t <- t(meth.var)
pca <- prcomp(meth.t,center = T,scale. = T)
sum.pca <- summary(pca)
write.table(pca$x[,1:5],paste0(table.dir,Sys.Date(),"_pca_5k_mvCpG.txt"),sep="\t",quote=F,row.names = T)

pdf(paste0(plot.dir,Sys.Date(),"_mv_5k_cpg_pca.pdf"))
ggplot(as.data.frame(pca$x),aes(PC1,PC2,colour=pheno.ct$CellType))+
  geom_point(size=5)+
  theme_classic()+
  scale_colour_manual(values=colours)+
  xlab(paste0("PC1 (",round(sum.pca$importance[2,1]*100,2),"% Variance)"))+
  ylab(paste0("PC2 (",round(sum.pca$importance[2,2]*100,2),"% Variance)"))+
  theme(legend.position = "none") 
ggplot(as.data.frame(pca$x),aes(PC3,PC4,colour=pheno.ct$CellType))+
  geom_point(size=5)+
  theme_classic()+
  scale_colour_manual(values=colours)+
  xlab(paste0("PC3 (",round(sum.pca$importance[2,3]*100,2),"% Variance)"))+
  ylab(paste0("PC4 (",round(sum.pca$importance[2,4]*100,2),"% Variance)"))+
  theme(legend.position = "none") 
ggplot(as.data.frame(pca$x),aes(PC1,PC2,colour=pheno.ct$CellType,shape=factor(pheno.ct$Replicate)))+
  geom_point(size=5)+
  theme_classic()+
  scale_colour_manual(values=colours)
ggplot(as.data.frame(pca$x),aes(PC1,PC2,label=pheno.ct$Sample_Name))+
  geom_text()+
  theme_classic()
ggplot(as.data.frame(pca$x),aes(PC2,PC3,label=pheno.ct$Sample_Name))+
  geom_text()+
  theme_classic()
dev.off()


# Heatmap -----------------------------------------------------------------
annot.df <- pheno.ct
rownames(annot.df) <- annot.df$Sample_Name
annot.df <- annot.df[,c("Replicate","CellType")]
anno.colour <- list("CellType"=colours,
                    "Replicate"=c("gray90","gray60","gray20"))

pdf(paste0(plot.dir,Sys.Date(),"_mv_5k_heatmap.pdf"))
pheatmap::pheatmap(meth.var[,c(1:12,22:27,13:21)],
                   show_rownames = F,
                   cluster_cols = F,
                   annotation_col = annot.df,
                   annotation_colors = anno.colour,
                   #color = colorRampPalette(c("steelblue", "white","darkred"))(50),
                   treeheight_row = 0
                   )
pheatmap::pheatmap(meth.var,
                   show_rownames = F,
                   cluster_cols = T,
                   annotation_col = annot.df,
                   annotation_colors = anno.colour,
                   clustering_distance_cols = "manhattan",
                   treeheight_row = 0
)
dev.off()
tiff(paste0(plot.dir,Sys.Date(),"_mv_5k_heatmap.tiff"))
pheatmap::pheatmap(meth.var[,c(1:12,22:27,13:21)],
                   show_rownames = F,
                   cluster_cols = F,
                   annotation_col = annot.df,
                   annotation_colors = anno.colour,
                   #color = colorRampPalette(c("steelblue", "white","darkred"))(50),
                   treeheight_row = 0
)
dev.off()


# Correlation replicates --------------------------------------------------
cor(meth.ct[,"lsk_rep1_bm"],meth.ct[,"lsk_rep2_bm"])%>%round(digits = 3)
cor(meth.ct[,"lsk_rep1_bm"],meth.ct[,"lsk_rep3_bm"])%>%round(digits = 3)
cor.test(meth.ct[,"lsk_rep1_bm"],meth.ct[,"lsk_rep2_bm"])
cor.test(meth.ct[,"lsk_rep1_bm"],meth.ct[,"lsk_rep3_bm"])

pdf(paste0(plot.dir,Sys.Date(),"_lsk_replicates.pdf"))
ggplot(as.data.frame(meth.ct),aes(lsk_rep1_bm,lsk_rep2_bm))+
  stat_pointdensity(geom = GeomPointRast) + 
  scale_color_viridis()+
  theme_classic()+
  ylim(0,1)+
  xlim(0,1)+
  stat_cor()+
  ggtitle("LSK TWGBS")
ggplot(as.data.frame(meth.ct),aes(lsk_rep1_bm,lsk_rep3_bm))+
  stat_pointdensity(geom = GeomPointRast) + 
  scale_color_viridis()+
  theme_classic()+
  ylim(0,1)+
  xlim(0,1)+
  stat_cor()+
  ggtitle("LSK TWGBS")
dev.off()