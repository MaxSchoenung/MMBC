###################################################################
#
# Script: 15_cluster_cres_v3.R
# Clustering of CREs
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
library(rGREAT)
library(grid)
library(RColorBrewer)

# Set Directories ---------------------------------------------------------
setwd("/icgc/analysis/OE0565_projects/mmbc_hierarchy/")
data.dir <- "/icgc/analysis/OE0565_projects/mmbc_hierarchy/data/"
save.dir <- "/icgc/analysis/OE0565_projects/mmbc_hierarchy/RDS/differential/"
plot.dir <- "/icgc/analysis/OE0565_projects/mmbc_hierarchy/plots/"
table.dir <- "/icgc/analysis/OE0565_projects/mmbc_hierarchy/tables/"


# Load the dataset --------------------------------------------------------
rnb.full <- readRDS("/icgc/analysis/OE0565_projects/mmbc_hierarchy/data/2021-03-15_rnbeadsset_betas_scaling.RDS")
rnb <- remove.samples(rnb.full,28:36)

## access the meth and pheno data
pheno.rnb <- pheno(rnb)
pheno.rnb$Sample_Name
meth.rnb <- meth(rnb,row.names=T)
colnames(meth.rnb) <- pheno.rnb$Sample_Name

# Load the annotation from stephen ----------------------------------------
annot <- read.delim(paste0(table.dir,"2022-05-25_MMBC_annotation_modified.bed"),stringsAsFactors = F)
rownames(annot) <- annot$IlmnID

# Load the HSC dmps -------------------------------------------------------
dmps <- read.delim(paste0(table.dir,"2021-05-14_diff_meth_lsk_p_0.05_delt_0.2.txt"),stringsAsFactors = F)
dmps.unique <- dmps[!duplicated(dmps$site),]

# Load the CREs -------------------------------------------------------
mdCRE <- readRDS(paste0(data.dir,"2022-03-08_overlapping_CRE.RDS"))
chromCRE <- readRDS(paste0(data.dir,"2022-03-08_previous_CRE.RDS"))
novel_mdCRE <- readRDS(paste0(data.dir,"2022-03-08_novel_mdCRE.RDS"))
cre.list <- list("univ_mdCRE"=mdCRE,"chromCRE"=chromCRE,"novel_mdCRE"=novel_mdCRE)

# Clustering with HSC DMPs ------------------------------------------------
cre.meth <- lapply(cre.list,function(x){meth.rnb[which(rownames(meth.rnb)%in%x),]})
lapply(cre.meth,dim)
lapply(cre.list,length)

## mean methylation
mean.meth_fun <- function(x){
  data.frame(
  "LSK"=rowMeans(x[,1:3]),
  "MEP"=rowMeans(x[,4:6]),
  "CMP"=rowMeans(x[,7:9]),
  "GMP"=rowMeans(x[,10:12]),
  "CD4"=rowMeans(x[,13:15]),
  "CD8"=rowMeans(x[,16:18]),
  "BCell"=rowMeans(x[,19:21]),
  "Mono"=rowMeans(x[,22:24]),
  "Neutro"=rowMeans(x[,25:27])
)}

cre.mean <- lapply(cre.meth,mean.meth_fun)

#focus on 1k random sites as the sets are too big
set.seed(1234)
cre.rand <- lapply(cre.mean,function(x){x[sample(rownames(x),1000,replace = F),]})
  
## calculate the z-score
fun_zscore <- function(x){(x - mean(x)) / sd(x)}
cre.z <- lapply(cre.rand,function(x){t(apply(x,1,fun_zscore))})
cre.dist <- lapply(cre.z,function(x)dist(x, method="euclidean"))
cre.clust <- lapply(cre.dist,function(x)hclust(x,method = "ward.D2"))

# Plot as a heatmap -------------------------------------------------------
pdf(paste0(plot.dir,Sys.Date(),"_cluster_CRE_sets.pdf"))
for(i in 1:length(cre.clust)){
print(pheatmap(cre.z[[i]],
         breaks=seq(-3,3,by = 0.1),
         color = colorRampPalette(c("steelblue", "white","darkred"))(length(seq(-3,3,by = 0.1))),
         show_rownames = F,
         cluster_rows=cre.clust[[i]],
         cluster_cols = F,
         main = paste0("Z-Score Clustering, 5k mvCpG, Euc. Dist., Ward.D2: ",names(cre.clust)[i])))
print(pheatmap(cre.rand[[i]],
         breaks=seq(0,1,by = 0.01),       
         show_rownames = F,
         cluster_rows=cre.clust[[i]],
         cluster_cols = F,
         main = paste0("Z-Score Clustering, 5k mvCpG, Euc. Dist., Ward.D2: ",names(cre.clust)[i])))
}
dev.off()


# Export clusters for cmdCREs ---------------------------------------------
novel_md.cre.z <- t(apply(cre.mean[[3]],1,fun_zscore))
novel_md.dist <- dist(novel_md.cre.z, method="euclidean")
novel_md.clust <- hclust(novel_md.dist,method = "ward.D2")

## row anno
clusts <- cutree(novel_md.clust,9)
anno.row <- data.frame("clusters"=factor(clusts),row.names = names(clusts))

## colouring annotation
row.colours <- list("clusters"=
                      c("4"="red4",
                        "1"="lightskyblue1",
                        "2"="dodgerblue4",
                        "9"="lightcyan1",
                        "3"="lightcyan3",
                        "8"="navy",
                        "5"="gold",
                        "6"="orange",
                        "7"="gray68"))
                  

## plot the clustering
pdf(paste0(plot.dir,Sys.Date(),"_cluster_cmdCREs.pdf"))
pheatmap(novel_md.cre.z[,c(1:4,8,9,7,6,5)],
         color = colorRampPalette(c("steelblue", "white","darkred"))(100),
         show_rownames = F,
         cluster_rows=novel_md.clust,
         cluster_cols = F,
         annotation_row = anno.row,
         annotation_colors = row.colours,
         main = "Euclidean Distance, WardD2, z-score Clustering")
pheatmap(cre.mean[[3]][,c(1:4,8,9,7,6,5)],
         show_rownames = F,
         cluster_rows=novel_md.clust,
         cluster_cols = F,
         annotation_row = anno.row,
         annotation_colors = row.colours,
         main = "Euclidean Distance, WardD2, z-score Clustering")
dev.off()

# Cluster within cluster --------------------------------------------------
set.seed(1234)
meth.plot <- novel_md.cre.z
for(i in c(4,1,2,9,3,8,5:7)){
  ids <- rownames(anno.row)[anno.row$clusters==i]
  inds <- which(rownames(meth.plot)%in%ids)
  inds <- sample(inds,size = 350,replace = F) #as smalles cluster has 371 sites (cluster 3)
  meth.sub <- meth.plot[inds,]
  meth.z.sub <- t(apply(meth.sub,1,fun_zscore))
  euc.dist.z.sub <- dist(meth.z.sub, method="euclidean")
  clustering.euc.z.sub <- hclust(euc.dist.z.sub,method = "ward.D2")
  if(i==4){
    out.meth <- meth.sub[clustering.euc.z.sub$order,]
  }else{
    out.meth <- rbind(out.meth,meth.sub[clustering.euc.z.sub$order,])
  }
}

pdf(paste0(plot.dir,Sys.Date(),"_cmdCREs_cluster_inCluster.pdf"))
pheatmap(out.meth,
         color = colorRampPalette(c(muted("blue"), "white",muted("red")))(100),
         show_rownames = F,
         cluster_rows=F,
         cluster_cols = F,
         annotation_row = anno.row,
         gaps_row=seq(350,nrow(out.meth),by = 350),
         annotation_colors = row.colours,
         main = "Euclidean Distance, WardD2, z-score Clustering")
pheatmap(cre.mean[[3]][rownames(out.meth),colnames(out.meth)],
         show_rownames = F,
         cluster_rows=F,
         cluster_cols = F,
         annotation_row = anno.row,
         gaps_row=seq(350,nrow(out.meth),by = 350),
         annotation_colors = row.colours,
         main = "Euclidean Distance, WardD2, z-score Clustering")
dev.off()

# Export the clusters -----------------------------------------------------
write.table(anno.row,paste0(table.dir,Sys.Date(),"_new_cre_9_cluster.txt"),sep="\t",quote=F,row.names = T)
test.df <- read.delim((paste0(table.dir,"2021-09-16_new_cre_9_cluster.txt")))
table(rownames(anno.row)==rownames(test.df))
table(anno.row$clusters==test.df$clusters)
export.df <- data.frame(annot[rownames(anno.row),c("Chromosome","Start","End","IlmnID")],
                        "cluster"=anno.row$clusters,
                        stringsAsFactors = F)
write.table(anno.row,paste0(table.dir,Sys.Date(),"_new_cre_9_cluster_annotated.txt"),sep="\t",quote=F,row.names = F)
