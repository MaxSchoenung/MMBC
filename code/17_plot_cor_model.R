###################################################################
#
# Script: 17_plot_cor_model.R
# Plot the correlation between gene expression and methylation
# Date: 2021-03-15
# Author: Maximilian Schönung
#
###################################################################

# Load Libraries ----------------------------------------------------------
library(dplyr)
library(pheatmap)
library(biomaRt)
library(RnBeads)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(DESeq2)
library(ChIPpeakAnno)
library(ggplot2)
library(VennDiagram)
library(caret)
library(ggpointdensity)
library(viridis)
library(lemon)

# Set Directories ---------------------------------------------------------
setwd("/icgc/analysis/OE0565_projects/mmbc_hierarchy/")
data.dir <- "/icgc/analysis/OE0565_projects/mmbc_hierarchy/data/"
plot.dir <- "/icgc/analysis/OE0565_projects/mmbc_hierarchy/plots/"
table.dir <- "/icgc/analysis/OE0565_projects/mmbc_hierarchy/tables/"


# Define the colours ------------------------------------------------------
colours <- c("LSK"="gray50",
             "CMP"="darkorange",
             "GMP"="darkgoldenrod2",
             "MEP"="orangered4",
             "Monocyte"="gold1",
             "Neutrophil"="darksalmon",
             "B-cell"="navy",
             "CD8_t-cell"="lightsteelblue1",
             "CD4_t-cell"="lightskyblue")

colours2 <- c("LSK"="gray50",
             "CMP"="darkorange",
             "GMP"="darkgoldenrod2",
             "MEP"="orangered4",
             "Mono"="gold1",
             "Neutro"="darksalmon",
             "BCell"="navy",
             "CD8"="lightsteelblue1",
             "CD4"="lightskyblue")

cluster.colours <- c("4"="red4",
                     "1"="lightskyblue1",
                     "2"="dodgerblue4",
                     "9"="lightcyan1",
                     "3"="lightcyan3",
                     "8"="navy",
                     "5"="gold",
                     "6"="orange",
                     "7"="gray68")

# Load the methylation data -----------------------------------------------
rnb.full <- readRDS("/icgc/analysis/OE0565_projects/mmbc_hierarchy/data/2021-03-15_rnbeadsset_betas_scaling.RDS")
rnb <- remove.samples(rnb.full,28:36)
rm(rnb.full)

## access the meth and pheno data
pheno.rnb <- pheno(rnb)
pheno.rnb$Sample_Name
meth.rnb <- meth(rnb,row.names=T)
colnames(meth.rnb) <- pheno.rnb$Sample_Name

# Load the HSC dmps and novel CRE programs -------------------------------------------------------
dmps <- read.delim(paste0(table.dir,"2021-05-14_diff_meth_lsk_p_0.05_delt_0.2.txt"),stringsAsFactors = F)
programs <- read.delim(paste0(table.dir,"2021-07-30_new_cre_9_cluster.txt"),stringsAsFactors = F)

# Load the annotation -----------------------------------------------------
annot.full <- read.delim(paste0(table.dir,"gene-annos_primary_one-row.bed"),stringsAsFactors = F)
colnames(annot.full)[1] <- "Chromosome"
annot.full <- annot.full[annot.full$Chromosome%in%paste0("chr",1:19),]
rownames(annot.full) <- annot.full$name
annot <- annot.full[,1:3]
annot.ranges <- GenomicRanges::makeGRangesFromDataFrame(annot)

## annotate the dmps
dmps.annot <- cbind(annot[dmps$site,],dmps)
programs.annot <- cbind(annot[rownames(programs),],programs)
dmp.ranges <- makeGRangesFromDataFrame(dmps.annot,keep.extra.columns = T)
programs.ranges <- makeGRangesFromDataFrame(programs.annot,keep.extra.columns = T)

# Load the RNA-Seq data ---------------------------------------------------
dds <- readRDS(paste0(data.dir,"haemopedia_dds_mmbc.RDS"))
dds_myeloid <- readRDS(paste0(data.dir,"haemopedia_dds_mmbc_myeloid.RDS"))
dds <- DESeq(dds)
dds_myeloid <- DESeq(dds_myeloid)
resultsNames(dds)
resultsNames(dds_myeloid)
rld <- vst(dds,blind = F)
mat <- assay(rld)

## rename the results and perform lfcShrink
res.list <- list()
for(i in 2:length(resultsNames(dds))){
  message(paste0("Processing sample:",i))  
  res.list[[i]] <- lfcShrink(dds,coef=i,type="ashr",quiet=T,format="DataFrame")%>%as.data.frame
}
res.list[[10]] <- lfcShrink(dds_myeloid,coef=7,type="ashr",quiet=T,format="DataFrame")%>%as.data.frame
res.list[[11]] <- lfcShrink(dds_myeloid,coef=2,type="ashr",quiet=T,format="DataFrame")%>%as.data.frame
res.list <- res.list[2:11]
names(res.list) <- c(resultsNames(dds)[2:9],resultsNames(dds_myeloid)[7],resultsNames(dds_myeloid)[2])
res.filt <- lapply(res.list,function(x)x[x$padj<0.01&x$log2FoldChange>2,])
lapply(res.filt,dim)

# Load the ENSEMBL Mart ----------------------------------------------------------------
mm10 <- useMart("ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl",host="http://nov2020.archive.ensembl.org") #biomart for annotation

# Mean Methylation per grop -----------------------------------------------
mean.meth <- data.frame(
  "LSK"=rowMeans(meth.rnb[,1:3]),
  "MEP"=rowMeans(meth.rnb[,4:6]),
  "CMP"=rowMeans(meth.rnb[,7:9]),
  "GMP"=rowMeans(meth.rnb[,10:12]),
  "CD4"=rowMeans(meth.rnb[,13:15]),
  "CD8"=rowMeans(meth.rnb[,16:18]),
  "BCell"=rowMeans(meth.rnb[,19:21]),
  "Mono"=rowMeans(meth.rnb[,22:24]),
  "Neutro"=rowMeans(meth.rnb[,25:27])
)
mean.exp <- data.frame(
  "LSK"=rowMeans(mat[,15:16]),
  "MEP"=rowMeans(mat[,9:10]),
  "CMP"=rowMeans(mat[,23:24]),
  "GMP"=rowMeans(mat[,25:26]),
  "CD4"=rowMeans(mat[,27:28]),
  "CD8"=rowMeans(mat[,29:30]),
  "BCell"=rowMeans(mat[,1:8]),
  "Mono"=rowMeans(mat[,11:14]),
  "Neutro"=rowMeans(mat[,17:22])
)


# Load the correlation model ----------------------------------------------
cor.distinct <- readRDS(paste0(data.dir,"2021-09-17_novelCRE_9programs_cpg_gene_corr.RDS"))
cor.sig <- cor.distinct[cor.distinct$p_adj<0.01&abs(cor.distinct$slope)>4,]


## add the methylation and expression data
merged.cor <- merge(cor.sig,mean.meth,by.x="cpg",by.y=0)
merged.cor.melt <- reshape2::melt(merged.cor,id.vars=c("cpg","gene","cor_value","p_value","p_adj","slope","dist_TSS","chr","TSS_start",
                                                       "cpgs_start","ind","gene_id"))
exp.melt <- reshape2::melt(as.matrix(mean.exp))
merged.cor.melt.full <- merge(merged.cor.melt,exp.melt,by.x=c("gene","variable"),by.y=c("Var1","Var2"),all.x=T,all.y=F)
colnames(merged.cor.melt.full)
colnames(merged.cor.melt.full)[14:15] <- c("methylation","expression")
head(merged.cor.melt.full)

# Check the correlation of gene expression and methylation for each cell type specific DMP set
clust.vec <- c(4,1,2,9,9,3,3,3,8,5,5,6)
pat.vec <- c("MEP","CD4","CD4","T.cell","T.cell","CD4","CD8","B","B","Myeloid","Myeloid","Neutro")
ct.vec <- c("MEP","CD4","CD4","CD4","CD8","CD4","CD8","BCell","BCell","Mono","Neutro","Neutro")

for(i in 1:12){
  message(i)
  cluster.select <- clust.vec[i]
  exp.select <- grep(pat.vec[i],names(res.list))
  ct.select <- ct.vec[i]


plot.set <- merged.cor.melt.full[which(merged.cor.melt.full$cpg%in%rownames(programs[which(programs$cluster==cluster.select),,drop=F])),]
plot.set <- plot.set[abs(plot.set$slope)>4,]
plot.set2 <- merge(plot.set,res.list[[exp.select]][,c("padj","log2FoldChange"),drop=F],by.x="gene",by.y="row.names",all.x=T,all.y=F)
tx.id2 <- getBM(attributes=c("external_gene_name","ensembl_gene_id"),filters="ensembl_gene_id",values=plot.set2$gene,mart=mm10)
plot.set3 <- merge(plot.set2,tx.id2,by.x="gene",by.y="ensembl_gene_id",all.x=T,all.y=F)
plot.set3$var2 <- factor(plot.set3$variable==ct.select,labels = c("Others",ct.select))
  
pdf(paste0(plot.dir,Sys.Date(),"_cor_cluster_",clust.vec[i],"_",ct.vec[i],".pdf"),height = 4)
print(ggplot(plot.set3,aes(expression,methylation,label=gene,colour=log2FoldChange))+
  geom_point()+
  theme_classic()+
  facet_wrap(~var2)+
  theme(legend.position="bottom")+
  xlab("Expression")+
  ylab("Beta-Value")+
  scale_colour_gradient(low = "white", high = "black"))
print(ggplot(plot.set3,aes(expression,methylation,label=gene,colour=as.factor(variable)))+
  geom_point(size=1)+
  theme_classic()+
  facet_wrap(~var2)+
  theme(legend.position="bottom")+
  xlab("Expression")+
  ylab("Beta-Value")+
  scale_colour_manual(values = colours2)
  )
print(ggplot(plot.set3,aes(expression,methylation,label=external_gene_name,colour=log2FoldChange))+
  geom_text(size=1)+
  theme_classic()+
  facet_wrap(~var2)+
  theme(legend.position="bottom")+
  xlab("Expression")+
  ylab("Beta-Value")+
  scale_colour_gradient(low = "white", high = "black"))
dev.off()
}

# Check the correlation of gene expression and methylation for each cell type specific DMP set
cluster.select <- 4
exp.select <- grep("MEP",names(res.list))
table(plot.set3$variable)
ct.select <- "MEP"

plot.set <- merged.cor.melt.full[which(merged.cor.melt.full$cpg%in%rownames(programs[which(programs$cluster==cluster.select),,drop=F])),]
plot.set <- plot.set[abs(plot.set$slope)>4,]
plot.set2 <- merge(plot.set,res.list[[exp.select]][,c("padj","log2FoldChange"),drop=F],by.x="gene",by.y="row.names",all.x=T,all.y=F)
tx.id2 <- getBM(attributes=c("external_gene_name","ensembl_gene_id"),filters="ensembl_gene_id",values=plot.set2$gene,mart=mm10)
plot.set3 <- merge(plot.set2,tx.id2,by.x="gene",by.y="ensembl_gene_id",all.x=T,all.y=F)
plot.set3$var2 <- factor(plot.set3$variable==ct.select,labels = c("Others",ct.select))
plot.set4 <- plot.set3%>%
  group_by(var2,ind)%>%
  dplyr::summarize(meth_mean=mean(methylation),
                   meth_sd=sd(methylation),
                   exp_mean=mean(expression),
                   exp_sd=sd(expression),
                   log2FC=mean(log2FoldChange),
                   gene=unique(gene_id))

pdf(paste0(plot.dir,Sys.Date(),"_cor_cluster4_MEP_plotting.pdf"),height = 4)
ggplot(plot.set4,aes(exp_mean,meth_mean,colour=log2FC))+
  geom_point(aes(size=log2FC))+
  geom_errorbar(aes(ymin=meth_mean-meth_sd, ymax=meth_mean+meth_sd))+
  geom_errorbarh(aes(xmin=exp_mean-exp_sd,xmax=exp_mean+exp_sd))+
  ggrepel::geom_label_repel(data=plot.set4[plot.set4$log2FC>7,],aes(label=gene),size=2,
                            min.segment.length = unit(0, 'lines'),
                            segment.linetype = 2,
                            max.overlaps=Inf,
                            nudge_x = 3,
                            nudge_y = .1)+
  theme_classic()+
  facet_wrap(~var2)+
  theme(legend.position="bottom")+
  xlab("Expression")+
  ylab("Beta-Value")+
  scale_size(range=c(.1,3))+
  scale_color_gradient2(low=muted("blue"),mid = "white",high=cluster.colours["4"],midpoint = 0)
dev.off()

# Check the correlation of gene expression and methylation for each cell type specific DMP set
cluster.select <- 8
exp.select <- grep("B",names(res.list))
table(plot.set3$variable)
ct.select <- "BCell"

plot.set <- merged.cor.melt.full[which(merged.cor.melt.full$cpg%in%rownames(programs[which(programs$cluster==cluster.select),,drop=F])),]
plot.set <- plot.set[abs(plot.set$slope)>4,]
plot.set2 <- merge(plot.set,res.list[[exp.select]][,c("padj","log2FoldChange"),drop=F],by.x="gene",by.y="row.names",all.x=T,all.y=F)
tx.id2 <- getBM(attributes=c("external_gene_name","ensembl_gene_id"),filters="ensembl_gene_id",values=plot.set2$gene,mart=mm10)
plot.set3 <- merge(plot.set2,tx.id2,by.x="gene",by.y="ensembl_gene_id",all.x=T,all.y=F)
plot.set3$var2 <- factor(plot.set3$variable==ct.select,labels = c("Others",ct.select))
plot.set4 <- plot.set3%>%
  group_by(var2,ind)%>%
  dplyr::summarize(meth_mean=mean(methylation),
                   meth_sd=sd(methylation),
                   exp_mean=mean(expression),
                   exp_sd=sd(expression),
                   log2FC=mean(log2FoldChange),
                   gene=unique(gene_id))

pdf(paste0(plot.dir,Sys.Date(),"_cor_cluster8_BCell_plotting.pdf"),height = 4)
ggplot(plot.set4,aes(exp_mean,meth_mean,colour=log2FC))+
  geom_point(aes(size=log2FC))+
  geom_errorbar(aes(ymin=meth_mean-meth_sd, ymax=meth_mean+meth_sd))+
  geom_errorbarh(aes(xmin=exp_mean-exp_sd,xmax=exp_mean+exp_sd))+
  ggrepel::geom_label_repel(data=plot.set4[plot.set4$log2FC>11,],aes(label=gene),size=2,
                            min.segment.length = unit(0, 'lines'),
                            segment.linetype = 2,
                            max.overlaps=Inf,
                            nudge_x = 3,
                            nudge_y = .1)+
  theme_classic()+
  facet_wrap(~var2)+
  theme(legend.position="bottom")+
  xlab("Expression")+
  ylab("Beta-Value")+
  scale_size(range=c(.1,3))+
  #scale_colour_gradient(low = "white", high = cluster.colours["8"])
  scale_color_gradient2(low=muted("blue"),mid = "white",high=cluster.colours["8"],midpoint = 0)
dev.off()


# Check the correlation of gene expression and methylation for each cell type specific DMP set
cluster.select <- 6
exp.select <- grep("Neutro",names(res.list))
table(plot.set3$variable)
ct.select <- "Neutro"

plot.set <- merged.cor.melt.full[which(merged.cor.melt.full$cpg%in%rownames(programs[which(programs$cluster==cluster.select),,drop=F])),]
plot.set <- plot.set[abs(plot.set$slope)>4,]
plot.set2 <- merge(plot.set,res.list[[exp.select]][,c("padj","log2FoldChange"),drop=F],by.x="gene",by.y="row.names",all.x=T,all.y=F)
tx.id2 <- getBM(attributes=c("external_gene_name","ensembl_gene_id"),filters="ensembl_gene_id",values=plot.set2$gene,mart=mm10)
plot.set3 <- merge(plot.set2,tx.id2,by.x="gene",by.y="ensembl_gene_id",all.x=T,all.y=F)
plot.set3$var2 <- factor(plot.set3$variable==ct.select,labels = c("Others",ct.select))
plot.set4 <- plot.set3%>%
  group_by(var2,ind)%>%
  dplyr::summarize(meth_mean=mean(methylation),
                   meth_sd=sd(methylation),
                   exp_mean=mean(expression),
                   exp_sd=sd(expression),
                   log2FC=mean(log2FoldChange),
                   gene=unique(gene_id))

pdf(paste0(plot.dir,Sys.Date(),"_cor_cluster6_Neutro_plotting.pdf"),height = 4)
ggplot(plot.set4,aes(exp_mean,meth_mean,colour=log2FC))+
  geom_point(aes(size=log2FC))+
  geom_errorbar(aes(ymin=meth_mean-meth_sd, ymax=meth_mean+meth_sd))+
  geom_errorbarh(aes(xmin=exp_mean-exp_sd,xmax=exp_mean+exp_sd))+
  ggrepel::geom_label_repel(data=plot.set4[plot.set4$log2FC>12.5,],aes(label=gene),size=2,
                            min.segment.length = unit(0, 'lines'),
                            segment.linetype = 2,
                            max.overlaps=Inf,
                            nudge_x = 5,
                            nudge_y = .1)+
  theme_classic()+
  facet_wrap(~var2)+
  theme(legend.position="bottom")+
  xlab("Expression")+
  ylab("Beta-Value")+
  scale_size(range=c(.1,3))+
  #scale_colour_gradient(low = "white", high = cluster.colours["8"])
  scale_color_gradient2(low=muted("blue"),mid = "white",high=cluster.colours["6"],midpoint = 0)
dev.off()

# Check the correlation of gene expression and methylation for each cell type specific DMP set
cluster.select <- 5
exp.select <- grep("Myeloid",names(res.list))
table(plot.set3$variable)
ct.select <- "Neutro"

plot.set <- merged.cor.melt.full[which(merged.cor.melt.full$cpg%in%rownames(programs[which(programs$cluster==cluster.select),,drop=F])),]
plot.set <- plot.set[abs(plot.set$slope)>4,]
plot.set2 <- merge(plot.set,res.list[[exp.select]][,c("padj","log2FoldChange"),drop=F],by.x="gene",by.y="row.names",all.x=T,all.y=F)
tx.id2 <- getBM(attributes=c("external_gene_name","ensembl_gene_id"),filters="ensembl_gene_id",values=plot.set2$gene,mart=mm10)
plot.set3 <- merge(plot.set2,tx.id2,by.x="gene",by.y="ensembl_gene_id",all.x=T,all.y=F)
plot.set3$var2 <- factor(plot.set3$variable==ct.select,labels = c("Others",ct.select))
levels(plot.set3$var2) <- c("Others","Neutro","GMP","Mono")
plot.set3$var2[plot.set3$variable=="GMP"] <- "GMP"
plot.set3$var2[plot.set3$variable=="Mono"] <- "Mono"
plot.set3$var2 <- factor(plot.set3$var2,levels=c("Others","GMP","Mono","Neutro"))
plot.set4 <- plot.set3%>%
  group_by(var2,ind)%>%
  dplyr::summarize(meth_mean=mean(methylation),
                   meth_sd=sd(methylation),
                   exp_mean=mean(expression),
                   exp_sd=sd(expression),
                   log2FC=mean(log2FoldChange),
                   gene=unique(gene_id))

pdf(paste0(plot.dir,Sys.Date(),"_cor_cluster_5_Myeloid_plotting.pdf"),height = 6,width=8)
ggplot(plot.set4,aes(exp_mean,meth_mean,colour=log2FC))+
  geom_point(aes(size=log2FC))+
  geom_errorbar(aes(ymin=meth_mean-meth_sd, ymax=meth_mean+meth_sd))+
  geom_errorbarh(aes(xmin=exp_mean-exp_sd,xmax=exp_mean+exp_sd))+
  ggrepel::geom_label_repel(data=plot.set4[plot.set4$log2FC>8,],aes(label=gene),size=2,
                            min.segment.length = unit(0, 'lines'),
                            segment.linetype = 2,
                            max.overlaps=Inf,
                            #nudge_x = 2,
                            nudge_y = .1)+
  theme_classic()+
  facet_wrap(~var2,ncol=2)+
  theme(legend.position="bottom")+
  xlab("Expression")+
  ylab("Beta-Value")+
  scale_size(range=c(.1,3))+
  #scale_colour_gradient(low = "white", high = cluster.colours["5"])
  scale_color_gradient2(low=muted("blue"),mid = "white",high="gold",midpoint = 0)
dev.off()

# Check the correlation of gene expression and methylation for each cell type specific DMP set
cluster.select <- 9
exp.select <- grep("T.cell",names(res.list))
table(plot.set3$variable)
ct.select <- "CD8"

plot.set <- merged.cor.melt.full[which(merged.cor.melt.full$cpg%in%rownames(programs[which(programs$cluster==cluster.select),,drop=F])),]
plot.set <- plot.set[abs(plot.set$slope)>4,]
plot.set2 <- merge(plot.set,res.list[[exp.select]][,c("padj","log2FoldChange"),drop=F],by.x="gene",by.y="row.names",all.x=T,all.y=F)
tx.id2 <- getBM(attributes=c("external_gene_name","ensembl_gene_id"),filters="ensembl_gene_id",values=plot.set2$gene,mart=mm10)
plot.set3 <- merge(plot.set2,tx.id2,by.x="gene",by.y="ensembl_gene_id",all.x=T,all.y=F)
plot.set3$var2 <- factor(plot.set3$variable==ct.select,labels = c("Others",ct.select))
levels(plot.set3$var2) <- c("Others","CD8","CD4")
plot.set3$var2[plot.set3$variable=="CD4"] <- "CD4"
plot.set3$var2 <- factor(plot.set3$var2,levels=c("Others","CD4","CD8"))
plot.set4 <- plot.set3%>%
  group_by(var2,ind)%>%
  dplyr::summarize(meth_mean=mean(methylation),
                   meth_sd=sd(methylation),
                   exp_mean=mean(expression),
                   exp_sd=sd(expression),
                   log2FC=mean(log2FoldChange),
                   gene=unique(gene_id))

pdf(paste0(plot.dir,Sys.Date(),"_cor_cluster_9_T-cell_plotting.pdf"),height = 4,width=8)
ggplot(plot.set4,aes(exp_mean,meth_mean,colour=log2FC))+
  geom_point(aes(size=log2FC))+
  geom_errorbar(aes(ymin=meth_mean-meth_sd, ymax=meth_mean+meth_sd))+
  geom_errorbarh(aes(xmin=exp_mean-exp_sd,xmax=exp_mean+exp_sd))+
  ggrepel::geom_label_repel(data=plot.set4[plot.set4$log2FC>12,],aes(label=gene),size=2,
                            min.segment.length = unit(0, 'lines'),
                            segment.linetype = 2,
                            max.overlaps=Inf,
                            nudge_x = 2)+
  theme_classic()+
  facet_wrap(~var2,ncol=3)+
  theme(legend.position="bottom")+
  xlab("Expression")+
  ylab("Beta-Value")+
  scale_size(range=c(.1,3))+
  #scale_colour_gradient(low = "white", high = cluster.colours["5"])
  scale_color_gradient2(low=muted("blue"),mid = "white",high="forestgreen",midpoint = 0)
dev.off()


# Cluster with more than 10 associations ----------------------------------
cor.sig.programs <- merge(cor.sig,programs,by.x="cpg",by.y=0,all.x=T,all.y=F)
head(cor.sig.programs)
table(cor.sig.programs$clusters)



# Generate a large plot CT x cluster --------------------------------------
merged.programs.full <- merge(merged.cor.melt.full,programs,by.x="cpg",by.y=0)
merged.programs.full$clusters <- factor(merged.programs.full$clusters,levels=c(4,1,2,9,3,8,5,6,7))
pdf(paste0(plot.dir,Sys.Date(),"_cor_clusters_overview.pdf"),height = 12,width=8)
ggplot(merged.programs.full,aes(expression,methylation,colour=variable))+
  geom_point(alpha=.4)+
  theme_classic()+
  facet_grid(clusters~variable)+
  theme(legend.position="bottom")+
  xlab("Expression")+
  ylab("Beta-Value")+
  scale_colour_manual(values = colours2)+
  coord_capped_cart(bottom='both', left='both')+
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf,size=1)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=1)
dev.off()
